/*
 *  DynFullLaplacianInverseSolver.cpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <harmonicResistance/DynFullLaplacianInverseSolver.hpp>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {
DynFullLaplacianInverseSolver::DynFullLaplacianInverseSolver(const Graph &G)
    : DynLaplacianInverseSolver(G), lpinv(G.numberOfNodes()), decomp(G) {}

void DynFullLaplacianInverseSolver::run() {
    decomp.run();
    lpinv = computeLpinv(G, decomp);
    hasRun = true;
}

double DynFullLaplacianInverseSolver::totalResistanceDifference(const GraphEvent &ev) const {
    assureFinished();

    const node i = ev.u;
    const node j = ev.v;
    const auto col_i = lpinv.row(i); // lpinv is symmetric, Take row instead of column, because
                                     // DenseMatrix is stored row-major
    const auto col_j = lpinv.row(j);
    const double R_ij = col_i[i] + col_j[j] - 2 * col_i[j];
    double w;
    if (ev.type == GraphEvent::EDGE_ADDITION)
        w = 1.0 / (1.0 + R_ij);
    else if (ev.type == GraphEvent::EDGE_REMOVAL)
        w = 1.0 / (1.0 - R_ij);
    else
        throw std::logic_error(
            "Trace difference cannot be computed for events other than edge addition or deletion!");
    const auto norm = (col_i - col_j).length();
    return norm * norm * w * G.numberOfNodes();
}

double computeTotalHarmonicResistance(const DenseMatrix &lpinv,
                                      const ComponentDecomposition &decomp) {
    // iterate all pairs of each component
    double totalHarmonicResistance = 0;
    for (auto &component : decomp.getComponents()) {
        for (auto u : component) {
            for (auto v : component) {
                if (u < v) {
                    const double R_ij = lpinv(u, u) + lpinv(v, v) - 2 * lpinv(u, v);
                    totalHarmonicResistance += 1.0 / R_ij;
                }
            }
        }
    }
    return totalHarmonicResistance;
}

double
DynFullLaplacianInverseSolver::totalHarmonicResistanceDifference(const GraphEvent &ev) const {
    assureFinished();

    const node i = ev.u;
    const node j = ev.v;
    const auto col_i = lpinv.row(i); // lpinv is symmetric, Take row instead of column, because
                                     // DenseMatrix is stored row-major
    const auto col_j = lpinv.row(j);

    // cache the total harmonic resistance of the base graph for reuse in subsequent calls of this
    // function
    if (!totalHarmonicResistance) {
        totalHarmonicResistance = computeTotalHarmonicResistance(lpinv, decomp);
    }

    // compute new total harmonic resistance
    // two cases: removed edge is bridge (resistance of incident node pair is 1) or it is not a
    // bridge. handle both cases differently

    // compute resistance of the removed edge to differentiate the cases
    const double R_ij = lpinv(i, i) + lpinv(j, j) - 2 * lpinv(i, j);

    // resulting new resistance
    double newTotalHarmonicResistance = std::numeric_limits<double>::signaling_NaN();

    // handle both cases:
    if ((1 - R_ij) < 1e-6) { // if close to 1, it is a bridge
        INFO("[DynFullLaplacianInverseSolver::totalHarmonicResistanceDifference] identified (", i,
             ", ", j, ") as bridge.");
        auto newG = G;
        newG.removeEdge(ev.u, ev.v);
        // the connected components changed - we know exactly how, but there is no dynamic connected
        // components (TODO)
        auto newDecomp = ParallelConnectedComponents(newG);
        newDecomp.run();
        auto newLpinv = computeLpinv(newG, newDecomp);
        newTotalHarmonicResistance = computeTotalHarmonicResistance(newLpinv, newDecomp);
    } else {
        // otherwise apply Sherman-Morrison-Formula to compute updated lpinv
        auto newLpinv =
            lpinv + Vector::outerProduct<DenseMatrix>(col_i - col_j, col_i - col_j) / (1.0 - R_ij);
        // we know that the connected components stay the same in this case, so there is no need to
        // recompute
        newTotalHarmonicResistance = computeTotalHarmonicResistance(newLpinv, decomp);
    }

    if (!std::isfinite(newTotalHarmonicResistance)) {
        WARN("total harmonic resistance is infinite! event is ", ev.u, ", ", ev.v);
    }
    // return difference
    return std::abs(totalHarmonicResistance.value() - newTotalHarmonicResistance);
}

DenseMatrix
DynFullLaplacianInverseSolver::computeLpinv(const Graph &G,
                                            const ComponentDecomposition &decomp) const {
    auto lap = CSRMatrix::laplacianMatrix(G);
    const count n = lap.numberOfColumns();
    Lamg<CSRMatrix> lamg;
    if (decomp.numberOfComponents() == 1)
        lamg.setupConnected(lap);
    else
        lamg.setup(lap);

    DenseMatrix result(n);

    // for each component, build the rhsbase and add to map of bases
    std::map<count, Vector> rhsbases;
    for (auto &component : decomp.getComponents()) {
        Vector base(n);
        for (node v : component)
            base[v] = 1;
        base *= -1.0 / component.size();
        rhsbases.insert({decomp.componentOfNode(component[0]), base});
    }

    const count maxThreads = static_cast<count>(omp_get_max_threads());

    // Solution vectors: one per thread
    std::vector<Vector> solutions(maxThreads, Vector(n));

    // Right hand side vectors: one per thread
    std::vector<Vector> rhss(maxThreads, Vector(n));

    const count iters = (n % maxThreads == 0) ? n / maxThreads : n / maxThreads + 1;
    for (count i = 0; i < iters; ++i) {
        // Index of the next vertex to process
        const index base = i * maxThreads;

#pragma omp parallel
        {
            // Each thread solves a linear system from `base` to `base + #threads - 1`
            // const index thread = 0;
            const index thread = omp_get_thread_num();
            const node v = base + thread;
            if (v < n) {
                // Reset solution and rhs vector of the current thread
                solutions[thread].fill(0.0);

                // Set up system to compute the diagonal column
                rhss[thread] = rhsbases[decomp.componentOfNode(v)];
                rhss[thread][v] += 1.;
            }
        }

        if (base + maxThreads >= n) {
            // Last iteration: some threads cannot be used.
            // Resize rhss and solutions to the number of vertices left to be processed.
            rhss.resize(n - base);
            solutions.resize(rhss.size());
        }

        // parallelSolve is broken currently
        lamg.parallelSolve(rhss, solutions);
        // for (int i = 0; i < rhss.size(); ++i)
        //     lamg.solve(rhss[i], solutions[i]);

        // TODO: change this to have the row access in the outer loop ?
        // Store the results
        for (index idx = 0; idx < solutions.size(); ++idx) {
            const node v = base + idx;
            if (v < n) {
                for (index row = 0; row < n; ++row) {
                    result.setValue(row, v, solutions[idx][row]);
                }
            } else
                break;
        }
    }
    return result;
}

DenseMatrix DynFullLaplacianInverseSolver::updateShermanMorrison(const GraphEvent &ev) const {
    DenseMatrix result(G.numberOfNodes());
    if (ev.type == GraphEvent::EDGE_ADDITION)
        assert(G.hasEdge(ev.u, ev.v));
    if (ev.type == GraphEvent::EDGE_REMOVAL)
        assert(!G.hasEdge(ev.u, ev.v));

    const auto i = ev.u;
    const auto j = ev.v;
    const double R_ij = lpinv(i, i) + lpinv(j, j) - 2 * lpinv(i, j);

    double w_negative;
    if (ev.type == GraphEvent::EDGE_ADDITION)
        w_negative = -1.0 / (1.0 + R_ij);
    else if (ev.type == GraphEvent::EDGE_REMOVAL)
        w_negative = -1.0 / (1.0 - R_ij);
    else
        throw std::logic_error("update does not support events other than "
                               "edge addition or deletion!");
    const auto v = lpinv.row(i) - lpinv.row(j);

    const auto n = lpinv.numberOfRows();
#pragma omp parallel for
    for (index i = 0; i < n; i++) {
        const auto updateVec = v[i] * w_negative * v;
        for (index j = 0; j < n; j++) {
            result.setValue(j, i, lpinv(j, i) + updateVec[j]);
        }
    }
    return result;
}

double DynFullLaplacianInverseSolver::totalForestDistanceDifference(const GraphEvent &ev) const {
    assureFinished();

    const node i = ev.u;
    const node j = ev.v;
    const auto col_i = lpinv.row(i); // lpinv is symmetric, Take row instead of column, because
                                     // DenseMatrix is stored row-major
    const auto col_j = lpinv.row(j);
    const double R_ij = col_i[i] + col_j[j] - 2 * col_i[j];
    double w;
    if (ev.type == GraphEvent::EDGE_ADDITION)
        w = 1.0 / (1.0 + R_ij);
    else if (ev.type == GraphEvent::EDGE_REMOVAL)
        w = 1.0 / (1.0 - R_ij);
    else
        throw std::logic_error(
            "Trace difference cannot be computed for events other than edge addition or deletion!");
    const auto norm = (col_i - col_j).length();
    const auto lastrowdiff = (col_i[G.numberOfNodes() - 1] - col_j[G.numberOfNodes() - 1]);
    return norm * norm * w * (G.numberOfNodes() - 1)
           - G.numberOfNodes() * w * lastrowdiff * lastrowdiff;
}

void DynFullLaplacianInverseSolver::update(GraphEvent ev) {
    assureFinished();
    assureUpdated(ev);

    INFO("[DynFullLaplacianInverseSolver::update] updating graph with event ", ev.u, ", ", ev.v);

    const node i = ev.u;
    const node j = ev.v;
    const double R_ij = lpinv(i, i) + lpinv(j, j) - 2 * lpinv(i, j);

    if ((1 - R_ij) < 1e-6) { // if close to 1, it is a bridge
        INFO("[DynFullLaplacianInverseSolver::update] identified (", i, ", ", j,
             ") as bridge. Update by recomputing");
        // update the componentDecomposition by running the algorithm again (it is not dynamic)
        decomp.run();
        lpinv = computeLpinv(G, decomp);
    } else {
        lpinv = updateShermanMorrison(ev);
    }

    totalHarmonicResistance.reset();
}

} // namespace NetworKit
