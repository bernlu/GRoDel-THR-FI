/*
 *  DynLaplacianInverseSolver.hpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */
#ifndef NETWORKIT_ROBUSTNESS_DYN_LAPLACIAN_INVERSE_SOLVER_HPP_
#define NETWORKIT_ROBUSTNESS_DYN_LAPLACIAN_INVERSE_SOLVER_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/base/DynAlgorithm.hpp>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/algebraic/DynamicMatrix.hpp>

namespace NetworKit {

/**
 * computes L+ for the graph G and supports updating L+ for edge insertions and deletions.
 */
class DynLaplacianInverseSolver : public Algorithm, public DynAlgorithm {
public:
    // from Algorithm
    void run() override = 0;

    // from DynAlgorithm
    // assumes that the event has happened before the call to update().
    void update(GraphEvent) override = 0;
    void updateBatch(const std::vector<GraphEvent> &) override {
        throw std::logic_error("updateBatch is not supported!");
    };

    // specific
    virtual double totalResistanceDifference(const GraphEvent &ev) const = 0;
    virtual double totalForestDistanceDifference(const GraphEvent &ev) const = 0;
    virtual double totalHarmonicResistanceDifference(const GraphEvent &ev) const = 0;

protected:
    DynLaplacianInverseSolver(const Graph &G) : G(G) {}
    const Graph &G;

    void assureUpdated(const GraphEvent &ev) const {
        if (ev.type == GraphEvent::EDGE_ADDITION)
            assert(G.hasEdge(ev.u, ev.v));
        if (ev.type == GraphEvent::EDGE_REMOVAL)
            assert(!G.hasEdge(ev.u, ev.v));
    }
};

inline std::ostream &operator<<(std::ostream &os, const DenseMatrix &M) {
    for (index i = 0; i < M.numberOfRows(); i++) {
        for (index j = 0; j < M.numberOfColumns(); j++)
            os << M(i, j) << ", ";
        os << std::endl;
    }
    return os;
}

inline std::ostream &operator<<(std::ostream &os, const CSRMatrix &M) {
    for (index i = 0; i < M.numberOfRows(); i++) {
        for (index j = 0; j < M.numberOfColumns(); j++)
            os << M(i, j) << ", ";
        os << std::endl;
    }
    return os;
}

inline std::ostream &operator<<(std::ostream &os, const DynamicMatrix &M) {
    for (index i = 0; i < M.numberOfRows(); i++) {
        for (index j = 0; j < M.numberOfColumns(); j++)
            os << M(i, j) << ", ";
        os << std::endl;
    }
    return os;
}

} // namespace NetworKit

#endif // NETWORKIT_ROBUSTNESS_DYN_LAPLACIAN_INVERSE_SOLVER_HPP_
