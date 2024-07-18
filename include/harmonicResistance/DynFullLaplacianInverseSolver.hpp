/*
 *  DynFullLaplacianInverseSolver.hpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */
#ifndef NETWORKIT_ROBUSTNESS_DYN_FULL_LAPLACIAN_INVERSE_SOLVER_HPP_
#define NETWORKIT_ROBUSTNESS_DYN_FULL_LAPLACIAN_INVERSE_SOLVER_HPP_

#include <harmonicResistance/DynLaplacianInverseSolver.hpp>
#include <networkit/algebraic/DenseMatrix.hpp>
#include <networkit/components/ParallelConnectedComponents.hpp>

#include <optional>

namespace NetworKit {

class DynFullLaplacianInverseSolver : public DynLaplacianInverseSolver {
public:
    DynFullLaplacianInverseSolver(const Graph &G);

    /// @brief note that this function currently expects the graph to be connected!
    void run() override;
    /// @brief  an update may disconnect the graph
    void update(GraphEvent ev) override;
    virtual double totalResistanceDifference(const GraphEvent &ev) const override;
    virtual double totalForestDistanceDifference(const GraphEvent &ev) const override;
    virtual double totalHarmonicResistanceDifference(const GraphEvent &ev) const override;

    const DenseMatrix &getLpinv() { return lpinv; }

public:
    DenseMatrix lpinv;
    mutable std::optional<double> totalHarmonicResistance;
    ParallelConnectedComponents decomp;

    DenseMatrix updateShermanMorrison(const GraphEvent &ev) const;
    DenseMatrix computeLpinv(const Graph &G, const ComponentDecomposition &decomp) const;
};

} // namespace NetworKit

#endif // NETWORKIT_ROBUSTNESS_DYN_FULL_LAPLACIAN_INVERSE_SOLVER_HPP_
