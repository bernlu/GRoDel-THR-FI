/*
 *  RobustnessGreedy.hpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */
#ifndef NETWORKIT_ROBUSTNESS_ROBUSTNESS_GREEDY_HPP_
#define NETWORKIT_ROBUSTNESS_ROBUSTNESS_GREEDY_HPP_

#include <networkit/base/Algorithm.hpp>
#include <networkit/graph/Graph.hpp>

#include <harmonicResistance/DynLaplacianInverseSolver.hpp>

#include <optional>

namespace NetworKit {

class RobustnessGreedy : public Algorithm {
public:
    double getResultValue() const;

    const std::vector<Edge> &getResultItems() const {
        assureFinished();
        return result;
    }

protected:
    RobustnessGreedy(Graph &G, count k);
    Graph &G;
    const count k;
    std::unique_ptr<DynLaplacianInverseSolver> lapSolver;

    std::vector<Edge> result;
    double resultValue;

    // precond:
    // - G is augmented with forest node if metric==forest
    std::vector<Edge> buildCandidateSet() const;

    template <class Solver>
    void setupSolver(double solverEpsilon) {
        lapSolver = std::make_unique<Solver>(G, solverEpsilon);
        lapSolver->run();
    }

    template <class Solver>
    void setupSolver() {

        lapSolver = std::make_unique<Solver>(G);
        lapSolver->run();
    }

    void prepareGraph();

    void restoreGraph();
};

} // namespace NetworKit

#endif // NETWORKIT_ROBUSTNESS_ROBUSTNESS_GREEDY_HPP_
