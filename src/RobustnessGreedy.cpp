/*
 *  RobustnessGreedy.cpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <harmonicResistance/RobustnessGreedy.hpp>

#include <networkit/algebraic/CSRMatrix.hpp>
#include <networkit/numerics/LAMG/Lamg.hpp>

namespace NetworKit {

RobustnessGreedy::RobustnessGreedy(Graph &G, count k) : G(G), k(k) {}

double RobustnessGreedy::getResultValue() const {
    assureFinished();
    return resultValue;
}

std::vector<Edge> RobustnessGreedy::buildCandidateSet() const {
    std::vector<Edge> items;

    G.forEdges([&](node u, node v) { items.push_back(Edge(u, v)); });
    return items;
}

void RobustnessGreedy::restoreGraph() {
    // reset G to original state
    for (auto &edge : result)
        G.addEdge(edge.u, edge.v);
}

} // namespace NetworKit
