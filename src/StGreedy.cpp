/*
 *  stGreedy.cpp
 *
 *  Created on: 27.06.2023
 *     Authors: Maria Predari <predarimaria@gmail.com>
 *              Lukas Berner <Lukas.Berner@hu-berlin.de>
 */

#include <harmonicResistance/StGreedy.hpp>

#include <harmonicResistance/DynFullLaplacianInverseSolver.hpp>

#include <networkit/auxiliary/Timer.hpp>

#include <harmonicResistance/SubmodularGreedy.hpp>

namespace NetworKit {

StGreedy::StGreedy(Graph &G, count k) : RobustnessGreedy(G, k) {}

void StGreedy::run() {
    // Compute pseudoinverse of laplacian

    result.clear();
    resultValue = 0;

    // setup solver
    Aux::StartedTimer timer;
    setupSolver<DynFullLaplacianInverseSolver>();
    timer.stop();
    INFO("solver setup time: ", timer.elapsedTag());
    // candidates
    std::vector<Edge> items = buildCandidateSet();

    SubmodularGreedy<Edge> greedy(items, k);

    uint64_t solverTime = 0;
    count gainCalls = 0;

    greedy.setGainFunction([&](const Edge &e) {
        GraphEvent ev(GraphEvent::EDGE_REMOVAL, e.u, e.v);
        gainCalls++;
        double gain = 0;
        Aux::StartedTimer timer;
        gain = lapSolver->totalHarmonicResistanceDifference(ev);
        timer.stop();
        solverTime += timer.elapsedMicroseconds();
        return gain;
    });
    greedy.setPickedItemCallback([&](const Edge &e) {
        INFO("Edge picked: (", e.u, ", ", e.v, ")");
        GraphEvent ev(GraphEvent::EDGE_REMOVAL, e.u, e.v);
        G.removeEdge(ev.u, ev.v);
        Aux::StartedTimer timer;
        lapSolver->update(ev);
        timer.stop();
        solverTime += timer.elapsedMicroseconds();
    });

    greedy.run();

    INFO("gain calls: ", gainCalls);
    INFO("time spend in solver (microseconds): ", solverTime);

    resultValue = greedy.getResultValue();
    result = greedy.getResultItems();

    // reset G to original state
    restoreGraph();

    this->hasRun = true;
}

} // namespace NetworKit
