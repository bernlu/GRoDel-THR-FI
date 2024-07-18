#include <harmonicResistance/lib.hpp>

#include <harmonicResistance/StGreedy.hpp>
#include <networkit/auxiliary/Timer.hpp>
#include <networkit/components/ConnectedComponents.hpp>
#include <networkit/graph/Graph.hpp>
#include <networkit/io/GraphReader.hpp>
#include <networkit/io/METISGraphReader.hpp>
#include <networkit/io/NetworkitBinaryReader.hpp>

using NetworKit::node;

int main(int argc, char *argv[]) {
    std::vector<std::string> args(argv, argv + argc);

    int k = -1;
    for (int i = 0; i < args.size(); ++i) {
        if (args[i] == "-k") {
            k = std::stoi(args[i + 1]);
        }
    }
    if (k == -1) {
        std::cout << "Error: missing argument -k";
        return 1;
    }

    std::cout << "harmonicResistanceGreedy called with:\n" << args[1] << "\nk = " << k << "\n";

    auto G = NetworKit::NetworkitBinaryReader().read(args[1]);

    NetworKit::ConnectedComponents comp(G);
    G = comp.extractLargestConnectedComponent(G, true);
    G.removeMultiEdges();
    G.removeSelfLoops();

    std::cout << "n: " << G.numberOfNodes() << ", m: " << G.numberOfEdges() << "\n";

    NetworKit::StGreedy solver(G, k);

    Aux::Timer t;
    t.start();
    solver.run();
    t.stop();

    std::cout << "Running time: " << t.elapsedMilliseconds() << "ms\n";

    std::cout << "Result edges: ";
    for (auto &edge : solver.getResultItems())
        std::cout << "(" << edge.u << ", " << edge.v << "), ";
    std::cout << "\n";

    return 0;
}
