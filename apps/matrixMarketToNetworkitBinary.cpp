#include <networkit/algebraic/MatrixTools.hpp>
#include <networkit/io/MatrixMarketReader.hpp>
#include <networkit/io/NetworkitBinaryWriter.hpp>

int main(int argc, char *argv[]) {
    std::vector<std::string> args(argv, argv + argc);

    auto matrix = NetworKit::MatrixMarketReader().read(args[1]);
    auto graph = MatrixTools::matrixToGraph(matrix);

    auto path_without_suffix = args[1].substr(0, args[1].rfind("."));

    NetworKit::NetworkitBinaryWriter().write(graph, path_without_suffix + ".nkb");

    return 0;
}
