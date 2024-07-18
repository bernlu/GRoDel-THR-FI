import argparse
import networkit as nk
import numpy as np
import pandas as pd
import sys
import itertools
from datetime import datetime


def forestIndex(graph):
    n = graph.numberOfNodes()
    L = nk.algebraic.laplacianMatrix(graph).toarray()
    Omega = np.linalg.inv(np.eye(n) + L)
    return n * np.trace(Omega) - n


def findOpt(graph, k=5):
    bestEdges = None
    bestFI = forestIndex(graph)

    goodSolutions = []
    for edges in itertools.combinations(graph.iterEdges(), k):
        # remove edges from graph
        for edge in edges:
            graph.removeEdge(*edge)
        newFI = forestIndex(graph)

        if (newFI - bestFI) > -0.1:
            goodSolutions.append({"edges": edges, "forestIndex": newFI})

        if newFI > bestFI:
            goodSolutions.clear()
            goodSolutions.append({"edges": edges, "forestIndex": newFI})
            bestEdges = edges
            bestFI = newFI
            print("current best:", bestFI, ", ", bestEdges)
        # add edges again
        for edge in edges:
            graph.addEdge(*edge)
    print("best edges:", bestEdges)
    print("best FI:", bestFI)
    df = pd.DataFrame(goodSolutions)
    df = df[df.forestIndex < (bestFI + 0.1)]

    return bestEdges, bestFI, df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--instance")
    parser.add_argument("-k", type=int)
    parser.add_argument("-o")
    args = parser.parse_args()

    if not args.instance:
        print("error: instance required")
        raise RuntimeError("no instance provided!")

    print("exhaustive search for optimal forest index solution.")
    print("called with:")
    print("instance: ", args.instance)
    print("k: ", args.k)
    print("o: ", args.o)

    G = nk.readGraph(args.instance)
    G = nk.components.ConnectedComponents(G).extractLargestConnectedComponent(G, True)
    G.removeMultiEdges()
    G.removeSelfLoops()

    nk.overview(G)

    before = datetime.now()
    bestEdges, bestFI, df = findOpt(G, args.k)
    after = datetime.now()

    elapsed = after - before
    print("time: ", elapsed)
    print("value: ", bestFI)
    print("items: ", bestEdges)
    df.to_csv(args.o + f"/{args.instance.split('/')[-1]}-FI.csv")
