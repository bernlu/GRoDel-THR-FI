import argparse
import networkit as nk
import numpy as np
import pandas as pd
import sys
import itertools
from datetime import datetime


def resistance(Lpinv, u, v):
    return Lpinv[u, u] + Lpinv[v, v] - 2 * Lpinv[u, v]


def totalHarmonicResistance(graph):
    L = nk.algebraic.laplacianMatrix(graph).toarray()
    lpinv = np.linalg.pinv(L)
    components = nk.components.ParallelConnectedComponents(graph).run().getComponents()
    H = 0
    for comp in components:
        for u, v in itertools.product(comp, comp):
            if u >= v:
                continue
            H = H + 1.0 / resistance(lpinv, u, v)
    return H


def findOpt(graph, k=5):
    bestEdges = None
    bestTHR = totalHarmonicResistance(graph)

    goodSolutions = []
    for edges in itertools.combinations(graph.iterEdges(), k):
        # remove edges from graph
        for edge in edges:
            graph.removeEdge(*edge)
        newTHR = totalHarmonicResistance(graph)

        if (bestTHR - newTHR) > -0.1:
            goodSolutions.append({"edges": edges, "totalHarmonicResistance": newTHR})

        if newTHR < bestTHR:
            goodSolutions.clear()
            goodSolutions.append({"edges": edges, "totalHarmonicResistance": newTHR})
            bestEdges = edges
            bestTHR = newTHR
            print("current best:", bestTHR, ", ", bestEdges)
        # add edges again
        for edge in edges:
            graph.addEdge(*edge)
    print("best edges:", bestEdges)
    print("best THR:", bestTHR)
    df = pd.DataFrame(goodSolutions)
    df = df[df.totalHarmonicResistance < (bestTHR + 0.1)]

    return bestEdges, bestTHR, df


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--instance")
    parser.add_argument("-k", type=int)
    parser.add_argument("-o")
    args = parser.parse_args()

    if not args.instance:
        print("error: instance required")
        raise RuntimeError("no instance provided!")

    print("exhaustive search for optimal total harmonic resistance solution.")
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
    bestEdges, bestTHR, df = findOpt(G, args.k)
    after = datetime.now()

    elapsed = after - before
    print("time: ", elapsed)
    print("value: ", bestTHR)
    print("items: ", bestEdges)
    df.to_csv(args.o + f"/{args.instance.split('/')[-1]}-THR.csv")
