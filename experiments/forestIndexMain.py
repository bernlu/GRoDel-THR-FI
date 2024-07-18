import networkit as nk

import argparse
from datetime import datetime

parser = argparse.ArgumentParser()
parser.add_argument("-i", "--instance")
parser.add_argument("-k", type=int)
parser.add_argument("-p", "--problem")
parser.add_argument("--loglevel")
args = parser.parse_args()

if args.loglevel:
    nk.engineering.setLogLevel(args.loglevel)
else:
    nk.engineering.setLogLevel("INFO")

if not args.instance:
    print("error: instance required")
    raise RuntimeError("no instance provided!")

print("called with:")
print("instance: ", args.instance)
print("k: ", args.k)

G = nk.readGraph(args.instance)
G = nk.components.ConnectedComponents(G).extractLargestConnectedComponent(G, True)
G.removeMultiEdges()
G.removeSelfLoops()

nk.overview(G)


alg = nk.robustness.StGreedy(
    G, args.k, nk.robustness.RobustnessProblem.GLOBAL_REDUCTION
)


before = datetime.now()
alg.run()
after = datetime.now()

elapsed = after - before
print("time: ", elapsed)
print("value: ", alg.getResultValue())
print("items: ", alg.getResultItems())
