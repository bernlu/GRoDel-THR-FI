# total harmonic resistance for edge deletions

code for the paper 'Introducing Total Harmonic Resistance for Graph Robustness under Edge Deletions' ([arxiv](https://arxiv.org/abs/2407.11521)).

This repository contains the greedy algorithm for THR and our experimental setup (via simexpal)

# how to run experiments

## generating / loading instances
The case study instances can be generated using the [generate-osm.py](apps/generate-osm.py) script.
Requirements for pip can be found in the top-level [requirements.txt](requirements.txt).
Additionally, the NetworKit version found in `extern/networkit` has to be installed via `pip install extern/networkit`.

All other instances can be found on snap, networkrepository, and/or konect. Before running the experiments, all files have to be converted to the networkit binary format (eg by just reading and saving all files once with networkit in python; for the MatrixMarket format, a converter tool is available in [apps](apps)).

## launching experiments

1. navigate to experiments dir
2. install python requirements: `pip install -r requirements.txt` (you should use a different venv to the one used for instance installing!)
3. build with simexpal `simex b make`
4. launch experiments `simex e launch`
5. evaluate: run the jupyter notebooks contained in this directory.