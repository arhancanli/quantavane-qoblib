# quantavane-qoblib

Search programs and the independent checker behind two submissions to QOBLIB, the Quantum Optimization
Benchmarking Library (https://github.com/ZIB-AOPT/QOBLIB):

- `10-topology/submissions/20261102_EdgeSwitchSA_Canli` (order/degree problem: smallest diameter for n nodes and
  maximum degree d), made with `src/qob_topo.c`;
- `07-independentset/submissions/20261102_NuMVC_Canli` (maximum stable set), made with `src/qob_numvc.c`.

The results, per-instance summaries and objective time series are in those submission folders. This repository
holds the code and the exact commands, so every run can be repeated.

## Contents

| Path | What it is |
| :--- | :--------- |
| `src/qob_topo.c` | Order/degree search: degree-preserving edge switches, simulated annealing, descending diameter target, Moore lower bound as cutoff |
| `src/qob_numvc.c` | Stable sets via vertex cover local search after NuMVC (Cai, Su, Luo and Sattar, JAIR 46, 2013), greedy clique-cover upper bound as cutoff |
| `check/check.py` | Independent checker (Python 3, standard library, integers only) |
| `runs/commands.txt` | Every run behind the submissions: instance, seed, limits and the exact command line |

## Build

```sh
cc -O3 -o qob_topo src/qob_topo.c -lm
cc -O3 -o qob_numvc src/qob_numvc.c
```

Tested with Apple clang 21.0.0 on macOS 26.6.2 (Apple M5 Pro). Both programs are single-threaded C with no
dependencies.

## Run

Paths to QOBLIB instance files assume a clone of QOBLIB in `QOBLIB/`.

```sh
# order/degree: n, d, CPU budget (s), stall limit (CPU s without a new incumbent), seed, graph out, time series out,
# then sweep (moves per cooling round), T0, T1, weight of far pairs
./qob_topo 30 6 150 60 1 topology_30_6.gph topology_30_6.json 10000000 8 0.5 1

# stable set: graph, CPU budget (s), seed, write-threshold (0 = always write), solution out,
# gamma (0 = n/2), rho (0 = 0.3), extra stop size (0 = none), time series out
./qob_numvc QOBLIB/07-independentset/instances/frb100-40.gph 30 1 0 frb100-40.sol 0 0 0 frb100-40.json
```

`qob_topo` starts from a scrambled circulant graph, takes its diameter as the first incumbent and then lowers the
target diameter by one each time it finds a graph that meets it. It computes the Moore lower bound (the smallest k
with 1 + d + d(d-1) + ... + d(d-1)^(k-1) >= n), never targets below it, and stops when the incumbent reaches it:
the graph is then proven optimal. Otherwise a run ends at its CPU budget or after the stall limit without a new
incumbent.

`qob_numvc` first computes a greedy clique cover. A stable set has at most one vertex per clique, so the number of
cliques is an upper bound, and the run stops when the stable set reaches it (proven optimal), or at its CPU budget.
Neither program uses a known value. Both print their bound and write the solution in QOBLIB's format, plus the
objective time series in QOBLIB's JSON format: one run, with the first incumbent, every improvement and the end of
the run.

## Check

```sh
python3 check/check.py topology 30 6 topology_30_6.gph 2        # recomputes the diameter by BFS from every node
python3 check/check.py mis QOBLIB/07-independentset/instances/frb100-40.gph frb100-40.sol 97
```

The checker is a separate code path from the searches. It reads only the instance and the solution file and
recomputes the objective from its definition. It prints `OK <kind> <value>`, and exits non-zero when the file is
invalid or the value differs from the optional claim. It also handles LABS sequences and network design solutions.

## Timing

Both programs are single-threaded. They measure and report the process's own CPU time (user plus system): budgets,
stall limits and the times in the time series are CPU seconds. That is the time a run takes with a core to itself.
The runs behind the submissions shared a busy 18-core machine (load about 70 to 150) at nice 10, so their wall-clock
times also contain waiting for a core; that waiting is excluded from every reported time.

## Earlier builds

The longer stable-set runs on frb59-26-2 (seeds 2 to 5, 1800 s wall limits) and frb100-40 (seeds 2 to 5, 3600 s, and
6 and 7, 7200 s wall limits) were made with earlier builds of `qob_numvc.c`. Those builds had the same search, but no
bound, a wall-clock limit and no time-series file. Their CPU times come from their logs, which record the CPU time of
every new best set; the greedy start, within the first 0.1 s, has only a wall-clock reading. The first seed-1 runs
on both instances had no CPU readings at all, so seed 1 was rerun with the current build, which has the same seed and
the same search and so follows the same path. The commands of all of these runs are in `runs/commands.txt`.

## AI disclosure

The code and this text were written with Claude (Anthropic), model Claude Opus 5.5 (`claude-opus-5-5`), in Claude
Code, and reviewed by Arhan Canli, who is responsible for them.

## License

MIT (see `LICENSE`).
