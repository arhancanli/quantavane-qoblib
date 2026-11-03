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
| `src/qob_topo.c` | Order/degree search: degree-preserving edge switches, simulated annealing, descending diameter target |
| `src/qob_numvc.c` | Stable sets via vertex cover local search after NuMVC (Cai, Su, Luo and Sattar, JAIR 46, 2013) |
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

# stable set: graph, wall-clock limit (s), seed, write-threshold (0 = always write), solution out,
# gamma (0 = n/2), rho (0 = 0.3), stop size (0 = none), time series out
./qob_numvc QOBLIB/07-independentset/instances/frb100-40.gph 30 1 0 frb100-40.sol 0 0 0 frb100-40.json
```

`qob_topo` starts from a scrambled circulant graph, takes its diameter as the first incumbent and then lowers the
target diameter by one each time it finds a graph that meets it. A run ends at its CPU budget or after the stall
limit without a new incumbent. `qob_numvc` runs to its wall-clock limit. Neither program uses a bound or a known
value to stop. Both write the solution in QOBLIB's format and the objective time series in QOBLIB's JSON format:
one run, with the first incumbent, every improvement and the end of the run.

## Check

```sh
python3 check/check.py topology 30 6 topology_30_6.gph 2        # recomputes the diameter by BFS from every node
python3 check/check.py mis QOBLIB/07-independentset/instances/frb100-40.gph frb100-40.sol 97
```

The checker is a separate code path from the searches. It reads only the instance and the solution file and
recomputes the objective from its definition. It prints `OK <kind> <value>`, and exits non-zero when the file is
invalid or the value differs from the optional claim. It also handles LABS sequences and network design solutions.

## Timing

The runs shared a busy 18-core machine with other jobs, at nice 10. A process was often ready to run but waiting
for a free core, so its CPU time is lower than its wall-clock time. `qob_topo` budgets are in CPU seconds and
`qob_numvc` limits are in wall-clock seconds, as stated in `runs/commands.txt`.

## Earlier builds

The longer stable-set runs (frb59-26-2: seeds 1-5, limit 3600 s for seed 1 and 1800 s for seeds 2-5; frb100-40:
3600 s for seeds 1-5 and 7200 s for seeds 6-7) were made with earlier builds of `qob_numvc.c`. Those builds had the same search but wrote no time-series file,
and the builds used for seed 1 of each instance did not report CPU time. Their time series were rebuilt from the
run logs, which record every new incumbent with its time. frb59-26-2 seed 1 had a 3600 s limit and was ended by hand
at 1186 s, after it had reached 59 at 1086 s.

One output fix was made during the 30-second stable-set runs. The first build of `qob_numvc.c` used for them did not
save the starting set when that set was already the run's best. The 54 runs where this happened, all on small
instances, were rerun with the fixed build, and the reruns are the ones reported. The search is the same in both
builds.

## AI disclosure

The code and this text were written with Claude (Anthropic), model Claude Opus 5.5 (`claude-opus-5-5`), in Claude
Code, and reviewed by Arhan Canli, who is responsible for them.

## License

MIT (see `LICENSE`).
