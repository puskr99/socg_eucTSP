# Code and data for "Euclidean TSP on Two Crossing Lines in O(n^5) Time, and a Linear Cut Lemma for k Lines"

Everything here is a computational check of statements in the paper; no proof depends on it (the signs in Table 3 are
also proved by hand in Appendix B; `exact_witnesses.py` re-checks them exactly).

## Environment

Tested with:

* Python 3.10.12 (CPython, Linux x86-64)
* packages (subset of `pip freeze`):
  ```
  numpy==2.2.6
  networkx==3.4.2      # k_lines/ only
  sympy==1.14.0        # two_lines/exact_witnesses.py only
  mpmath==1.3.0        # dependency of sympy
  ```
  `pip install numpy==2.2.6 networkx==3.4.2 sympy==1.14.0`
* C compiler: gcc 11.4.0 (Ubuntu 22.04). The two optional shared libraries are built in this directory with
  ```
  gcc -O2 -shared -fPIC -o libhk.so hk.c         # Held-Karp
  gcc -O2 -shared -fPIC -o libtwo.so two_opt.c   # 2-opt local search
  ```
  `tsp_exact.py` loads them if present and otherwise falls back to pure Python (fine up to n = 10, slow at n = 12-14).
  All run times below are with the libraries, on one core of an ordinary desktop PC.
* Optional: a Concorde binary (not included), only for `two_lines/concorde_two_lines.py`; set
  `CONCORDE_BIN=/path/to/concorde`. The stored outputs of these runs are in `*/results/`.

All scripts use paths relative to their own location and can be run from any directory. Commands below are run from
this directory.

## Generators

All instances are generated inside the scripts with `numpy.random.default_rng(SEED)` (the first argument), except
`concorde_two_lines.py`, which uses `random.Random(SEED)`.

* `two_lines/validate_random.py`: n uniform in [NMIN, NMAX], angle uniform in [0.02, pi/2] (every 7th instance exactly
  pi/2); signed positions on the two lines cycle through five modes: uniform in [-1, 1], log-radial
  (+-exp(U[-4, 0])), clustered (+-(1, 2 or 3) + U[0, 0.05]), one-sided (U[0, 1]), integer radii 1..5 (ties across rays);
  each point on a uniformly random line (both lines used).
* `two_lines/validate_adversarial.py`: the 7 stored adversarial 12-point instances
  (`two_lines/data/adversarial_instances.json`, produced by `wedge_search.py` / `cut_search.py`) and, per instance,
  PERTURBATIONS random perturbations (angle + N(0, 0.02), radii * exp(N(0, 0.05)), every 5th: one point to the
  opposite ray).
* `two_lines/cut_search.py`, `two_lines/wedge_search.py`: hill climbing from random instances (angle U[0.05, pi/2],
  positions U[-1, 1], random lines), moves: angle, one position, one line label.
* `k_lines/halfslot_check.py`, `halfslot_2opt.py`, `halfslot_hill.py`: random arrangements (pencils of 2-4 concurrent
  lines, triangles with the net of Lemma 27, nonobtuse triangles, two parallel lines with a transversal, 2-4 parallel
  lines) and random points on them (function `make`, `sample_points` in `halfslot_check.py`).
* `k_lines/concorde_segments.py`: no generator; reads the stored instances and Concorde optima in
  `k_lines/data/segments_concorde_optima.json`.

## Seeds

The random seeds of the original runs behind the numbers in the paper were not recorded. We therefore fixed the seed
2027 for every randomised experiment and reran them with this seed (except the Concorde runs, see below). All rows except the Concorde ones were rerun and confirm the
numbers stated in the paper.

## Experiments: exact commands, expected output, run time

| # | command | paper | expected output (seed 2027) | time |
|---|---|---|---|---|
| 1 | `python two_lines/exact_witnesses.py` | Table 3 | `all 14 witnesses certified exactly: True` (exit status 0); G <= -0.754 for all 14 | 2 s |
| 2 | `python two_lines/pair_witnesses.py 200000 2027` | Tables 2, 3, Prop. 12 | 10 excluded / 14 allowed pairs, `triangles in the allowed graph: []`, `largest same-orientation set satisfying (F1)-(F4): 2`, all 14 witnesses `OK`, `witnesses cover exactly the allowed pairs: True`; sampled min G > 0 for every excluded pair, < 0 for every allowed pair | 2 s |
| 3 | `python two_lines/cap_test.py` (deterministic) | Sec. 3.2 remarks | `max wedge cut edges / level 4`; cap 4 = Held-Karp optimum 2.036277995; cap 3: 2.047357472 (+0.54%) | 70 s |
| 4 | `python two_lines/validate_random.py 2027 260 5 12` | Sec. 3.2 remarks (260 of the 372) | `done 260 mismatches 0` | about 30 min (rerun with seed 2027: `done 260 mismatches 0`, 1819 s, peak memory 1.9 GB) |
| 5 | `python two_lines/validate_adversarial.py 2027 15` | Sec. 3.2 remarks (7 x 16 = 112 of the 372) | `TOTAL 112 mismatches 0 max wedge cut edges at one level 4` | about 100 min (rerun with seed 2027: `TOTAL 112 mismatches 0 max wedge cut edges at one level 4`, 5875 s, peak memory 1.9 GB) |
| 6 | `python two_lines/cut_search.py 2027 N 300 50` for N = 12, 13, 14 | Sec. 1, Sec. 7 (largest radial cut 6) | `best max cut 6` for each N (rerun with seed 2027) | 16 s / 35 s / 74 s for N = 12 / 13 / 14 |
| 7 | `python k_lines/halfslot_check.py 2027 3200 13` | Sec. 5, Sec. 7 (3,200 exact optima) | all violation counters 0, `maxcutN` 1.33 | 32 s |
| 8 | `python k_lines/halfslot_2opt.py 2027 320` | Sec. 5, Sec. 7 (320 2-optimal tours) | all violation counters 0, `maxcutN` 1.5 | 160 s |
| 9 | `python k_lines/halfslot_hill.py 2027 24 11` | Sec. 7 (adversarial) | no violation reported, `best 0.17` (largest cut / (4(N+Lambda_t)-8)) | 90 s |
| 10 | `python k_lines/concorde_segments.py 1000 4 out.json` | Sec. 7 (Concorde optima on segments) | 54 instances, all `ok`, `maxcross <= pieces` on all 54, `improving_pairs` 0 | 25 s |

The bound "at most 1.5N" for k lines in Section 7 is the largest value of `maxcutN` over all runs of rows 7-9 (and of
the original runs); the bound 8N - 8 is never reached (`maxratio` = cut / (4(N + Lambda_t) - 8) stays far below 1).

`max states/level` in rows 4-5 is the largest number of DP states at one level (about 1.2 million at n = 12); peak
memory is below 300 MB.

### Concorde runs (not rerun)

`python two_lines/concorde_two_lines.py SEED OUT.json 100,200,500,1000` needs Concorde (via `concorde_tsp.py`,
coordinates scaled by 1e7, TSPLIB EUC_2D). Its original seed was not recorded; the outputs used in Section 7
(largest radial cut 6 on Concorde optima, n = 100-1000, angles 3-90 degrees) are stored in
`two_lines/results/concorde_maxcut.json` and `two_lines/results/concorde_wedge_stats.json`. Row 10 recomputes the
cuts of the stored Concorde optima on segments without Concorde.

## Layout

| file | what it checks | paper |
|---|---|---|
| `tsp_exact.py`, `hk.c`, `two_opt.c` | Held-Karp optimum, 2-opt local search | used by all |
| `concorde_tsp.py` | optional Concorde wrapper | Sec. 7 |
| `two_lines/dp.py` | the radial-sweep dynamic program (`solve(theta, pts, ct=4)`; caps: 2 along stubs per ray, 2 wedge ends per wedge, `ct` wedge ends in total) | Prop. 15, App. B |
| `two_lines/radial.py` | instances on two lines, radial cut statistics | |
| `two_lines/validate_random.py` | DP = Held-Karp on random instances | Sec. 3.2 remarks |
| `two_lines/validate_adversarial.py`, `two_lines/data/adversarial_instances.json` | DP = Held-Karp on 7 adversarial 12-point instances and random perturbations of them | Sec. 3.2 remarks |
| `two_lines/cap_test.py` | the 12-point instance at theta = pi/2 with 4 wedge cut edges at one level; DP with total cap 3 is 0.54% too long | Sec. 3.2 remarks |
| `two_lines/pair_witnesses.py` | Table 2 (10 excluded pairs, 14 allowed, triangle-free, Prop. 12 by enumeration), the 14 witnesses of Table 3 in floating point, random sampling (excluded pairs: min G > 0) | Tables 2, 3 |
| `two_lines/exact_witnesses.py` | the 14 witnesses of Table 3 exactly: closed forms of G (sympy radicals), max(x,u) < min(y,v) in rationals, G < 0 by rigorous rational interval bounds (integer square roots) and independently by exact root isolation of the minimal polynomial of G | Table 3 |
| `two_lines/cut_search.py` | adversarial search (exact optima) for large radial prefix cuts | Sec. 1, Sec. 7 |
| `two_lines/wedge_search.py` | adversarial search on wedge cut edges; objective `sumw` produced the instance of `cap_test.py` | Sec. 3.2 remarks |
| `two_lines/concorde_two_lines.py` (needs Concorde), `two_lines/results/*.json` | radial cuts of Concorde optima, n = 100-1000, angles 3-90 degrees | Sec. 7 |
| `k_lines/halfslot_check.py` | Lemmas 22-24 and Theorem 21 through the half-slot graph H on Held-Karp optima; reports max cut / N | Sec. 5, Sec. 7 |
| `k_lines/halfslot_hill.py` | adversarial hill climbing on the same checks | Sec. 7 |
| `k_lines/halfslot_2opt.py` | the same checks on 2-optimal tours, n = 100-400 | Sec. 5, Sec. 7 |
| `k_lines/concorde_segments.py`, `k_lines/faithful_q.py`, `k_lines/data/segments_concorde_optima.json`, `k_lines/results/concorde_segments_cuts.json` | prefix cuts of stored Concorde optima (n = 40-1000, 2-4 random segments) sorted by the distance to an exact rational faithful net; cut <= N on all 54 instances | Sec. 7 |

## Quick smoke test (about 3 minutes)

```
python two_lines/exact_witnesses.py
python two_lines/pair_witnesses.py 100000 0
python two_lines/validate_random.py 1 10 5 8          # done 10 mismatches 0
python two_lines/cap_test.py
python two_lines/validate_adversarial.py 0 0 3,6      # TOTAL 2 mismatches 0
python k_lines/halfslot_check.py 0 40 10
python k_lines/halfslot_2opt.py 0 2
python k_lines/concorde_segments.py 80 3
```
