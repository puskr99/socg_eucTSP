# Code and data for "Euclidean TSP on Two Crossing Lines in O(n^5) Time, and a Linear Cut Lemma for k Lines"

Everything here is a computational check of statements in the paper; no proof depends on it.

## Requirements

* Python 3.10 (tested with 3.10.12), `numpy` (tested 2.2.6); `networkx` (tested 3.4.2) for `k_lines/`.
* Optional, for speed: a C compiler. Without the two shared libraries below, pure-Python Held-Karp and 2-opt are used
  (fine up to n = 10, slow at n = 12-13).

```
pip install numpy networkx
gcc -O2 -shared -fPIC -o libhk.so hk.c         # Held-Karp
gcc -O2 -shared -fPIC -o libtwo.so two_opt.c   # 2-opt local search
```

* Optional: a Concorde binary (not included), only for the two `concorde_*` scripts; set `CONCORDE_BIN=/path/to/concorde`.
  The stored outputs of these runs are in `*/results/`.

All scripts use paths relative to their own location and can be run from any directory, e.g.
`python two_lines/validate_random.py 1 20 5 10`.

## Layout

| file | what it checks | paper |
|---|---|---|
| `tsp_exact.py`, `hk.c`, `two_opt.c` | Held-Karp optimum, 2-opt local search | used by all |
| `two_lines/dp.py` | the radial-sweep dynamic program (`solve(theta, pts, ct=4)`; caps: 2 along stubs per ray, 2 wedge ends per wedge, `ct` wedge ends in total) | Prop. 15, App. B |
| `two_lines/radial.py` | instances on two lines, radial cut statistics | |
| `two_lines/validate_random.py` | DP = Held-Karp on random instances (uniform, log-radial, clustered, one-sided, integer radii with ties; every 7th at the right angle) | Sec. 3.2 remarks |
| `two_lines/validate_adversarial.py`, `two_lines/data/adversarial_instances.json` | DP = Held-Karp on 7 adversarial 12-point instances and random perturbations of them | Sec. 3.2 remarks |
| `two_lines/cap_test.py` | the 12-point instance at theta = pi/2 with 4 wedge cut edges at one level; DP with total cap 3 is 0.54% too long | Sec. 3.2 remarks |
| `two_lines/pair_witnesses.py` | Table 1 (10 excluded pairs, 14 allowed, triangle-free, Prop. 12 by enumeration), the 14 witnesses of Table 2 (G < 0), random sampling (excluded pairs: min G > 0) | Table 1, Table 2 |
| `two_lines/cut_search.py` | adversarial search (exact optima) for large radial prefix cuts; largest found 6 | Sec. 7 |
| `two_lines/wedge_search.py` | adversarial search on wedge cut edges; objective `sumw` produced the instance of `cap_test.py` | Sec. 3.2 remarks |
| `two_lines/concorde_two_lines.py` (needs Concorde), `two_lines/results/*.json` | radial cuts of Concorde optima, n = 100-1000, angles 3-90 degrees | Sec. 7 |
| `k_lines/halfslot_check.py` | Lemmas 22-24 and Theorem 21 through the half-slot graph H on Held-Karp optima (pencils, triangles with the nets of Lemma 27, nonobtuse triangles, parallel lines with and without a transversal); reports max cut / N | Sec. 5, Sec. 7 |
| `k_lines/halfslot_hill.py` | adversarial hill climbing on the same checks | Sec. 7 |
| `k_lines/halfslot_2opt.py` | the same checks on 2-optimal tours, n = 100-400 | Sec. 5, Sec. 7 |
| `k_lines/concorde_segments.py`, `k_lines/faithful_q.py`, `k_lines/data/segments_concorde_optima.json`, `k_lines/results/concorde_segments_cuts.json` | prefix cuts of stored Concorde optima (n = 40-1000, 2-4 random segments) sorted by the distance to an exact rational faithful net; cut <= N on all 54 instances | Sec. 7 |

## Runs behind the numbers in the paper

The original random seeds were not recorded; any seed reproduces the checks.

* 372 DP vs Held-Karp comparisons (n <= 12): `python two_lines/validate_random.py SEED 260 5 12` (260 instances) and
  `python two_lines/validate_adversarial.py SEED 15` (7 instances x 16 = 112). Expected: 0 mismatches. With the C
  Held-Karp the first takes about 20 minutes; states per level reach about 1.2 million at n = 12.
* Cap 3 vs 4: `python two_lines/cap_test.py` (about 20 s).
* Tables 1 and 2: `python two_lines/pair_witnesses.py 200000 0` (seconds).
* Largest radial cut 6: `python two_lines/cut_search.py SEED 12 300 50` (and n = 13, 14), several seeds.
* k lines, max cut / N <= 1.5: `python k_lines/halfslot_check.py SEED 3200 13`, `python k_lines/halfslot_2opt.py SEED 320`,
  `python k_lines/halfslot_hill.py SEED 24 11`. Expected: all violation counters 0.
* Concorde optima on segments, cut <= N: `python k_lines/concorde_segments.py 1000 4 out.json` (reads the stored optima;
  no Concorde needed).

## Quick smoke test (a few minutes)

```
python two_lines/pair_witnesses.py 100000 0
python two_lines/validate_random.py 1 10 5 8
python two_lines/cap_test.py
python two_lines/validate_adversarial.py 0 0 3,6
python k_lines/halfslot_check.py 0 40 10
python k_lines/halfslot_2opt.py 0 2
python k_lines/concorde_segments.py 80 3
```
