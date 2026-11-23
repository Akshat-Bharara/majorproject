# Towards Claim 2 of Singireddy-Basappa-Mitchell (the two-sided search tree)

Tools behind Section 6.3 of `../paper/paper.pdf`.

| File | Content |
|---|---|
| `probe.cpp` | statistics of the two-sided search: successful roots, dead states, closest-pair roots |
| `strat.cpp`, `strat2.cpp` | the witness-preserving search (only safe steps); root rules |
| `stuck.cpp` | classification of stuck states (both steps unsafe) |
| `hunt_cp.cpp`, `hunt_sa.cpp`, `anneal_stuck.cpp`, `anneal_global.cpp` | simulated annealing against the root rules |
| `firststep.cpp` | sanity test of Lemma 18 / Theorem 20 (attachment, safe first step) |
| `check_roots.py`, `data/root_rule_instances.json` | independent exact check of Proposition 22 and of a vertex with no safe first step |

Status: Claim 2 is proved for k <= 3 and under the arc property; the first
step can always be made safe (Theorem 20). The general claim and the sharper
Conjecture 23 remain open; natural root rules (closest-pair endpoint, smallest
angle, globally smallest angle) provably do not work (Proposition 22).
