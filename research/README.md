# research/ -- exploratory, NOT production code

Everything in this folder is a standalone development/diagnostic script
from real attempts to fully close two open items in `LIMITATIONS.md`:

- `oguchi_pair_cluster_prototype.py` -- an Oguchi/Bethe-Peierls two-spin
  cluster correction to the mean-field Gd model (Items 1.1 / 1.1b).
  Internally validated (passes its own 3-check validation plan) but
  **not adopted**: checked against Dan'kov et al.'s own published
  DeltaT_ad values, and does not resolve the target discrepancy any
  better than the existing mean-field model. See the file's own
  docstring and LIMITATIONS.md Item 1.1c for the full finding, including
  why a *second*, independent approach (Franco-style critical-exponent
  "universal curve" scaling) was also ruled out before being built, using
  data already in this repo.

- `regenerator_1d_v2_coupled_prototype.py` -- a properly-coupled
  (single implicit tridiagonal solve per timestep, instead of
  operator-split) 1-D AMR regenerator solver, built to test whether
  operator-splitting order explains the span-undershoot in
  `core/regenerator_1d.py` (Item 1.3 / issue #8). It does not: span is
  still ~0 across a realistic NTU range, and (independently re-verified
  here) the script isn't even numerically stable over many cycles in
  its current form.

**Neither file is imported by anything in `core/` or `main.py`, and
neither is covered by `tests/`.** They are kept here, not deleted,
because a documented dead end is worth more than a re-walked one --
each file's docstring records exactly what was tried, what was ruled
out, and why, so future work on these two items doesn't have to
rediscover the same four-plus eliminated hypotheses from scratch.

See `LIMITATIONS.md` Items 1.1c and 1.3 for the narrative summary.
