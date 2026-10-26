#!/usr/bin/env python3
# Copyright (c) 2026 Raban Heller. MIT licence; see ../../../LICENSE.md.
"""Exact tie classes of the Boschloo ordering statistics.

For every table (x1, x2) of a design (n1, n2) the Fisher conditional p-value
given the total t = x1 + x2 is computed as an exact fraction a / C(N, t):

  lower:   P(X1 <= x1 | t)
  upper:   P(X1 >= x1 | t)
  minlike: sum of P(X1 = y | t) over y with P(X1 = y | t) <= P(X1 = x1 | t)

using Python integers (math.comb) and fractions.Fraction. Tables are ranked
from the most extreme value (rank 1) upwards; exactly equal values share a
rank. Nothing here uses floating-point arithmetic, so the ranks are the
exact tie classes against which the engine's tolerance-based levels are
checked.

Usage:  exact_classes.py OUTDIR n1:n2 [n1:n2 ...]

Writes OUTDIR/ranks_<n1>_<n2>.csv (columns x1, x2, rank_lower, rank_upper,
rank_minlike; x1 varies fastest) and appends one row per design and side to
OUTDIR/gaps.csv with the number of distinct values, the number of tables
tied with another table, and the smallest relative gap between consecutive
distinct values (on the value scale for values <= 1/2 and on the complement
scale 1 - value for values > 1/2).
"""

import csv
import os
import sys
from fractions import Fraction
from math import comb


def ordering_values(n1, n2):
    """Exact ordering statistics of all tables, x1 varying fastest."""
    big_n = n1 + n2
    tables = [(x1, x2) for x2 in range(n2 + 1) for x1 in range(n1 + 1)]
    weight = {(x1, x2): comb(n1, x1) * comb(n2, x2) for (x1, x2) in tables}
    strata = {}
    for (x1, x2) in tables:
        strata.setdefault(x1 + x2, []).append((x1, x2))
    values = {"lower": {}, "upper": {}, "minlike": {}}
    for t, members in strata.items():
        denom = comb(big_n, t)
        members = sorted(members)                       # ascending x1
        ws = [weight[m] for m in members]
        assert sum(ws) == denom                         # Vandermonde
        run = 0
        for m, w in zip(members, ws):
            run += w
            values["lower"][m] = Fraction(run, denom)
        run = 0
        for m, w in zip(reversed(members), reversed(ws)):
            run += w
            values["upper"][m] = Fraction(run, denom)
        for m, w in zip(members, ws):
            values["minlike"][m] = Fraction(sum(v for v in ws if v <= w), denom)
    return tables, values


def dense_ranks(tables, vals):
    distinct = sorted(set(vals[m] for m in tables))
    position = {v: i + 1 for i, v in enumerate(distinct)}
    return [position[vals[m]] for m in tables], distinct


def min_relative_gaps(distinct):
    half = Fraction(1, 2)
    low = [v for v in distinct if v <= half]
    high = [1 - v for v in distinct if v > half]      # complement scale, descending in v
    high.sort()
    def smallest(seq):
        best = None
        for a, b in zip(seq, seq[1:]):
            if a > 0:
                g = (b - a) / a
                best = g if best is None or g < best else best
        return best
    return smallest(low), smallest(high)


def main(argv):
    if len(argv) < 3:
        sys.stderr.write(__doc__)
        return 2
    outdir = argv[1]
    os.makedirs(outdir, exist_ok=True)
    gap_path = os.path.join(outdir, "gaps.csv")
    new_gap_file = not os.path.exists(gap_path)
    with open(gap_path, "a", newline="") as gf:
        gw = csv.writer(gf)
        if new_gap_file:
            gw.writerow(["n1", "n2", "side", "tables", "distinct_values", "tables_tied",
                         "min_rel_gap_low", "min_rel_gap_high_complement"])
        for spec in argv[2:]:
            n1, n2 = (int(z) for z in spec.split(":"))
            tables, values = ordering_values(n1, n2)
            ranks = {}
            for side in ("lower", "upper", "minlike"):
                r, distinct = dense_ranks(tables, values[side])
                ranks[side] = r
                counts = {}
                for v in values[side].values():
                    counts[v] = counts.get(v, 0) + 1
                tied = sum(c for c in counts.values() if c > 1)
                g_low, g_high = min_relative_gaps(distinct)
                gw.writerow([n1, n2, side, len(tables), len(distinct), tied,
                             "" if g_low is None else "%.6e" % float(g_low),
                             "" if g_high is None else "%.6e" % float(g_high)])
            with open(os.path.join(outdir, "ranks_%d_%d.csv" % (n1, n2)), "w", newline="") as f:
                w = csv.writer(f)
                w.writerow(["x1", "x2", "rank_lower", "rank_upper", "rank_minlike"])
                for i, (x1, x2) in enumerate(tables):
                    w.writerow([x1, x2, ranks["lower"][i], ranks["upper"][i], ranks["minlike"][i]])
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv))
