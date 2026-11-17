# v33 — design record

## What v33 is

v33 is a combine-tuning round on top of the verified v32 cell
(practice 0.117532, LOO 0.208604/0.210005). It adds two genuinely
new levers that v3–v32 never tried, and both pay:

1. **A fifth member (mv).** The per-hour median across four views
   (OLS, console shape, all shape, clone forecast) was built in v22
   as an experiment and never used as a combine member. This round
   validates it against its own cached rows (0.000000 on practice
   AND LOO) and uses it as the mid-band median source:
   median(Y5, clone, mv). The mid median moves the median-only
   practice loss from 0.182041 to 0.180454 on its own.
2. **A range disagreement measure.** The disagreement sd is now
   (max − min) / 2 of the five member medians, not the population
   sd. Hours where one member is an outlier get a bigger shift.
   A pure rescale of the population sd (×1.3) scores 0.117032 —
   worse than the range measure — so the profile, not just the
   scale, is doing work. Range/2.5 (a smaller effective sd) is the
   most LOO-robust variant (0.207718/0.208890) but costs practice.

The packaged cell also deepens the mid ramp to −0.18 (steepest ever
used), raises the early ramp to +0.055, and re-weights the spread
blend to 45% Y5 + 45% Z0 + 10% clone (less clone weight helps both
gates). The mid ramp wall is at −0.18: −0.19/−0.20 flatten practice
(0.116468/0.116460) while p01 — the most LOO-typical practice case —
worsens sharply (0.1371/0.1379).

## The search (all through `combine` on the cached member rows)

Cached rows: rows29/prows29 (4 members) + mv rows built this round
(work/mv-rows.json, mv-prows.json — validated 0.000000 vs
rows22.json). Exact practice lnV1 are used in the sweep (the rounded
10.059/11.666/12.504/11.688 in earlier harnesses shift cached-row
scores ~1e-6 from direct grades; the direct grade is the authority,
and it matched v32 exactly: 0.11753273644054554).

- **Phase A — new combine levers on the 4-member rows.** 10 cells
  pass. Best: sd range/2 (0.117334 / 0.208302 / 0.209509), sd sample
  (0.117443), spread_w (0.45, 0.45, 0.10) (0.117462), early base
  0.70 (0.117503), sqrt-variance spread blend (0.117523), late
  slope 0.25 (0.117524). Dead: per-band caps, W5 lo-only, w1 0,
  mid_w 0.65/0.70, early y1z0mix, late med4 (all fail LOO).
- **Phase B — mv as a median source.** Mid med3c = median(Y5,
  clone, mv): **0.117105 / 0.208184 / 0.209462**, median-only
  0.180454. Mid 0.5·Y5+0.5·mv: 0.117244 / 0.207729 / 0.209029
  (best LOO margins). Dead: early mv (0.119890), early y1mv,
  late med4mv, late med5 (all fail LOO); mv in the spread blend
  hurts LOO in every routing tried.
- **Phase C — combine the winners.** 25 cells pass. med3c+range2:
  0.117036 / 0.207996 / 0.209037. med3c+range2+mid −0.16:
  0.116751 / 0.208094 / 0.209225. med3c+range2+spread
  (0.45,0.45,0.10): 0.117020 / 0.207827 / 0.208885 (best LOO
  margins at that point).
- **Phase D — refinement + range-vs-scale check.** mid −0.18 +
  spread (0.45,0.45,0.10): 0.116505 / 0.208049 / 0.209289.
  Pure-scale ×1.3: 0.117032; range/2.5: 0.117107 / 0.207718 /
  0.208890; range/1.5 fails LOO-large. Early ramp 0.055:
  0.116884 on the anchor.
- **Phase E — final micro-round.** mid −0.18 + spread
  (0.45,0.45,0.10) + early 0.055: **0.116353 / 0.208240 /
  0.209473, median-only 0.179039** — packaged. The more
  LOO-robust alternate (mid −0.17 + spread + early 0.055):
  0.116440 / 0.208166 / 0.209353.

## The packaged cell (v33)

Members y5/z0/y1/clone/mv (mv = per-hour median across OLS,
console shape, all shape, clone; spread = mean of the four views'
spreads). Shift (0.60 + 0.20·size), late base 0.45; size shift
−0.04; early ramp +0.055, mid ramp −0.18, late ramp +0.05; early
median Y1, mid median median(Y5, clone, mv), late median
median(Y5, Z0, Y1); k_dis 0.85; spread 0.975×(0.45·sigY5 +
0.45·sigZ0 + 0.10·sigClone); late spread ramp 0.975→0.84
(REPLACES the 0.975 base factor in the late band); W5 spread *=
1 + 1.25·size·sd; W1 sig_lo *= 1 + 0.05·size; asym
(0.90,1.10)/(1.10,0.80)/(0.95,0.85); cap 0.315; disagreement sd
= (max − min)/2 of the five member medians.

## Results (direct-graded, exact RULES scorer, "Average viewers")

- Practice **0.116353** (p00 0.1059 / p01 0.1357 / p02 0.0775 /
  p03 0.1463) vs v32 0.117532 — **−0.001179**. All four practice
  cases improve individually (v32: 0.1070/0.1366/0.0778/0.1488).
- LOO-large **0.208240** vs 0.208604 — **−0.000364**.
- LOO-mega **0.209473** vs 0.210005 — **−0.000532**.
- Median-only practice **0.179039** (v32: 0.182041) — the ninth
  round in a row the median-only loss moves. Median-only LOO-large
  0.314632 (v32: 0.316268).

Two direct grades were used (the v32-cell reproduction, then this
candidate — identical to its cached-row score to six decimals).
Validation (work/h33.py): all five members reproduce the cached
rows to **0.000000** max |diff| on practice (92 rows) AND on LOO
(851 rows), mu and sig. The combine then scores 0.116353 /
0.208240 / 0.209473 / median-only 0.179039 from the rows, and the
direct grade reproduces 0.116353 exactly.

## Honest caveats

1. **Selection bias.** Every parameter was chosen on the same four
   practice launches and the same 37-launch LOO that all previous
   rounds selected on. The margins (−0.001179 practice, −0.0004/−0.0005
   LOO) are thin, and part of them is selection noise. The practice
   board is only a proxy: placements are held until finalize, and the
   final is graded on GTA VI's real release day.
2. **The practice board is flooded.** As of 17 Nov ~21:08 WIB, 198
   entries score below 0.122315 — the top is dominated by replay /
   memorization entries (67 entries at exactly 1e-06). v33 at 0.116353
   would sit around rank #199 on the practice board, nowhere near the
   top 20, let alone the top 10 the owner now requires for submission.
   The final rewards honest generalization (v33's LOO), not memorized
   practice curves — but whether prizes follow the practice rank or the
   final rank is decided by the platform's finalize step, and that is
   not in this repo's control.
3. **Extrapolation.** The final launch (lnV1 ~ 14+) sits beyond the
   library's largest launch (13.68). The size shift (−0.04·size) and
   the size-indexed disagreement shift (0.60 + 0.20·size) both
   extrapolate hardest exactly there, and no gate in this repo can
   validate them at that size. The mid ramp is now −0.18, the
   steepest ever used — it buys practice through the two biggest
   practice launches while p01 (the most LOO-typical case) is the
   first to worsen when the ramp goes deeper (−0.19/−0.20).
4. **The range measure is mostly a rescale.** The D3 check shows a
   pure population-sd rescale (×1.3) reaches 0.117032, close to the
   range/2 anchor (0.117036) — the range profile adds only a small
   LOO edge on top of the rescale. Treat the range lever as a
   slightly-better rescale, not a new mechanism.
5. **0.100 is practically unreachable honestly.** The median is the
   bottleneck (best median-only practice ever is now 0.179039), the
   per-hour OLS residual sd is 0.35–0.45 log, and the practice answers
   are 2–2.5× rougher than the comparables. 0.100 would need a near-
   perfect smooth model. v33 gains −0.001179, the biggest single-round
   practice gain since v32 — but it comes from a fifth member and a
   disagreement redefinition, not from a new median family; no median
   family remains unbuilt or unmeasured after 33 rounds except ideas
   that need data the case file does not have.
6. **Not submitted.** The owner's rule (17 Nov ~21:21 WIB) is to submit
   only a candidate that reaches the top 10. v33 does not. It is a
   verified fallback, better than v32 on all three gates.
