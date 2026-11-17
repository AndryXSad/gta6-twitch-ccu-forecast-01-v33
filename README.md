# GTA VI launch-day Twitch CCU forecaster — v33

[![make.fun · gta6-twitch-ccu-forecast-01](https://make.fun/badge/gta6-twitch-ccu-forecast-01.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![rank](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v33/rank.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01/leaderboard?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)
[![reward](https://make.fun/badge/gta6-twitch-ccu-forecast-01/AndryXSad/gta6-twitch-ccu-forecast-01-v33/reward.svg)](https://make.fun/bounty/gta6-twitch-ccu-forecast-01?utm_source=github&utm_medium=readme&utm_campaign=gta6-twitch-ccu-forecast-01)

A make.fun `gta6-twitch-ccu-forecast-01` entry. It forecasts, for each of
hours 2–24 of Grand Theft Auto VI's release day (2026-11-19 UTC), five
quantiles (5%, 25%, 50%, 75%, 95%) of the concurrent Twitch viewers in the
"Grand Theft Auto VI" category, from the first hour's print and the
launch-day curves of past big game launches shipped inside each case.

## How it works

`main.py` (Python 3.13, standard library only, no network) builds five
member predictors in log space from the 37 comparable past launches:

1. **Y5** — v4-family blend, tier EA filter (drop tier == ea), shape
   regressions on `[lnV1, g12]`. `g12` is the "unlock momentum" feature:
   lnV1 minus log viewers 12 hours before release day.
2. **Z0** — v4-family blend, naive EA filter (drop launches whose unlock
   note contains "early access"), shape regressions on `[lnV1, g12]`.
3. **Y1** — v4-family blend, naive EA filter, shape regressions on
   `[lnV1, lnPREM]`.
4. **Clone** — v9-clone baseline: whole-library shape views, sample-sd
   spreads, weighted-mean blend (0.6 OLS / 0.2 console-shape / 0.2
   all-shape) with spread 0.81× the weighted mean of the view spreads.
5. **MV** — per-hour median across four views (per-hour OLS, console
   shape, all shape, clone forecast); spread = mean of the four views'
   spreads. New in v33: it is the mid-band median source together with
   Y5 and clone, and the fifth member in the disagreement range.

A v4-family member mixes three views per hour (0.5 OLS, 0.2
console-only shape, 0.3 all-launch shape): the median is the weighted
mean, the spread is the **sqrt-variance blend**
`sqrt(sum(w * sig^2) / sum(w))`. The shape views differ per member
(fit_hours 10 vs 20, min_fits 6 vs 5, windowed vs whole-library
pools, population vs sample sd, OLS pool counted before vs after the
hour filter). All five members reproduce the cached rows to 0.000000
on practice AND on LOO.

The **log-median is chosen per hour band** — the median is the bottleneck:

| band | hours | median source | ramp |
| ---- | ----- | ------------- | ---- |
| early | 2–6 | Y1 log-median | 0.00 → +0.055 |
| mid | 7–14 | median of Y5 / clone / MV log-medians | 0.00 → −0.18 |
| late | 15–24 | median of Y5 / Z0 / Y1 log-medians | 0.00 → +0.05 |

On top of that, a **size-dependent shift**: above lnV1 11.5 the median is
pulled down by 0.04 per log unit of size, in every band, and a **deeper,
size-indexed disagreement shift** `mu += (0.60 + 0.20 * size) * sd` of the
five member log-medians (late band base 0.45). The disagreement measure
is new in v33: **(max − min) / 2** of the member medians — a range
measure, not a population sd. Hours where one member is an outlier get
a bigger shift, which the LOO gate rewards. A pure rescale of the
population sd (×1.3) scores 0.117032, worse than the range measure, so
the profile — not just the scale — is doing work.

The **log-spread** is 0.975× the weighted mean of the member spreads
(45% Y5 + 45% Z0 + 10% clone), widened by a **disagreement factor**
`1 + 0.85 * sd(member medians)`, shaped per band with separate
lower/upper scales: early (0.90, 1.10), mid (1.10, 0.80), late
(0.95, 0.85), narrowed across the late band (spread scale 0.975 → 0.84
from hour 15 to hour 24 — this REPLACES the 0.975 base factor in the
late band), and finally **capped at 0.315**. Two size levers from v31:
**W1** grows only the lower spread with size (`slo *= 1 + 0.05 * size`)
and **W5** lets size amplify the disagreement effect on both sides
(`*= 1 + 1.25 * size * sd`). Quantiles are lognormal.

## Results

Direct-graded with the exact RULES scorer on the four practice cases
("Average viewers" answers): **practice pinball loss 0.116353**
(p00 0.1059 / p01 0.1357 / p02 0.0775 / p03 0.1463),
LOO-large **0.208240**, LOO-mega **0.209473**. That beats v32
(0.117532 / 0.208604 / 0.210005) on practice by −0.001179, with
LOO-large −0.000364 and LOO-mega −0.000532 — and all four practice
cases improve individually. Median-only practice is **0.179039**
(LOO-large 0.314632) — the median moved again, the ninth round in a
row. See `DESIGN.md` for the full reasoning record and the honest
caveats. Not submitted: the owner's rule is to submit only a
candidate that reaches the top 10.

## Run it

```sh
MAKEFUN_INPUT_DIR=<dir with *.in> MAKEFUN_OUTPUT_DIR=<dir> python3 main.py
```

Each `*.in` case produces a matching `*.out`:
`{"hours": [{"hour": h, "q": [q05, q25, q50, q75, q95]}]}` for h = 2..24.
