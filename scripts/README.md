# Astronomical Transit Scripts

Python scripts that generate the astronomical backbone of the Samavi Calendar:
precise Moon meridian transit moments computed with the DE440 planetary
ephemeris via [Skyfield](https://rhodesmill.org/skyfield/).

| Script | Reference Event | Description |
|---|---|---|
| `moon_transits_skyfield_ha00.py` | Upper transit (HA = 00h) | Moon crossing the local meridian (upper culmination) |
| `moon_transits_skyfield_ha12.py` | Lower transit (HA = 12h) | Moon crossing the anti-meridian (lower culmination) |

## Key design points

- **Uniform 4-hour scan step** with bracket validation — no skipped transits,
  no heuristic time jumps.
- **Brent's method root refinement** with a residual tolerance of `1.0e-5` hours.
- **Historic calendar output** (Julian before 1582-10-15, Gregorian after),
  compatible with Stellarium's historic calendar display.
- **Automatic ephemeris handling:** on first run, the script looks for
  `de440.bsp` in the local Skyfield data directory and **downloads it
  automatically if it is not already present**. No manual download is required.
- **Built-in validation engine:** every run produces an interval audit,
  a solver-issue log, and a run summary with statistical checks
  (monotonic sequence, gap detection, expected-count comparison).

## Validation reports (`validation/`)

Run summaries from verification passes over three sensitive historical
windows — 1552 (pre-Gregorian epoch), 1582 (Julian→Gregorian transition),
and 2000 (modern J2000 era). All windows: `PASS`, 0 solver issues,
0 skipped transits.
`