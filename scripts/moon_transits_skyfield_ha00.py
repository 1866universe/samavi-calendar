"""
Moon HA=00 Upper-Meridian Transit Generator (Skyfield / DE440)
Historical Julian/Gregorian Calendar Display — Fully Audited & Validated

Calculates precise UTC/UT1 moments of Moon Upper Meridian Transit (HA = 0h)
for the Samavi Calendar Dataset and global astronomical research.
"""

import os
import sys
import numpy as np
import pandas as pd
import skyfield_data
from skyfield.api import wgs84
from skyfield.iokit import Loader
from skyfield.jpllib import SpiceKernel
from scipy.optimize import brentq

# ==============================================================================
# Configuration & Observation Parameters
# ==============================================================================
LAT_DEG = 35 + 39 / 60 + 32 / 3600      # 35° 39' 32" N (Tehran Meridian)
LON_DEG = 51 + 3 / 60 + 28 / 3600       # 51° 3' 28" E
ELEV_M = 0

START_YEAR  = 1552
START_MONTH = 3
START_DAY   = 24
START_HOUR  = 0

TEST_YEARS  = 660                       # Total years to process

SCAN_STEP_HOURS = 4.0                   # Robust uniform scan window
ROOT_XTOL = 1.0e-11                     # Solver precision in days (~0.86 microseconds)

# ==============================================================================
# Historic Calendar Utilities (Julian / Gregorian Switch on 1582-10-15)
# ==============================================================================

def cal_to_jd(year, month, day, hour=0, minute=0, second=0.0):
    """Convert historic calendar date to Julian Date (JD)."""
    y = year
    m = month
    d = day

    if m <= 2:
        y -= 1
        m += 12

    is_gregorian = (year > 1582) or (year == 1582 and (month > 10 or (month == 10 and day >= 15)))

    if is_gregorian:
        a = y // 100
        b = 2 - a + (a // 4)
    else:
        b = 0

    day_fraction = (hour + minute / 60.0 + second / 3600.0) / 24.0
    jd = int(365.25 * (y + 4716)) + int(30.6001 * (m + 1)) + d + b - 1524.5 + day_fraction
    return jd


def jd_to_cal(jd):
    """Convert Julian Date (JD) to historic calendar date (Julian before 1582-10-15; Gregorian after)."""
    z = int(jd + 0.5)
    f = (jd + 0.5) - z

    if z < 2299161:
        a = z
    else:
        alpha = int((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - (alpha // 4)

    b = a + 1524
    c = int((b - 122.1) / 365.25)
    d = int(365.25 * c)
    e = int((b - d) / 30.6001)

    day_value = b - d - int(30.6001 * e) + f
    month = e - 1 if e < 14 else e - 13
    year = c - 4716 if month > 2 else c - 4715

    day = int(day_value)
    seconds_of_day = int(round((day_value - day) * 86400.0))
    if seconds_of_day >= 86400:
        return jd_to_cal(jd + 1e-9)

    hour, rem = divmod(seconds_of_day, 3600)
    minute, second = divmod(rem, 60)

    return year, month, day, hour, minute, second


# ==============================================================================
# Celestial Mathematics & HA=00 Solver
# ==============================================================================

def actual_hour_angle_hours(observer, target, t):
    """Returns actual apparent Hour Angle in hours [0, 24)."""
    apparent = observer.at(t).observe(target).apparent()
    ha_hours = apparent.hadec()[0].hours % 24.0
    return ha_hours


def signed_upper_transit_offset_hours(observer, target, t):
    """
    Returns signed angular distance in hours from Upper Meridian (HA = 0h).
    Domain: (-12.0, +12.0] hours.
    Upper transit occurs when this value crosses zero from negative to positive.
    """
    ha_hours = observer.at(t).observe(target).apparent().hadec()[0].hours % 24.0
    if ha_hours > 12.0:
        return ha_hours - 24.0
    return ha_hours


def is_real_upper_transit_bracket(val_left, val_right):
    """
    Validates genuine crossing of HA = 0h.
    Requires crossing from negative to positive within smooth progression.
    """
    return (val_left < 0.0 < val_right) and (val_right - val_left < 8.0)


def refine_upper_transit(ts, observer, target, jd_left, jd_right, tol=ROOT_XTOL):
    """Refine Upper Transit root using Brent's method."""
    def function(jd_val):
        t_eval = ts.ut1_jd(jd_val)
        return signed_upper_transit_offset_hours(observer, target, t_eval)

    try:
        f_left = function(jd_left)
        f_right = function(jd_right)

        if not (f_left < 0.0 < f_right):
            return None, f"invalid_bracket_signs: left={f_left:.6e}, right={f_right:.6e}"

        root_jd = brentq(function, jd_left, jd_right, xtol=tol)
        root_val = function(root_jd)

        if abs(root_val) > 1.0e-5:
            return None, f"root_residual_too_large: {root_val:.6e}"

        return root_jd, None

    except Exception as exc:
        return None, f"exception: {str(exc)}"


# ==============================================================================
# Scan Engine & Validation
# ==============================================================================

def scan_all_upper_transits(ts, observer, target, jd_start, jd_end, step_hours=SCAN_STEP_HOURS):
    """Uniformly scan the timeline for all Upper Meridian Transits."""
    step_days = step_hours / 24.0
    total_steps = int(np.ceil((jd_end - jd_start) / step_days))

    raw_roots = []
    solver_issues = []

    cur_jd = jd_start
    cur_t = ts.ut1_jd(cur_jd)
    cur_val = signed_upper_transit_offset_hours(observer, target, cur_t)

    print(f"Processing interval:")
    print(f"  Start: {START_YEAR:04d}-{START_MONTH:02d}-{START_DAY:02d} {START_HOUR:02d}:00:00")
    print(f"  End:   {START_YEAR + TEST_YEARS:04d}-{START_MONTH:02d}-{START_DAY:02d} {START_HOUR:02d}:00:00")
    print(f"  Start JD: {jd_start:.10f}")
    print(f"  End JD:   {jd_end:.10f}")
    print(f"  Scan step: {step_hours:.3f} hours")

    for step_index in range(1, total_steps + 1):
        next_jd = min(jd_start + step_index * step_days, jd_end)
        next_t = ts.ut1_jd(next_jd)
        next_val = signed_upper_transit_offset_hours(observer, target, next_t)

        if is_real_upper_transit_bracket(cur_val, next_val):
            root_jd, failure_reason = refine_upper_transit(ts, observer, target, cur_jd, next_jd)

            if root_jd is not None:
                raw_roots.append(root_jd)
            else:
                solver_issues.append({
                    "step_index": step_index,
                    "jd_left": cur_jd,
                    "jd_right": next_jd,
                    "val_left": cur_val,
                    "val_right": next_val,
                    "reason": failure_reason,
                })

        if step_index % 5000 == 0:
            progress = (step_index / total_steps) * 100.0
            print(f"Progress: {progress:6.2f}% | JD {next_jd:.5f} | roots found: {len(raw_roots):,}")

        cur_jd, cur_val = next_jd, next_val

    return raw_roots, solver_issues


def deduplicate_transit_roots(roots_list, merge_threshold_days=0.25):
    """Remove duplicate root findings within physical threshold."""
    if not roots_list:
        return [], []

    sorted_roots = sorted(roots_list)
    unique_roots = [sorted_roots[0]]
    duplicates_removed = []

    for item in sorted_roots[1:]:
        if (item - unique_roots[-1]) < merge_threshold_days:
            duplicates_removed.append(item)
        else:
            unique_roots.append(item)

    return unique_roots, duplicates_removed


def build_audit_dataframe(unique_roots):
    """Construct interval audit dataframe."""
    records = []
    for i in range(len(unique_roots)):
        curr = unique_roots[i]
        prev = unique_roots[i - 1] if i > 0 else np.nan
        next_ = unique_roots[i + 1] if i < len(unique_roots) - 1 else np.nan

        dt_prev_h = (curr - prev) * 24.0 if not np.isnan(prev) else np.nan
        dt_next_h = (next_ - curr) * 24.0 if not np.isnan(next_) else np.nan

        records.append({
            "transit_id": i + 1,
            "jd": curr,
            "interval_from_prev_hours": dt_prev_h,
            "interval_to_next_hours": dt_next_h,
        })
    return pd.DataFrame(records)


def validate_transit_sequence(audit_df, total_days):
    """Perform mathematical and physical anomaly checks."""
    total_count = len(audit_df)
    if total_count < 2:
        return {"status": "FAIL", "reason": "Insufficient root count."}

    diffs = audit_df["interval_to_next_hours"].dropna()

    is_monotonic = bool((audit_df["jd"].diff().dropna() > 0).all())
    short_intervals = int((diffs < 23.5).sum())
    long_intervals = int((diffs > 26.5).sum())
    severe_gaps = int((diffs > 36.0).sum())

    expected_count = int(round(total_days / (24.8412 / 24.0)))
    count_difference = total_count - expected_count

    stats = {
        "total_roots": total_count,
        "strictly_increasing": is_monotonic,
        "short_interval_count": short_intervals,
        "long_interval_count": long_intervals,
        "severe_gap_count": severe_gaps,
        "minimum_interval_hours": float(diffs.min()),
        "maximum_interval_hours": float(diffs.max()),
        "mean_interval_hours": float(diffs.mean()),
        "median_interval_hours": float(diffs.median()),
        "estimated_expected_count": expected_count,
        "count_difference_from_estimate": count_difference,
    }

    if (not is_monotonic) or (severe_gaps > 0) or (abs(count_difference) > 10):
        stats["status"] = "FAIL"
    else:
        stats["status"] = "PASS"

    return stats


# ==============================================================================
# Main Program Entry Point
# ==============================================================================

def main():
    print("=" * 78)
    print("Moon HA=00 Upper-Meridian Transit Generator")
    print("Historical Julian/Gregorian calendar — Stellarium-compatible display")
    print("=" * 78)
    print()

    data_dir = os.path.join(os.path.dirname(skyfield_data.__file__), "data")
    print("Skyfield data directory:\n" + data_dir)
    print()

    load = Loader(data_dir)
    ts = load.timescale()

    de440_path = os.path.join(data_dir, "de440.bsp")
    print("Ephemeris:\n" + de440_path)
    print()

    eph = SpiceKernel(de440_path)
    earth = eph["earth"]
    moon = eph["moon"]
    observer = earth + wgs84.latlon(LAT_DEG, LON_DEG, elevation_m=ELEV_M)

    jd_start = cal_to_jd(START_YEAR, START_MONTH, START_DAY, START_HOUR, 0, 0.0)
    jd_end = cal_to_jd(START_YEAR + TEST_YEARS, START_MONTH, START_DAY, START_HOUR, 0, 0.0)

    raw_roots, solver_issues = scan_all_upper_transits(
        ts, observer, moon, jd_start, jd_end, step_hours=SCAN_STEP_HOURS
    )

    unique_roots, duplicates_removed = deduplicate_transit_roots(raw_roots)
    audit_df = build_audit_dataframe(unique_roots)
    validation_summary = validate_transit_sequence(audit_df, jd_end - jd_start)

    main_records = []
    for i, root_jd in enumerate(unique_roots, start=1):
        y, m, d, hh, mm, ss = jd_to_cal(root_jd)
        transit_str = f"{y:04d}-{m:02d}-{d:02d} {hh:02d}:{mm:02d}:{ss:02d}"

        t_eval = ts.ut1_jd(root_jd)
        ha_actual = actual_hour_angle_hours(observer, moon, t_eval)
        offset_val = signed_upper_transit_offset_hours(observer, moon, t_eval)

        main_records.append({
            "transit_id": i,
            "stellarium_transit_utc": transit_str,
            "jd": f"{root_jd:.10f}",
            "ha_actual": f"{ha_actual:.12f}",
            "ha00_residual_hours": f"{offset_val:.12e}",
        })

    main_df = pd.DataFrame(main_records)

    base_name = f"moon_transits_HA00_{START_YEAR}"
    script_dir = os.path.dirname(os.path.abspath(__file__))

    main_csv = os.path.join(script_dir, f"{base_name}_Stellarium.csv")
    audit_csv = os.path.join(script_dir, f"{base_name}_interval_audit.csv")
    issues_csv = os.path.join(script_dir, f"{base_name}_solver_issues.csv")
    dups_csv = os.path.join(script_dir, f"{base_name}_duplicates_removed.csv")
    summary_txt = os.path.join(script_dir, f"{base_name}_run_summary.txt")

    main_df.to_csv(main_csv, index=False)
    audit_df.to_csv(audit_csv, index=False)
    pd.DataFrame(solver_issues).to_csv(issues_csv, index=False)
    pd.DataFrame({"jd": duplicates_removed}).to_csv(dups_csv, index=False)

    with open(summary_txt, "w", encoding="utf-8") as f:
        f.write("Moon Upper-Meridian Transit (HA=00) Dataset — Run Summary\n")
        f.write("=" * 64 + "\n\n")
        f.write(f"Start civil date: {START_YEAR:04d}-{START_MONTH:02d}-{START_DAY:02d} {START_HOUR:02d}:00:00\n")
        f.write(f"End civil date:   {START_YEAR + TEST_YEARS:04d}-{START_MONTH:02d}-{START_DAY:02d} {START_HOUR:02d}:00:00\n")
        f.write(f"Detected unique transits: {len(unique_roots):,}\n")
        f.write(f"Solver issues: {len(solver_issues):,}\n")
        f.write(f"Duplicate detections removed: {len(duplicates_removed):,}\n\n")
        f.write("Validation Metrics\n")
        f.write("-" * 64 + "\n")
        for k, v in validation_summary.items():
            f.write(f"{k}: {v}\n")

    print()
    print("=" * 78)
    print("CALCULATION COMPLETED")
    print("=" * 78)
    print(f"Raw root detections:          {len(raw_roots):,}")
    print(f"Unique transit records:       {len(unique_roots):,}")
    print(f"Duplicate detections removed: {len(duplicates_removed):,}")
    print(f"Solver issues:                {len(solver_issues):,}")
    print(f"Minimum interval: {validation_summary['minimum_interval_hours']:.9f} hours")
    print(f"Maximum interval: {validation_summary['maximum_interval_hours']:.9f} hours")
    print(f"Mean interval:    {validation_summary['mean_interval_hours']:.9f} hours")
    print(f"Median interval:  {validation_summary['median_interval_hours']:.9f} hours")
    print()
    print("Files saved:")
    print(f"  Main dataset:     {main_csv}")
    print(f"  Interval audit:   {audit_csv}")
    print(f"  Solver issues:    {issues_csv}")
    print(f"  Run summary:      {summary_txt}")
    print()

    if len(main_df) > 0:
        print("First five generated records:")
        print(main_df.head(5).to_string(index=False))
        print("\nLast five generated records:")
        print(main_df.tail(5).to_string(index=False))

    print()
    print(f"VALIDATION RESULT: {validation_summary['status']}")
    if validation_summary['status'] == "PASS":
        print("No skipped-transit pattern was detected.")
    else:
        print("WARNING: Sequence anomalies detected.")


if __name__ == "__main__":
    main()
