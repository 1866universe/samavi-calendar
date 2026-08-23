"""
Moon HA=12 Upper-Meridian Transit Generator (Skyfield / DE440)
Historical Julian/Gregorian Calendar Display — Fully Audited & Validated

Calculates precise UTC/UT1 moments of Moon Upper Meridian Transit (HA = 0h)
for the Samavi Calendar Dataset and global astronomical research.
"""

import os
import math
import traceback
from datetime import datetime

import pandas as pd
import skyfield_data
from scipy.optimize import brentq
from skyfield.api import wgs84
from skyfield.iokit import Loader
from skyfield.jpllib import SpiceKernel


# =============================================================================
# SETTINGS
# =============================================================================

# Observer location:
# 35° 39' 32" N
# 51° 03' 28" E
LAT_DEG = 35.0 + 39.0 / 60.0 + 32.0 / 3600.0
LON_DEG = 51.0 + 3.0 / 60.0 + 28.0 / 3600.0
ELEV_M = 0.0

# Historical civil-calendar start date:
# Julian calendar before 1582-10-15
# Gregorian calendar from 1582-10-15 onward
START_YEAR = 1552
START_MONTH = 3
START_DAY = 24
START_HOUR = 0
START_MINUTE = 0
START_SECOND = 0.0

# Calculation span
TEST_YEARS = 660

# Coarse scan interval.
#
# Six hours is comfortably shorter than half of a lunar transit cycle.
# It provides a robust sign-change bracket around every lower-meridian transit
# without relying on the previously found transit.
SCAN_STEP_HOURS = 4.0

# Brent root-finder tolerance in Julian days.
# 1e-11 day is much smaller than one displayed second.
ROOT_XTOL_DAYS = 1.0e-11

# Roots closer than this threshold are considered duplicate detections.
DUPLICATE_TOLERANCE_SECONDS = 2.0

# Expected interval range between consecutive lower-meridian lunar transits.
#
# This range is intentionally broad. It is used for validation and reporting,
# not for event detection.
MIN_NORMAL_INTERVAL_HOURS = 20.0
MAX_NORMAL_INTERVAL_HOURS = 30.0

# A second category indicating a severe skipped-transit pattern.
SEVERE_GAP_HOURS = 36.0

# Hour-angle distance limit used to reject the false discontinuity at +/-12 h.
#
# Around a genuine lower-meridian transit, the normalized target function
# crosses zero and both endpoint values are near zero.
#
# At the unrelated wrap boundary, values jump between approximately +12 and
# -12. Requiring both endpoints to remain inside this limit rejects that wrap.
REAL_ROOT_ENDPOINT_LIMIT_HOURS = 6.0

# Progress message interval
PROGRESS_EVERY_DAYS = 365.25

# Output numeric precision
JD_DECIMAL_PLACES = 10
HA_DECIMAL_PLACES = 12


# =============================================================================
# HISTORICAL CALENDAR CONVERSION
# =============================================================================

def is_gregorian_civil_date(year, month, day):
    """
    Return True if the supplied historical civil date belongs to the
    Gregorian calendar according to the Stellarium-compatible transition:

        Julian through 1582-10-04
        Gregorian from 1582-10-15

    The omitted dates 1582-10-05 through 1582-10-14 are invalid.
    """
    if (year, month, day) >= (1582, 10, 15):
        return True

    return False


def validate_historical_date(year, month, day):
    """
    Reject the ten civil dates omitted during the Gregorian reform.
    """
    if year == 1582 and month == 10 and 5 <= day <= 14:
        raise ValueError(
            "Historical civil dates 1582-10-05 through 1582-10-14 "
            "do not exist in the selected Julian/Gregorian calendar."
        )


def cal_to_jd(year, month, day, hour=0, minute=0, second=0.0):
    """
    Convert a historical civil-calendar date to Julian Day.

    Calendar rule:
      - Julian calendar before 1582-10-15
      - Gregorian calendar from 1582-10-15 onward

    This convention is intended to remain synchronized with Stellarium's
    historical civil-calendar display.
    """
    validate_historical_date(year, month, day)

    original_year = year
    original_month = month
    original_day = day

    if month <= 2:
        year -= 1
        month += 12

    if is_gregorian_civil_date(
        original_year,
        original_month,
        original_day
    ):
        a = math.floor(year / 100.0)
        b = 2 - a + math.floor(a / 4.0)
    else:
        b = 0

    day_fraction = (
        float(hour)
        + float(minute) / 60.0
        + float(second) / 3600.0
    ) / 24.0

    jd = (
        math.floor(365.25 * (year + 4716))
        + math.floor(30.6001 * (month + 1))
        + original_day
        + b
        - 1524.5
        + day_fraction
    )

    return float(jd)


def jd_to_calendar_components(jd):
    """
    Convert Julian Day to historical civil-calendar components without
    rounding the time to a whole second.

    Returns:
        year, month, day, fractional_day
    """
    shifted = float(jd) + 0.5
    z = math.floor(shifted)
    f = shifted - z

    if z < 2299161:
        a = z
    else:
        alpha = math.floor((z - 1867216.25) / 36524.25)
        a = z + 1 + alpha - math.floor(alpha / 4.0)

    b = a + 1524
    c = math.floor((b - 122.1) / 365.25)
    d = math.floor(365.25 * c)
    e = math.floor((b - d) / 30.6001)

    day_value = b - d - math.floor(30.6001 * e) + f

    if e < 14:
        month = e - 1
    else:
        month = e - 13

    if month > 2:
        year = c - 4716
    else:
        year = c - 4715

    day = math.floor(day_value)
    fractional_day = day_value - day

    return int(year), int(month), int(day), float(fractional_day)


def jd_to_cal(jd):
    """
    Convert Julian Day to the selected historical civil calendar.

    The returned time is rounded to the nearest displayed second.
    Midnight rollover is handled iteratively and does not use recursion.
    """
    working_jd = float(jd)

    for _ in range(3):
        year, month, day, fractional_day = (
            jd_to_calendar_components(working_jd)
        )

        seconds_of_day = int(round(fractional_day * 86400.0))

        if seconds_of_day < 86400:
            hour, remainder = divmod(seconds_of_day, 3600)
            minute, second = divmod(remainder, 60)

            return (
                year,
                month,
                day,
                int(hour),
                int(minute),
                int(second),
            )

        # Rounding reached the next midnight. Convert the following instant
        # instead of recursively calling this function.
        working_jd = math.floor(working_jd + 0.5) + 0.5

    raise RuntimeError(
        f"Could not resolve civil-calendar midnight rollover for JD {jd!r}"
    )


def format_historical_datetime(jd):
    """
    Format a Julian Day in the selected historical civil calendar.
    """
    year, month, day, hour, minute, second = jd_to_cal(jd)

    return (
        f"{year:04d}-{month:02d}-{day:02d} "
        f"{hour:02d}:{minute:02d}:{second:02d}"
    )


# =============================================================================
# HOUR ANGLE
# =============================================================================

def normalize_ha_hours(ha_h):
    """
    Normalize an hour angle to [-12, +12) hours.
    """
    return ((float(ha_h) + 12.0) % 24.0) - 12.0


def signed_lower_transit_offset_hours(observer, target, t):
    """
    Returns the signed angular distance in hours from the lower meridian (HA = 12h).
    Domain: (-12.0, +12.0] hours.
    A genuine lower transit occurs when this value crosses zero from negative to positive.
    """
    ha_hours = observer.at(t).observe(target).apparent().hadec()[0].hours % 24.0

    offset = ha_hours - 12.0
    if offset <= -12.0:
        offset += 24.0
    elif offset > 12.0:
        offset -= 24.0
    return offset


def actual_hour_angle_hours(observer, target, t):
    """
    Return the apparent hour angle normalized to [0, 24) hours.

    At a correctly solved lower-meridian transit, this value should be close
    to 12 hours.
    """
    apparent = observer.at(t).observe(target).apparent()
    return apparent.hadec()[0].hours % 24.0


# =============================================================================
# ROOT DETECTION
# =============================================================================

def is_real_lower_transit_bracket(left_value, right_value):
    """
    A true lower transit occurs when Moon's HA increases continuously through 12h:
    left_value <= 0 and right_value >= 0, and right_value > left_value.

    The wrap discontinuity at upper transit (HA=0) goes from +12 to -12 (downwards),
    which is strictly impossible for a true transit.
    """
    if not (math.isfinite(left_value) and math.isfinite(right_value)):
        return False

    # True transit crosses zero with positive slope: left <= 0 and right >= 0
    if left_value <= 0.0 and right_value >= 0.0:
        # Continuity check: over a 4h step, delta is around +3.5 to +4.5 hours.
        # It must never jump wildly.
        if 0.0 < (right_value - left_value) < 8.0:
            return True

    return False


def refine_lower_transit(ts, observer, target, left_jd, right_jd, left_value=None, right_value=None):
    """
    Refine one bracketed lower-meridian transit using Brent's method.
    """
    def function(jd):
        return signed_lower_transit_offset_hours(observer, target, ts.ut1_jd(jd))

    if left_value is None:
        left_value = function(left_jd)
    if right_value is None:
        right_value = function(right_jd)

    if not is_real_lower_transit_bracket(left_value, right_value):
        return None, "not_a_real_lower_transit_bracket"

    if abs(left_value) < 1.0e-13:
        return float(left_jd), "left_endpoint_root"
    if abs(right_value) < 1.0e-13:
        return float(right_jd), "right_endpoint_root"

    try:
        root_jd = brentq(
            function,
            left_jd,
            right_jd,
            xtol=ROOT_XTOL_DAYS,
            rtol=1.0e-14,
            maxiter=100,
        )
    except Exception as exc:
        return None, f"brentq_failure: {exc}"

    root_val = function(root_jd)
    if abs(root_val) > 1.0e-5:
        return None, f"root_residual_too_large: {root_val:.6e}"

    return float(root_jd), "ok"


def scan_all_lower_transits(
    ts,
    observer,
    target,
    start_jd,
    end_jd,
):
    """
    Scan the entire requested time span in fixed independent intervals.

    This design deliberately avoids chaining each search to the previously
    detected transit. Therefore, a single failure cannot shift the search
    window and cause a cascade of multi-day data loss.

    Returns:
        roots
        solver_issues
    """
    scan_step_days = SCAN_STEP_HOURS / 24.0

    roots = []
    solver_issues = []

    current_jd = float(start_jd)
    current_value = signed_lower_transit_offset_hours(
        observer,
        target,
        ts.ut1_jd(current_jd),
    )

    total_days = end_jd - start_jd
    next_progress_day = PROGRESS_EVERY_DAYS

    while current_jd < end_jd:
        next_jd = min(current_jd + scan_step_days, end_jd)

        try:
            next_value = signed_lower_transit_offset_hours(
                observer,
                target,
                ts.ut1_jd(next_jd),
            )
        except Exception as exc:
            solver_issues.append({
                "issue_type": "sample_evaluation_failure",
                "left_jd": current_jd,
                "right_jd": next_jd,
                "message": repr(exc),
            })

            # Restart cleanly at the next sample.
            current_jd = next_jd
            try:
                current_value = signed_lower_transit_offset_hours(
                    observer,
                    target,
                    ts.ut1_jd(current_jd),
                )
            except Exception:
                current_value = float("nan")
            continue

        if is_real_lower_transit_bracket(
            current_value,
            next_value,
        ):
            root_jd, message = refine_lower_transit(
                ts=ts,
                observer=observer,
                target=target,
                left_jd=current_jd,
                right_jd=next_jd,
                left_value=current_value,
                right_value=next_value,
            )

            if root_jd is None:
                solver_issues.append({
                    "issue_type": "root_refinement_failure",
                    "left_jd": current_jd,
                    "right_jd": next_jd,
                    "left_value_hours": current_value,
                    "right_value_hours": next_value,
                    "message": message,
                })
            elif start_jd <= root_jd < end_jd:
                roots.append(root_jd)

        current_jd = next_jd
        current_value = next_value

        elapsed_days = current_jd - start_jd
        if elapsed_days >= next_progress_day:
            percentage = 100.0 * elapsed_days / total_days

            print(
                f"Progress: {percentage:6.2f}% | "
                f"JD {current_jd:.5f} | "
                f"raw roots: {len(roots):,}"
            )

            next_progress_day += PROGRESS_EVERY_DAYS

    return roots, solver_issues


def deduplicate_roots(roots):
    """
    Sort detected roots and remove duplicate detections.

    Duplicate roots can occur only when a root lies exactly on a scan-block
    boundary. Such duplicates are expected to be extremely rare.
    """
    if not roots:
        return [], []

    tolerance_days = DUPLICATE_TOLERANCE_SECONDS / 86400.0
    sorted_roots = sorted(float(root) for root in roots)

    unique_roots = [sorted_roots[0]]
    duplicate_records = []

    for root in sorted_roots[1:]:
        previous = unique_roots[-1]
        difference_seconds = (root - previous) * 86400.0

        if abs(root - previous) <= tolerance_days:
            duplicate_records.append({
                "kept_jd": previous,
                "discarded_jd": root,
                "difference_seconds": difference_seconds,
            })
        else:
            unique_roots.append(root)

    return unique_roots, duplicate_records


# =============================================================================
# DATASET CONSTRUCTION AND VALIDATION
# =============================================================================

def build_output_rows(ts, observer, target, roots):
    """
    Build the primary CSV rows while preserving the original column names.
    """
    rows = []

    for transit_id, root_jd in enumerate(roots, start=1):
        t = ts.ut1_jd(root_jd)

        actual_ha = actual_hour_angle_hours(
            observer,
            target,
            t,
        )

        target_residual = normalize_ha_hours(
            actual_ha - 12.0
        )

        rows.append({
            "transit_id": transit_id,
            "stellarium_transit_utc":
                format_historical_datetime(root_jd),
            "jd": f"{root_jd:.{JD_DECIMAL_PLACES}f}",
            "ha_actual": f"{actual_ha:.{HA_DECIMAL_PLACES}f}",
            "ha12_residual_hours":
                f"{target_residual:.12e}",
        })

    return rows


def classify_interval(interval_hours):
    """
    Classify an interval for the independent audit report.
    """
    if interval_hours >= SEVERE_GAP_HOURS:
        return "SEVERE_GAP_POSSIBLE_SKIPPED_TRANSIT"

    if interval_hours > MAX_NORMAL_INTERVAL_HOURS:
        return "GAP_POSSIBLE_SKIPPED_TRANSIT"

    if interval_hours < MIN_NORMAL_INTERVAL_HOURS:
        return "SHORT_INTERVAL_OR_DUPLICATE"

    return "OK"


def build_interval_audit(roots):
    """
    Build an interval-by-interval validation report.
    """
    rows = []

    for index in range(1, len(roots)):
        previous_jd = roots[index - 1]
        current_jd = roots[index]

        interval_days = current_jd - previous_jd
        interval_hours = interval_days * 24.0
        interval_minutes = interval_days * 1440.0

        rows.append({
            "previous_transit_id": index,
            "current_transit_id": index + 1,
            "previous_transit":
                format_historical_datetime(previous_jd),
            "current_transit":
                format_historical_datetime(current_jd),
            "previous_jd": f"{previous_jd:.10f}",
            "current_jd": f"{current_jd:.10f}",
            "interval_minutes": f"{interval_minutes:.6f}",
            "interval_hours": f"{interval_hours:.9f}",
            "classification":
                classify_interval(interval_hours),
        })

    return rows


def validate_root_sequence(roots, start_jd, end_jd):
    """
    Perform independent structural validation of the detected sequence.
    """
    report = {
        "total_roots": len(roots),
        "strictly_increasing": True,
        "duplicate_or_reverse_count": 0,
        "short_interval_count": 0,
        "long_interval_count": 0,
        "severe_gap_count": 0,
        "minimum_interval_hours": None,
        "maximum_interval_hours": None,
        "mean_interval_hours": None,
        "median_interval_hours": None,
        "estimated_expected_count": None,
        "count_difference_from_estimate": None,
    }

    if len(roots) < 2:
        return report

    intervals = [
        (roots[i] - roots[i - 1]) * 24.0
        for i in range(1, len(roots))
    ]

    duplicate_or_reverse_count = sum(
        interval <= 0.0
        for interval in intervals
    )

    short_count = sum(
        interval < MIN_NORMAL_INTERVAL_HOURS
        for interval in intervals
    )

    long_count = sum(
        interval > MAX_NORMAL_INTERVAL_HOURS
        for interval in intervals
    )

    severe_count = sum(
        interval >= SEVERE_GAP_HOURS
        for interval in intervals
    )

    sorted_intervals = sorted(intervals)
    middle = len(sorted_intervals) // 2

    if len(sorted_intervals) % 2 == 1:
        median_interval = sorted_intervals[middle]
    else:
        median_interval = (
            sorted_intervals[middle - 1]
            + sorted_intervals[middle]
        ) / 2.0

    mean_interval = sum(intervals) / len(intervals)

    # Count estimate is informational only. It is not used to create or reject
    # any event.
    estimated_count = int(round(
        (end_jd - start_jd) * 24.0 / mean_interval
    ))

    report.update({
        "strictly_increasing":
            duplicate_or_reverse_count == 0,
        "duplicate_or_reverse_count":
            duplicate_or_reverse_count,
        "short_interval_count":
            short_count,
        "long_interval_count":
            long_count,
        "severe_gap_count":
            severe_count,
        "minimum_interval_hours":
            min(intervals),
        "maximum_interval_hours":
            max(intervals),
        "mean_interval_hours":
            mean_interval,
        "median_interval_hours":
            median_interval,
        "estimated_expected_count":
            estimated_count,
        "count_difference_from_estimate":
            len(roots) - estimated_count,
    })

    return report


def write_summary_file(
    path,
    start_jd,
    end_jd,
    roots,
    solver_issues,
    duplicate_records,
    validation,
):
    """
    Write a human-readable run summary.
    """
    with open(path, "w", encoding="utf-8") as handle:
        handle.write(
            "Moon Lower-Meridian Transit Dataset — Run Summary\n"
        )
        handle.write("=" * 64 + "\n\n")

        handle.write(
            f"Generated at: {datetime.now().isoformat(timespec='seconds')}\n"
        )
        handle.write(
            f"Start civil date: "
            f"{format_historical_datetime(start_jd)}\n"
        )
        handle.write(
            f"End civil date:   "
            f"{format_historical_datetime(end_jd)}\n"
        )
        handle.write(f"Start JD: {start_jd:.10f}\n")
        handle.write(f"End JD:   {end_jd:.10f}\n")
        handle.write(
            f"Scan step: {SCAN_STEP_HOURS:.3f} hours\n"
        )
        handle.write(
            f"Root tolerance: {ROOT_XTOL_DAYS:.3e} day\n"
        )
        handle.write(
            f"Detected unique transits: {len(roots):,}\n"
        )
        handle.write(
            f"Solver issues: {len(solver_issues):,}\n"
        )
        handle.write(
            f"Duplicate detections removed: "
            f"{len(duplicate_records):,}\n\n"
        )

        handle.write("Validation\n")
        handle.write("-" * 64 + "\n")

        for key, value in validation.items():
            handle.write(f"{key}: {value}\n")

        handle.write("\nInterpretation\n")
        handle.write("-" * 64 + "\n")

        if (
            validation["strictly_increasing"]
            and validation["short_interval_count"] == 0
            and validation["long_interval_count"] == 0
            and len(solver_issues) == 0
        ):
            handle.write(
                "PASS: No structural transit-sequence anomaly was detected.\n"
            )
        else:
            handle.write(
                "REVIEW REQUIRED: Inspect the audit and solver-issue files.\n"
            )


# =============================================================================
# MAIN
# =============================================================================

def main():
    print("=" * 78)
    print("Moon HA=12 Lower-Meridian Transit Generator")
    print("Historical Julian/Gregorian calendar — Stellarium-compatible display")
    print("=" * 78)

    script_dir = os.path.dirname(os.path.abspath(__file__))

    output_file = os.path.join(
        script_dir,
        f"moon_transits_HA12_{START_YEAR}_Stellarium.csv",
    )

    interval_audit_file = os.path.join(
        script_dir,
        f"moon_transits_HA12_{START_YEAR}_interval_audit.csv",
    )

    solver_issues_file = os.path.join(
        script_dir,
        f"moon_transits_HA12_{START_YEAR}_solver_issues.csv",
    )

    duplicates_file = os.path.join(
        script_dir,
        f"moon_transits_HA12_{START_YEAR}_duplicates_removed.csv",
    )

    summary_file = os.path.join(
        script_dir,
        f"moon_transits_HA12_{START_YEAR}_run_summary.txt",
    )

    try:
        # ---------------------------------------------------------------------
        # Load offline ephemeris and time-scale data
        # ---------------------------------------------------------------------
        data_dir = os.path.join(
            os.path.dirname(skyfield_data.__file__),
            "data",
        )

        de440_path = os.path.join(data_dir, "de440.bsp")

        if not os.path.exists(de440_path):
            raise FileNotFoundError(
                f"DE440 ephemeris file was not found:\n{de440_path}"
            )

        print(f"\nSkyfield data directory:\n{data_dir}")
        print(f"\nEphemeris:\n{de440_path}")

        load = Loader(data_dir)
        ts = load.timescale()
        eph = SpiceKernel(de440_path)

        earth = eph["earth"]
        moon = eph["moon"]

        observer = earth + wgs84.latlon(
            LAT_DEG,
            LON_DEG,
            elevation_m=ELEV_M,
        )

        # ---------------------------------------------------------------------
        # Define the exact calculation interval
        # ---------------------------------------------------------------------
        start_jd = cal_to_jd(
            START_YEAR,
            START_MONTH,
            START_DAY,
            START_HOUR,
            START_MINUTE,
            START_SECOND,
        )

        end_jd = cal_to_jd(
            START_YEAR + TEST_YEARS,
            START_MONTH,
            START_DAY,
            START_HOUR,
            START_MINUTE,
            START_SECOND,
        )

        print(
            "\nProcessing interval:"
            f"\n  Start: {format_historical_datetime(start_jd)}"
            f"\n  End:   {format_historical_datetime(end_jd)}"
            f"\n  Start JD: {start_jd:.10f}"
            f"\n  End JD:   {end_jd:.10f}"
            f"\n  Scan step: {SCAN_STEP_HOURS:.3f} hours"
        )

        # ---------------------------------------------------------------------
        # Independent full-range scan
        # ---------------------------------------------------------------------
        raw_roots, solver_issues = scan_all_lower_transits(
            ts=ts,
            observer=observer,
            target=moon,
            start_jd=start_jd,
            end_jd=end_jd,
        )

        roots, duplicate_records = deduplicate_roots(raw_roots)

        # ---------------------------------------------------------------------
        # Construct outputs
        # ---------------------------------------------------------------------
        output_rows = build_output_rows(
            ts=ts,
            observer=observer,
            target=moon,
            roots=roots,
        )

        interval_audit_rows = build_interval_audit(roots)

        validation = validate_root_sequence(
            roots=roots,
            start_jd=start_jd,
            end_jd=end_jd,
        )

        # ---------------------------------------------------------------------
        # Save all files
        # ---------------------------------------------------------------------
        pd.DataFrame(output_rows).to_csv(
            output_file,
            index=False,
            encoding="utf-8-sig",
        )

        pd.DataFrame(interval_audit_rows).to_csv(
            interval_audit_file,
            index=False,
            encoding="utf-8-sig",
        )

        pd.DataFrame(solver_issues).to_csv(
            solver_issues_file,
            index=False,
            encoding="utf-8-sig",
        )

        pd.DataFrame(duplicate_records).to_csv(
            duplicates_file,
            index=False,
            encoding="utf-8-sig",
        )

        write_summary_file(
            path=summary_file,
            start_jd=start_jd,
            end_jd=end_jd,
            roots=roots,
            solver_issues=solver_issues,
            duplicate_records=duplicate_records,
            validation=validation,
        )

        # ---------------------------------------------------------------------
        # Console report
        # ---------------------------------------------------------------------
        print("\n" + "=" * 78)
        print("CALCULATION COMPLETED")
        print("=" * 78)

        print(
            f"\nRaw root detections:          {len(raw_roots):,}"
        )
        print(
            f"Unique transit records:       {len(roots):,}"
        )
        print(
            f"Duplicate detections removed: "
            f"{len(duplicate_records):,}"
        )
        print(
            f"Solver issues:                {len(solver_issues):,}"
        )

        if len(roots) > 1:
            print(
                f"Minimum interval: "
                f"{validation['minimum_interval_hours']:.9f} hours"
            )
            print(
                f"Maximum interval: "
                f"{validation['maximum_interval_hours']:.9f} hours"
            )
            print(
                f"Mean interval:    "
                f"{validation['mean_interval_hours']:.9f} hours"
            )
            print(
                f"Median interval:  "
                f"{validation['median_interval_hours']:.9f} hours"
            )

        print(
            f"\nShort intervals:              "
            f"{validation['short_interval_count']:,}"
        )
        print(
            f"Long intervals:               "
            f"{validation['long_interval_count']:,}"
        )
        print(
            f"Severe gaps:                  "
            f"{validation['severe_gap_count']:,}"
        )

        print("\nFiles saved:")
        print(f"  Main dataset:\n    {output_file}")
        print(f"  Interval audit:\n    {interval_audit_file}")
        print(f"  Solver issues:\n    {solver_issues_file}")
        print(f"  Removed duplicates:\n    {duplicates_file}")
        print(f"  Run summary:\n    {summary_file}")

        if output_rows:
            print("\nFirst five generated records:")
            print(
                pd.DataFrame(output_rows)
                .head(5)
                .to_string(index=False)
            )

            print("\nLast five generated records:")
            print(
                pd.DataFrame(output_rows)
                .tail(5)
                .to_string(index=False)
            )

        if (
            validation["strictly_increasing"]
            and validation["short_interval_count"] == 0
            and validation["long_interval_count"] == 0
            and validation["severe_gap_count"] == 0
            and len(solver_issues) == 0
        ):
            print(
                "\nVALIDATION RESULT: PASS"
                "\nNo skipped-transit pattern was detected."
            )
        else:
            print(
                "\nVALIDATION RESULT: REVIEW REQUIRED"
                "\nInspect the interval audit, solver issues, and summary files."
            )

    except Exception:
        print("\nFATAL ERROR")
        print("=" * 78)
        traceback.print_exc()
        raise


if __name__ == "__main__":
    main()
