# Samavi Calendar

## An Astronomical Research Model Based on Meridian Transit and Lunar Angular Cycles

**Official Live Application:** [https://1866universe.github.io/samavi-calendar/](https://1866universe.github.io/samavi-calendar/)

Drawing upon a data repository containing **232,849 precise astronomical points in time**, this system enables the calculation, tracking, and display of calendar data over a range equivalent to **693 astronomical years**, spanning from `1552-03-24` to `2212-03-24`.

The **Samavi Calendar** is an astronomical research model designed to analyze time on the basis of actual celestial events. In this calendar, time is not defined merely as a conventional division such as midnight, fixed days, or administrative months. Rather, its principal basis is the Moon's actual position relative to the observer's meridian and the Moon's angular cycles.

---

## Core Astronomical Principles

### 1. The Reference Starting Point

In the Samavi Calendar, the starting point of the time cycle is defined as the moment when the Moon is in a meridian-transit position and its hour angle reaches the reference value of `12h 00m 00s`. This moment constitutes the first reference point in the structure of the Samavi Calendar.

### 2. Month Structure & Orbital Phase Switch

A Samavi month is defined on the basis of a 28-point interval:

- **Points 1 through 27:** Represent the main path of the Moon's angular cycle.
- **Point 28:** Functions as a transition point, or **orbital phase switch**, to the following phase. At Point 28, the Moon becomes invisible to the observer; at this time, an orbital phase change occurs.

The subsequent 27-day cycle then begins, until the next 28th point is reached and the orbital phase switch is repeated continuously for every astronomical month.

### 3. Annual Cycle (336 Days)

A Samavi year consists of 12 months of 28 days each:
$$\text{Length of Year} = 12 \times 28 = 336 \text{ days}$$

This 336-day structure is not merely an administrative or conventional division. Rather, it is an intelligent natural system of 12 complete lunar cycles that occurs within the Moon’s annual orbital cycle. In other words, each astronomical year includes 12 cycles of 28-point cycles which, in total, complete a quasi-sinusoidal wave pattern over this period.

---

## Principle of Periodicity and the 18-Year Lunar Orbital Recurrence Cycle Based on Samavi Calendar Documentation
Based on the extracted empirical dataset and the behavioral trends recorded in the astronomical day duration charts of the Samavi Calendar, the annual orbital movement and behavior of the Moon exhibit a fully reproducible and closed cyclical pattern occurring in 18-year intervals (equivalent to 216 Samavi months or 6,048 Samavi days).

### Independence of the Pattern from the Starting Point (Arbitrary Baseline):
This cycle is not tied to any specific year. Whichever point in time is chosen as the baseline, exactly upon the completion of a full 18-year period, the orbital behavior of the 19th year identically mirrors that of the baseline year.

### For instance:
If year 1405 is selected as an arbitrary starting point:
The cycle concludes over an 18-year span (from 1405 to the end of 1422).
In year 1423 (the 19th year), the behavioral pattern of year 1405 identically repeats for another 18-year period.
The very same rule holds true for any arbitrary year (such as 1300, 1356, 2070, etc.).
This cycle is recorded as a fixed and proven periodic rule in the Samavi Calendar documentation, demonstrating the mathematical stability and absolute reproducibility of the database records without the slightest deviation.

📜 The recorded documentation is available at the bottom of the Astronomical_Day_Duration_Trend sheet in the Samavi Calendar Excel file. Furthermore, exploration and research across daily, monthly, and annual periods are facilitated through the integrated charting system within the same sheet.

## Overview of Month Structures in Conventional Calendars

To better understand the structure of the astronomical calendar and its differences from commonly used calendars, the conventional structure and month lengths of three calendar systems are reviewed below:

### 1. Gregorian Calendar

The Gregorian calendar is a solar calendar with 12 months and approximately 365 days per year.

- **31 Days:** January (1), March (3), May (5), July (7), August (8), October (10), December (12)
- **30 Days:** April (4), June (6), September (9), November (11)
- **February (2):** 28 days in common years and 29 days in leap years.
- **Year Length:** 365 days in a common year and 366 days in a leap year.

### 2. Solar Hijri (Jalali) Calendar

The Solar Hijri calendar is a solar calendar with 12 months, whose beginning is aligned with the vernal equinox (Nowruz).

- **31 Days:** The first half of the year, including Farvardin (1) through Shahrivar (6)
- **30 Days:** Mehr (7) through Bahman (11)
- **Esfand (12):** 29 days in common years and 30 days in leap years.
- **Year Length:** 365 days in a common year and 366 days in a leap year.

### 3. Hijri Qamari Calendar

The Hijri Qamari calendar is an observation-based lunar calendar based on the Moon’s orbital motion, consisting of 12 months and approximately 354 days per year. In computational systems (conventional/tabular), months alternate between 30 and 29 days.

- **30 Days:** Muharram (1), Rabi al-Awwal (3), Jumada al-Awwal (5), Rajab (7), Ramadan (9), Dhu al-Qidah (11)
- **29 Days:** Safar (2), Rabi al-Thani (4), Jumada al-Thani (6), Shaban (8), Shawwal (10)
- **Dhu al-Hijjah (12):** 29 days in common years and 30 days in Hijri Qamari leap years.
- **Year Length:** 354 days in a common year and 355 days in a Hijri Qamari leap year.

---

## Computational Differences in Samavi Calendar

One of the important differences between the Samavi Calendar and conventional calendars is that, in this system, calculations are performed on the basis of precise astronomical points in time. Therefore, the length of time intervals is not determined solely by simple day counting; rather, the actual recorded time of meridian transits, as astronomical reference points, is incorporated into the calculations.

For this reason, two intervals with the same number of days may not be exactly equal in their total number of minutes. This difference results from the dynamic nature of celestial motion and minor variations in the timing of astronomical events.

> **Scope Note:** The Samavi Calendar represents lunar astronomical events and differs from official, civil, or religious calendars. The purpose of this calendar is to provide a research model for examining time on the basis of the actual behavior of celestial bodies.

---

## Initial Reference Epoch

The Samavi Calendar begins from a specific astronomical reference moment. This moment is the zero point of the day-counting system, and all dates are calculated relative to it:

| System | Date / Time Representation |
| --- | --- |
| **Samavi Epoch** | `0001-01-01 12:44:35` (Tuesday, Iran Local Time) |
| **Gregorian Equivalent** | `1950-04-18 09:14:35 UTC` |
| **Solar Hijri Equivalent** | `1329-01-29 12:44:35 IRST` |
| **Hijri Qamari Equivalent** | `1369-06-29 09:14:35 UTC` |

This moment is regarded as the reference point for all temporal calculations in the Samavi Calendar, and all subsequent days, months, and years are determined based on the astronomical time interval from this point.

---

## Computational Engine and Data Transparency (Astronomical Scripts)

All temporal reference points and meridian transits recorded in the Samavi Calendar are derived directly from NASA’s orbital data (JPL Ephemeris) via the Skyfield astronomical computation engine. To maintain scientific transparency and provide researchers with the ability to perform independent verification, the Python scripts used to generate this database are provided in the [`scripts/`](scripts/) directory:

* Upper Transit Script (HA = 00h): Precise tracking of the Moon’s meridian transit moments from the observer’s local meridian (moon_transits_skyfield_ha00.py).
* Lower Transit Script (HA = 12h): Precise tracking of the Moon’s transit moments at the anti-meridian (180 degrees from the observer’s local meridian) (moon_transits_skyfield_ha12.py).

## Implementation and Technical Prerequisites:

* Base Ephemeris: Calculations are performed using the Skyfield software package and the DE440 planetary positioning file (de440.bsp) with the highest floating-point precision.
* Automaic Download and Deployment: The script structure is designed to automatically download and store the de440.bsp file from the authoritative astronomical servers upon the first execution if it is missing from the local directory; no manual intervention is required.
* Verification Reports: Detailed documentation of stability tests and successful error-free passage through critical historical junctures (including the transition from the Julian to the Gregorian calendar in 1582 AD) is available to researchers in the [`scripts/validation/`](scripts/validation/) directory.

---

## Download Guide & Data Files

The official research datasets and Excel models are available in the Assets section of the GitHub releases:

**Direct Link:** [Download Samavi Calendar Datasets & Releases](https://github.com/1866universe/samavi-calendar/releases)

[راهنمای فارسی](README_FA.md)

---

1866Universe 2026 ©
