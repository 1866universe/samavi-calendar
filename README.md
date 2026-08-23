# Samavi Calendar

## An Astronomical Research Model Based on Meridian Transit and Lunar Angular Cycles

**Official Live Application:** [https://1866universe.github.io/samavi-calendar/](https://1866universe.github.io/samavi-calendar/)

Drawing upon a data repository containing **232,849 precise astronomical points in time**, this system enables the calculation, tracking, and display of calendar data over a range equivalent to **693 astronomical years**, spanning from `1950-04-18` to `2610-02-25`.

The **Samavi Calendar** is an astronomical research model designed to analyze time on the basis of actual celestial events. In this calendar, time is not defined merely as a conventional division such as midnight, fixed days, or administrative months. Rather, its principal basis is the Moon's actual position relative to the observer's meridian and the Moon's angular cycles.

---

## Core Astronomical Principles

### 1. The Reference Starting Point
In the Samavi Calendar, the starting point of the time cycle is defined as the moment when the Moon is in a meridian-transit position and its hour angle reaches the reference value of `00h 00m 00s`. This moment constitutes the first reference point in the structure of the Samavi Calendar.

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
|---|---|
| **Samavi Epoch** | `0001-01-01 12:44:35` (Tuesday, Iran Local Time) |
| **Gregorian Equivalent** | `1950-04-18 09:14:35 UTC` |
| **Solar Hijri Equivalent** | `1329-01-29 12:44:35 IRST` |
| **Hijri Qamari Equivalent** | `1369-06-29 09:14:35 UTC` |

This moment is regarded as the reference point for all temporal calculations in the Samavi Calendar, and all subsequent days, months, and years are determined based on the astronomical time interval from this point.

---

## Download Guide & Data Files

The official research datasets and Excel models are available in the Assets section of the GitHub releases:

**Direct Link:** [Download Samavi Calendar Datasets & Releases](https://github.com/1866universe/samavi-calendar/releases)

[راهنمای فارسی](README_FA.md)

---
*1866Universe 2026 ©*