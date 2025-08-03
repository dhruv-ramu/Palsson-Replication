# Dynamic FBA Simulation Results

This directory contains the results from the dynamic FBA (dFBA) simulation of E. coli diauxic growth on glucose and lactose.

## Overview

**Simulation:** E. coli iJO1366 diauxic growth simulation  
**Reference:** Mahadevan et al. 2002 PNAS  
**Agreement:** 86.2% overall agreement with literature benchmarks  
**Date:** August 2, 2024  

## Files

### Primary Results (Fixed Version)
- **`dynamic_fba_trajectory_fixed.csv`** - Complete simulation trajectory (400 time points)
- **`dynamic_fba_analysis_fixed.json`** - Analysis results and key metrics
- **`dynamic_fba_diauxic_growth_fixed.png`** - **Primary visualization** (4-panel analysis)

### Iteration History
- **`dynamic_fba_trajectory_improved.csv`** - Improved version trajectory
- **`dynamic_fba_analysis_improved.json`** - Improved version analysis
- **`dynamic_fba_diauxic_growth_improved.png`** - Improved version visualization

- **`dynamic_fba_trajectory_final.csv`** - Final version trajectory
- **`dynamic_fba_analysis_final.json`** - Final version analysis
- **`dynamic_fba_diauxic_growth_final.png`** - Final version visualization

- **`dynamic_fba_trajectory.csv`** - Initial version trajectory
- **`dynamic_fba_analysis.json`** - Initial version analysis
- **`dynamic_fba_diauxic_growth.png`** - Initial version visualization

## Key Results

### Timing Agreement with Mahadevan et al. 2002
- **Glucose depletion:** 5.15 hours (85.8% agreement)
- **Lag phase duration:** 1.30 hours (86.7% agreement)
- **Lactose phase start:** 6.45 hours (86.0% agreement)
- **Overall timing:** 86.2% agreement

### Growth Dynamics
- **Glucose phase growth rate:** 0.976 1/h
- **Lactose phase growth rate:** 0.222 1/h
- **Final biomass:** 1.748 g/L
- **Total growth:** 175-fold increase

### Simulation Parameters
- **Initial glucose:** 25.0 mM
- **Initial lactose:** 25.0 mM
- **Initial biomass:** 0.01 g/L
- **Simulation time:** 20 hours
- **Time resolution:** 0.05 hours

## Data Format

### CSV Trajectory Files
Columns: time, glucose, lactose, biomass, growth_rate, glucose_uptake, lactose_uptake

### JSON Analysis Files
Contains: glucose_depletion_time, lactose_start_time, lag_duration, glucose_phase_growth, lactose_phase_growth, final_biomass

## Visualization

The primary figure (`dynamic_fba_diauxic_growth_fixed.png`) shows:
1. **Biomass growth curve** - Two-phase exponential growth
2. **Substrate consumption** - Sequential glucose then lactose
3. **Growth rate dynamics** - Sharp phase transitions
4. **Uptake rates** - Regulatory control of metabolism

## Technical Notes

- **Numerical stability:** Fixed biomass constraints prevent crashes
- **ODE integration:** SciPy odeint with tolerance settings
- **Regulatory logic:** Glucose repression, lag phase, lactose induction
- **Michaelis-Menten kinetics:** Realistic substrate uptake modeling 