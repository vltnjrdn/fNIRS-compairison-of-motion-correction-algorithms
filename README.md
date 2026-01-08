# fNIRS: Comparison of Motion Correction Algorithms

This repository contains scripts and datasets for analyzing fNIRS data with different motion correction algorithms.

---

## Scripts

There are **4 main scripts**:

1. **Whole dataset analysis**
   - Processes the entire dataset.
   - Available **with and without SSR**.

2. **Time & channel specific analysis**
   - Analyzes a user-defined time window and channel.
   - Available **with and without SSR**.

> Naming convention:
> - Scripts with `_time_channel` operate on a specific time window and channel.
> - Scripts with `_SSR` apply SSR correction.

---

## Datasets

Two datasets are included:

- `RLOESS.mat`
- `RLOESS_SSR.mat`

> These were saved separately because the algorithm processing is time-consuming.

---

## Functions

- `compute_PSD_oxyHb.m`:  
  Computes Power Spectral Density (PSD) metrics for oxygenated hemoglobin (HbO) across different frequency bands.

---
