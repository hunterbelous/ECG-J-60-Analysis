# ECG ST-Segment Analysis & Visualization

**Project overview (≈150 words).**  
This project analyzes and visualizes electrocardiogram (ECG) data from an Excel file to explore **ventricular repolarization** behavior. The pipeline loads time-series signals for two patients, band-pass filters them (0.5–40 Hz), detects **R-peaks**, estimates the **J-point** (end of QRS), and computes **ST-segment deviation relative to the PR baseline** at both **J** and **J+60 ms**—a commonly used measurement time in stress-test interpretation. Beats are **gated** to keep only physiologic detections (plausible J timing and J+60 safely before the next R) so averages are robust. Signals are then **normalized** so the median R-peak amplitude ≈ **1 mV** (heuristic), which makes magnitudes easier to interpret when true device gain is unknown. The first 10 windows (6 s, 50% overlap) are rendered **side-by-side** (Patient 1 | Patient 2) into an animated GIF for quick review. The output is intended for **exploratory signal analysis and education**—it is **not** a diagnostic tool.

---

## Visualization

![ST Windows (J and J+60 ms)](figures/st_windows_st_clinical_gated_norm.gif)  
If the GIF does not render, [open it directly](figures/st_windows_st_clinical_gated_norm.gif).

---

## Representative Window Results (mV, normalized)

- **Patient 1** — accepted beats: **12**  
  - **ST@J:** **−0.058 mV**  
  - **ST@J+60:** **−0.060 mV**  
  - *Interpretation:* mild, fairly flat ST depression in this lead.

- **Patient 2** — accepted beats: **5**  
  - **ST@J:** **−0.268 mV**  
  - **ST@J+60:** **−0.185 mV**  
  - *Interpretation:* more pronounced ST depression at J with partial recovery by 60 ms.

---

## Unofficial, Non-Diagnostic Summary

- In this representative window, **Patient 2 shows substantially greater ST depression** than Patient 1 (≈3–4× larger at J).  
- **J→J+60 trend:**  
  - **P1:** ST remains about the same (near-horizontal ST segment around −0.06 mV).  
  - **P2:** ST moves **toward baseline** by 60 ms, consistent with upsloping/horizontal depression in this lead.  
- Taken together, P2’s tracing shows a **stronger repolarization abnormality** in this analysis compared with P1 **in the displayed lead and window**.

---

## Caveats

- Values are in **normalized mV** (median R ≈ 1 mV); without true calibration, thresholds (e.g., ≥0.1 mV) are **approximate**.  
- Conclusions depend on **lead placement, filtering, and baseline selection**; only a subset of windows/leads is shown here.  
- Clinical decisions require **contiguous leads**, **consistent changes over time**, and **clinical context** (symptoms, workload).  
- This repository is for **exploratory analysis/education** and is **not a medical device** or diagnostic tool.