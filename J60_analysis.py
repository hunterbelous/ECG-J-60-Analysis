# J60_analysis.py
# Arrow keys to step windows.
# Works with: J60_ecg_analysis_project/ecgPatientData.xlsx

import os
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from scipy.signal import butter, filtfilt, find_peaks

# ----------------- USER INPUT -----------------
patient = int(input("Which patient would you like to view (1 or 2)? ").strip())
if patient not in (1, 2):
    raise ValueError("Please enter 1 or 2.")

# ----------------- PATHS / COLUMNS -----------------
HERE = os.path.dirname(os.path.abspath(__file__))
EXCEL_PATH = os.path.join(HERE, "ecgPatientData.xlsx")   # <— your file here
TIME_COL = "Time"
COL_P1   = "ECG Data Patient 1"
COL_P2   = "ECG Data Patient 2"

# ----------------- WINDOW / DISPLAY -----------------
WIN_SIZE_SEC = 6.0
PCT_OVERLAP  = 0.50
N_WINDOWS    = 10
FIGSIZE      = (12, 5)
DPI          = 130

# ----------------- CLINICAL SETTINGS -----------------
BP_LOW_HZ        = 0.5
BP_HIGH_HZ       = 40.0
R_MIN_DIST_S     = 0.24
R_PROM_FRAC_PTPT = 0.15
PR_START_MS      = 200
PR_END_MS        = 120
J_MIN_MS         = 30
J_MAX_MS         = 100
J_FALLBACK_MS    = 60
J60_MS           = 60
ST_BEFORE_NEXT_R_MS = 60
R_TARGET_MV      = 1.0   # median R amplitude -> 1 mV (heuristic)

# ----------------- HELPERS -----------------
def butter_bandpass(y, fs, low, high, order=3):
    nyq = 0.5 * fs
    b, a = butter(order, [low/nyq, high/nyq], btype="band")
    return filtfilt(b, a, y)

def median_dt(tvec):
    diffs = np.diff(tvec); diffs = diffs[diffs > 0]
    return float(np.median(diffs)) if diffs.size else 1.0

def window_indices(tvec, window_sec, pct_overlap):
    if len(tvec) < 2: return []
    dt = median_dt(tvec)
    w_n  = max(2, int(round(window_sec / dt)))
    step = max(1, int(round(w_n * (1.0 - pct_overlap))))
    frames = []; start = 0
    while start + w_n <= len(tvec):
        frames.append((start, start + w_n))
        start += step
    if frames and frames[-1][1] < len(tvec):
        frames.append((len(tvec) - w_n, len(tvec)))
    return frames

def detect_r_peaks(ecg, fs):
    ptp = np.ptp(ecg) if np.ptp(ecg) > 0 else (np.max(np.abs(ecg)) + 1e-9)
    prom = max(1e-6, R_PROM_FRAC_PTPT * ptp)
    dist = max(1, int(R_MIN_DIST_S * fs))
    r_idx, _ = find_peaks(ecg, prominence=prom, distance=dist)
    return r_idx

def estimate_j_point(ecg, r_idx, fs):
    # J ~ first sustained low-slope region after R (<=160 ms window); fallback R+60 ms
    n = ecg.size
    search = int(0.16 * fs)
    hold   = max(1, int(0.012 * fs))
    deriv  = np.abs(np.diff(ecg, prepend=ecg[0]))
    j_est = []
    for r in r_idx:
        a, b = r, min(n - 1, r + search)
        if b <= a + 2:
            j_est.append(min(n - 1, r + int(J_FALLBACK_MS * fs / 1000.0))); continue
        thr = max(1e-9, 0.15 * np.max(deriv[a:b+1]))
        found = None
        for k in range(a + 1, b - hold):
            if np.all(deriv[k:k+hold] < thr):
                found = k; break
        if found is None:
            found = r + int(J_FALLBACK_MS * fs / 1000.0)
        j_est.append(min(n - 1, max(0, int(found))))
    return np.array(j_est, dtype=int)

def qualify_beats(r_idx, j_idx, fs, n):
    # keep physiologic J timing and require J+60 to be safely before next R
    r_acc, j_acc, j60_acc = [], [], []
    for i, (r, j) in enumerate(zip(r_idx, j_idx)):
        j_min = r + int(J_MIN_MS * fs / 1000.0)
        j_max = r + int(J_MAX_MS * fs / 1000.0)
        if not (j_min <= j <= min(j_max, n - 1)): continue
        j60 = j + int(J60_MS * fs / 1000.0)
        r_next = r_idx[i + 1] if i + 1 < len(r_idx) else n + 10**9
        if j60 >= r_next - int(ST_BEFORE_NEXT_R_MS * fs / 1000.0): continue
        if j60 >= n: continue
        r_acc.append(r); j_acc.append(j); j60_acc.append(j60)
    return np.array(r_acc, dtype=int), np.array(j_acc, dtype=int), np.array(j60_acc, dtype=int)

def measure_st_vs_pr(ecg_mV, t, r_idx, j_idx, j60_idx, fs):
    st_j, st_j60 = [], []
    for r, j, j60 in zip(r_idx, j_idx, j60_idx):
        a = r - int(PR_START_MS * fs / 1000.0)
        b = r - int(PR_END_MS   * fs / 1000.0)
        a = max(0, a); b = max(0, b)
        if a >= b:
            a = max(0, r - int(0.25 * fs)); b = max(0, r - int(0.12 * fs))
        if a >= b: continue
        baseline = float(np.median(ecg_mV[a:b]))
        st_j.append(float(ecg_mV[j]   - baseline))
        st_j60.append(float(ecg_mV[j60] - baseline))
    return np.array(st_j), np.array(st_j60)

def scale_to_1mV(ecg, fs):
    # global heuristic: median |R - PRbaseline| -> 1 mV
    r_raw = detect_r_peaks(ecg, fs)
    if r_raw.size < 3: return 1.0
    j_raw = estimate_j_point(ecg, r_raw, fs)
    r, j, _ = qualify_beats(r_raw, j_raw, fs, len(ecg))
    amps = []
    for ri in r:
        a = ri - int(PR_START_MS * fs / 1000.0)
        b = ri - int(PR_END_MS   * fs / 1000.0)
        a = max(0, a); b = max(0, b)
        if a >= b:
            a = max(0, ri - int(0.25 * fs)); b = max(0, ri - int(0.12 * fs))
        if a >= b: continue
        baseline = float(np.median(ecg[a:b]))
        amps.append(abs(float(ecg[ri] - baseline)))
    med = np.median(amps) if amps else 0.0
    return float(R_TARGET_MV / med) if med > 1e-9 else 1.0

# ----------------- LOAD -----------------
df = pd.read_excel(EXCEL_PATH)
if TIME_COL not in df.columns:
    for c in ("Time", "time", "t"):
        if c in df.columns:
            TIME_COL = c; break

col = COL_P1 if patient == 1 else COL_P2
for c in (TIME_COL, col):
    if c not in df.columns:
        raise ValueError(f"Missing required column: {c}")

t = pd.to_numeric(df[TIME_COL], errors="coerce").to_numpy()
y = pd.to_numeric(df[col],      errors="coerce").to_numpy()
mask = ~np.isnan(t) & ~np.isnan(y); t, y = t[mask], y[mask]
order = np.argsort(t); t, y = t[order], y[order]

dt = np.diff(t); dt = dt[dt > 0]
fs = int(round(1.0 / np.mean(dt))) if dt.size else 1

yf  = butter_bandpass(y, fs, BP_LOW_HZ, BP_HIGH_HZ)
scale = scale_to_1mV(yf, fs)
ymV = yf * scale

frames = window_indices(t, WIN_SIZE_SEC, PCT_OVERLAP)
frames = frames[:N_WINDOWS] if frames else []
if not frames:
    raise RuntimeError("Not enough data for chosen window size.")

# ----------------- PRECOMPUTE WINDOWS -----------------
win = []
for (a, b) in frames:
    tw, yw = t[a:b], ymV[a:b]
    xw = tw - tw[0]
    r_raw = detect_r_peaks(yw, fs)
    j_raw = estimate_j_point(yw, r_raw, fs)
    r, j, j60 = qualify_beats(r_raw, j_raw, fs, len(yw))
    stJ, stJ60 = measure_st_vs_pr(yw, tw, r, j, j60, fs)
    medJ   = float(np.median(stJ))   if stJ.size   else np.nan
    medJ60 = float(np.median(stJ60)) if stJ60.size else np.nan
    win.append(dict(xw=xw, yw=yw, j=j, j60=j60, beats=j.size, medJ=medJ, medJ60=medJ60))

# ----------------- VIEWER (arrow keys) -----------------
fig, ax = plt.subplots(1, 1, figsize=FIGSIZE, dpi=DPI)
idx = 0

def draw(i):
    d = win[i]
    ax.cla()
    ax.plot(d["xw"], d["yw"], lw=1.5, label="ECG (0.5–40 Hz), ~mV")
    if d["beats"] > 0:
        ax.plot(d["xw"][d["j"]],   d["yw"][d["j"]],   "kx", ms=6, label="J")
    if len(d["j60"]) > 0:
        ax.plot(d["xw"][d["j60"]], d["yw"][d["j60"]], "mD", ms=5, label="J+60 ms")
    title = (f"Patient {patient} — Window {i+1}/{len(win)} (beats:{d['beats']})   "
             f"ST@J:{d['medJ']:.3f} mV   ST@J+60:{d['medJ60']:.3f} mV"
             if d["beats"] else
             f"Patient {patient} — Window {i+1}/{len(win)} (no qualified beats)")
    ax.set_title(title)
    ax.set_xlabel("Time (s)"); ax.set_ylabel("Voltage (mV, normalized)")
    ax.legend(loc="upper right")
    ax.set_xlim(0, d["xw"][-1])
    ypad = 0.1 * (np.max(d["yw"]) - np.min(d["yw"]) + 1e-9)
    ax.set_ylim(np.min(d["yw"]) - ypad, np.max(d["yw"]) + ypad)
    fig.canvas.draw_idle()

def on_key(event):
    global idx
    if event.key in ("right", " ", "enter"):
        idx = (idx + 1) % len(win); draw(idx)
    elif event.key == "left":
        idx = (idx - 1) % len(win); draw(idx)
    elif event.key == "home":
        idx = 0; draw(idx)
    elif event.key == "end":
        idx = len(win) - 1; draw(idx)

fig.canvas.mpl_connect("key_press_event", on_key)
draw(idx)
plt.suptitle("ST vs PR baseline at J and J+60 ms — use ← / → / Space", y=0.99)
plt.tight_layout(rect=[0, 0, 1, 0.97])
print(f"[INFO] Matplotlib backend: {matplotlib.get_backend()}")
print("[INFO] Use arrow keys (←/→) or Space to move through windows.")
plt.show()