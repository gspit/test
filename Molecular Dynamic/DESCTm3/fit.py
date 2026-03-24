import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import linregress

# ==========================================================
# FILE
# ==========================================================

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_FILE = os.path.join(SCRIPT_DIR, "H_T.dat")

# ==========================================================
# SETTINGS
# ==========================================================

T_MIN = 300
T_MAX = 900
MIN_POINTS = 6

# ==========================================================
# READ DATA
# ==========================================================

data = []
with open(DATA_FILE) as f:
    for line in f:
        if line.startswith("#") or not line.strip():
            continue
        T, _, V = map(float, line.split()[:3])
        if T_MIN <= T <= T_MAX:
            data.append([T, V])

df = pd.DataFrame(data, columns=["T", "V"]).sort_values("T")
T = df["T"].values
V = df["V"].values

# ==========================================================
# CLASSIC GEOMETRIC BILINEAR FIT
# ==========================================================

best_Tm = None
best_err = np.inf
best_fit = None

for i in range(MIN_POINTS, len(T) - MIN_POINTS):

    Tl, Vl = T[:i], V[:i]
    Th, Vh = T[i:], V[i:]

    sl, il, *_ = linregress(Tl, Vl)
    sh, ih, *_ = linregress(Th, Vh)

    if sh <= sl:
        continue

    # intersection temperature
    Tm = (ih - il) / (sl - sh)

    if Tm < Tl[-1] or Tm > Th[0]:
        continue

    err = np.sum((Vl - (sl * Tl + il))**2) + \
          np.sum((Vh - (sh * Th + ih))**2)

    if err < best_err:
        best_err = err
        best_Tm = Tm
        best_fit = (sl, il, sh, ih)

sl, il, sh, ih = best_fit

# ==========================================================
# OUTPUT
# ==========================================================

print("\n=========== CLASSIC MELTING POINT RESULT ===========")
print(f"Tm (V–T intersection) = {best_Tm:.1f} K")
print("====================================================\n")

# ==========================================================
# PLOT
# ==========================================================

plt.figure(figsize=(7,5))
plt.scatter(T, V, s=40, label="MD data")

Tfit = np.linspace(T.min(), T.max(), 400)
plt.plot(Tfit, sl*Tfit + il, "--", lw=2, label="Solid fit")
plt.plot(Tfit, sh*Tfit + ih, "--", lw=2, label="Liquid fit")

plt.scatter(best_Tm, sl*best_Tm + il, color="red", zorder=5)
plt.axvline(best_Tm, color="red", linestyle=":", lw=2,
            label=f"Tm = {best_Tm:.0f} K")

plt.xlabel("Temperature (K)")
plt.ylabel("Volume (Å³)")
plt.legend()
plt.tight_layout()
plt.show()
