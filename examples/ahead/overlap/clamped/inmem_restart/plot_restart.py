#!/usr/bin/env python3
"""Animate the in-memory restart march.  Run after restart_example1.jl.

sd1 switches FOM -> ROM at step 40 and back at step 81; sd2 stays FOM. The
shading follows the live mode, so the switch is visible as the left subdomain
turning blue and then red again.

Outputs  restart_switch.gif        one frame per control step
         restart_switch_strip.png  four frames spanning the ROM window
"""
import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from PIL import Image

H = os.path.dirname(os.path.abspath(__file__))
p = lambda f: os.path.join(H, f)

SD = [(-0.50, -0.20), (-0.30, 0.50)]       # sd1, sd2
TFINAL, FRAME_MS = 4.0e-4, 55
FOM_C, ROM_C, INK = "#c0392b", "#2c6fbb", "#1a1a1a"

z = np.loadtxt(p("out_zgrid.csv"), delimiter=",")
F = np.loadtxt(p("out_field.csv"), delimiter=",")
with open(p("out_modes.csv")) as fh:
    modes = [r[1] for r in list(csv.reader(fh))[1:]]
n = min(F.shape[0], len(modes))
amp = np.abs(F).max()
print(f"{n} steps, {len(z)} z-points")


def draw(ax, k):
    ax.clear()
    m1 = modes[k]                            # sd1 is the one that switches
    ax.axvspan(*SD[0], color=ROM_C if m1 == "ROM" else FOM_C, alpha=0.22, lw=0)
    ax.axvspan(*SD[1], color=FOM_C, alpha=0.22, lw=0)
    ax.plot(z, F[k], lw=1.8, color=INK, zorder=3)
    ax.axhline(0, lw=0.5, color="0.75")
    ax.set(xlim=(z[0], z[-1]), ylim=(-1.15 * amp, 1.15 * amp),
           xlabel="z", ylabel=r"$u_z$")
    ax.set_title(f"step {k+1}/{n}      sd1 = {m1}", fontsize=11)
    ax.text(0.985, 0.965, "red = FOM   blue = ROM", transform=ax.transAxes,
            ha="right", va="top", fontsize=7.5, color="0.35")
    ax.grid(alpha=0.25)


fig, ax = plt.subplots(figsize=(7.6, 3.9))
frames = []
for k in range(n):
    draw(ax, k)
    fig.tight_layout(); fig.canvas.draw()
    frames.append(Image.frombytes("RGBA", fig.canvas.get_width_height(),
                                  fig.canvas.buffer_rgba().tobytes()).convert("P"))
frames[0].save(p("restart_switch.gif"), save_all=True, append_images=frames[1:],
               duration=FRAME_MS, loop=0, optimize=True)
plt.close(fig)
print("wrote restart_switch.gif")

# static strip: just before the switch, inside the ROM window, and after
picks = [30, 55, 75, 100]
fig, axes = plt.subplots(2, 2, figsize=(11.0, 5.6))
for a, k in zip(axes.ravel(), picks):
    draw(a, min(k, n - 1))
fig.tight_layout(); fig.savefig(p("restart_switch_strip.png"), dpi=140)
print("wrote restart_switch_strip.png   frames", picks)
