#!/usr/bin/env python3
"""
3D FeRAM CNA floorplan schematic (double-sided staircase + side periphery).

Single change vs. the original figure: the staircase->periphery connection is
no longer a horizontal wire at WBL-plane level. It is now drawn as the
physically correct path:
        via DOWN (staircase contact -> substrate metal)
     -> lateral routing along the substrate
     -> via UP (substrate metal -> periphery cell)

Every geometric / styling quantity is exposed as a parameter at the top so the
figure can be retuned without touching the drawing logic.
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Rectangle
from matplotlib.lines import Line2D

# ----------------------------------------------------------------------------
# Style
# ----------------------------------------------------------------------------
plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica", "Liberation Sans"],
    "svg.fonttype": "none",
})

BLUE     = "#3b54d4"   # WBL planes + WBL label
PILLAR   = "#1b1b1b"   # vertical channel pillars + staircase contacts
SUBST    = "#e5e7eb"   # silicon substrate fill
SUBST_EC = "#9aa0a8"   # substrate border
GREEN    = "#27a36b"   # SA / mux block
ORANGE   = "#e8731c"   # drivers / decoders block
DARK     = "#2b2b2b"   # dimension lines, routing, generic text
ROUTE    = "#1b1b1b"   # routing conductor

# ----------------------------------------------------------------------------
# Geometry parameters
# ----------------------------------------------------------------------------
N_PLANES   = 9         # number of stacked WBL planes (n layers)
ARRAY_L    = 0.70      # left edge of the bottom (widest) plane
ARRAY_R    = 7.50      # right edge of the bottom (widest) plane
STEP       = 0.22      # horizontal inset per step, per side (double-sided)
PLANE_DY   = 0.17      # vertical pitch between planes
PLANE_Y0   = 2.92      # y of the bottom plane
PLANE_LW   = 8.0       # WBL plane line width (thicker -> less whitespace)

N_PILLARS  = 7         # thick vertical channel pillars across W_array
PILLAR_LW  = 5.5

# substrate band
SUB_L, SUB_R = 0.35, 13.20
SUB_B, SUB_T = 1.55, 2.75

# periphery blocks (sit on top of the substrate, beside the array)
PERI_B, PERI_T = 2.75, 3.62
GREEN_L, GREEN_R  = 9.70, 11.35
ORANGE_L, ORANGE_R = 11.35, 12.75

# routing: down-via x, lateral y, up-via x
VIA_DN_X = ARRAY_R            # via drops at the right staircase edge
VIA_UP_X = GREEN_L - 0.40     # TAV rises just left of the periphery
LAT_Y    = 0.5 * (SUB_B + SUB_T) - 0.05   # lateral run inside the substrate

# ----------------------------------------------------------------------------
# Figure
# ----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(13.0, 4.6), dpi=200)
ax.set_xlim(0, 13.55)
ax.set_ylim(0, 5.0)
ax.set_aspect("equal")
ax.axis("off")

PLANE_YT = PLANE_Y0 + (N_PLANES - 1) * PLANE_DY      # top plane y
TOP_L = ARRAY_L + (N_PLANES - 1) * STEP              # W_array left
TOP_R = ARRAY_R - (N_PLANES - 1) * STEP              # W_array right

# --- substrate -------------------------------------------------------------
ax.add_patch(FancyBboxPatch(
    (SUB_L, SUB_B), SUB_R - SUB_L, SUB_T - SUB_B,
    boxstyle="round,pad=0,rounding_size=0.12",
    linewidth=1.3, edgecolor=SUBST_EC, facecolor=SUBST, zorder=1))
ax.text(4.05, 0.5 * (SUB_B + SUB_T), "Si substrate (one die)",
        ha="center", va="center", fontsize=12.5, color="#3a3f47", zorder=6)

# --- WBL planes (double-sided staircase) -----------------------------------
for k in range(N_PLANES):
    xl = ARRAY_L + k * STEP
    xr = ARRAY_R - k * STEP
    y  = PLANE_Y0 + k * PLANE_DY
    ax.add_line(Line2D([xl, xr], [y, y], color=BLUE, lw=PLANE_LW,
                       solid_capstyle="round", zorder=3))

# --- staircase contacts (the little "lollipop" vias on each step) ----------
def stair_contact(x, y):
    # stem descends THROUGH the stack down to the substrate (TAV)
    ax.add_line(Line2D([x, x], [SUB_T, y + 0.05], color=PILLAR, lw=1.8,
                        solid_capstyle="butt", zorder=4))
    # landing pad / bulb on top of the step
    ax.add_patch(Rectangle((x - 0.045, y + 0.03), 0.09, 0.09,
                           facecolor=PILLAR, edgecolor=PILLAR, zorder=5))

for k in range(N_PLANES - 1):                 # skip top plane (it has pillars)
    y = PLANE_Y0 + k * PLANE_DY
    stair_contact(ARRAY_L + k * STEP + 0.07, y)        # left step
    if k != 0:                                         # bottom-right shares the routing TAV
        stair_contact(ARRAY_R - k * STEP - 0.07, y)    # right step

# --- vertical channel pillars through W_array ------------------------------
for i in range(N_PILLARS):
    x = TOP_L + (i + 0.5) * (TOP_R - TOP_L) / N_PILLARS
    ax.add_line(Line2D([x, x], [PLANE_Y0 - 0.04, PLANE_YT + 0.10],
                       color=PILLAR, lw=PILLAR_LW, solid_capstyle="round",
                       zorder=4))

# --- periphery blocks ------------------------------------------------------
def block(xl, xr, color, top, sub):
    ax.add_patch(FancyBboxPatch(
        (xl, PERI_B), xr - xl, PERI_T - PERI_B,
        boxstyle="round,pad=0,rounding_size=0.07",
        linewidth=0, facecolor=color, zorder=5))
    cx = 0.5 * (xl + xr)
    ax.text(cx, PERI_B + 0.52, top, ha="center", va="center",
            fontsize=11.5, fontweight="bold", color="white", zorder=6)
    ax.text(cx, PERI_B + 0.26, sub, ha="center", va="center",
            fontsize=9.5, color="white", zorder=6)

block(GREEN_L, GREEN_R, GREEN, "SA \u00b7 mux", "accumulation")
block(ORANGE_L, ORANGE_R, ORANGE, "drivers", "decoders")

# --- routing: via DOWN -> lateral -> via UP --------------------------------
def via_pad(x, y):
    ax.add_patch(Rectangle((x - 0.055, y - 0.055), 0.11, 0.11,
                           facecolor=ROUTE, edgecolor=ROUTE, zorder=8))

# down via (staircase contact into the substrate)
ax.add_line(Line2D([VIA_DN_X, VIA_DN_X], [PLANE_Y0, LAT_Y],
                   color=ROUTE, lw=3.0, solid_capstyle="round", zorder=7))
# lateral run inside the substrate
ax.add_line(Line2D([VIA_DN_X, VIA_UP_X], [LAT_Y, LAT_Y],
                   color=ROUTE, lw=2.4, solid_capstyle="round", zorder=7))
# up via (TAV: substrate metal up toward the periphery)
ax.add_line(Line2D([VIA_UP_X, VIA_UP_X], [LAT_Y, PERI_B],
                   color=ROUTE, lw=3.0, solid_capstyle="round", zorder=7))
# small connector line from the TAV into the periphery block
ax.add_line(Line2D([VIA_UP_X, GREEN_L], [PERI_B, PERI_B],
                   color=ROUTE, lw=2.4, solid_capstyle="round", zorder=7))
for (x, y) in [(VIA_DN_X, PLANE_Y0), (VIA_DN_X, LAT_Y),
               (VIA_UP_X, LAT_Y), (VIA_UP_X, PERI_B)]:
    via_pad(x, y)

ax.text(0.5 * (VIA_DN_X + VIA_UP_X), LAT_Y + 0.20, "lateral routing",
        ha="center", va="bottom", fontsize=10.5, color="#3a3f47", zorder=8)

# ----------------------------------------------------------------------------
# Labels
# ----------------------------------------------------------------------------
ax.text(0.55, PLANE_YT + 0.28, "double-sided\nstaircase",
        ha="left", va="bottom", fontsize=10.5, color="#3a3f47")
ax.text(0.5 * (TOP_L + TOP_R), PLANE_YT + 0.40,
        "stacked WBL planes (n layers)", ha="center", va="bottom",
        fontsize=12.5, color=BLUE, fontweight="bold")
ax.text(0.5 * (GREEN_L + ORANGE_R), PERI_T + 0.22,
        "periphery sits BESIDE the array", ha="center", va="bottom",
        fontsize=11, color="#3a3f47")

# ----------------------------------------------------------------------------
# Dimension arrows
# ----------------------------------------------------------------------------
def dim(x0, x1, y, label, ly_off=-0.24):
    ax.annotate("", xy=(x0, y), xytext=(x1, y),
                arrowprops=dict(arrowstyle="<->", color=DARK, lw=1.1))
    for x in (x0, x1):
        ax.add_line(Line2D([x, x], [y - 0.07, y + 0.07], color=DARK, lw=1.0))
    ax.text(0.5 * (x0 + x1), y + ly_off, label, ha="center", va="top",
            fontsize=11, color=DARK)

DIM_Y1 = 1.18
dim(ARRAY_L, TOP_L,  DIM_Y1, "D_stair")
dim(TOP_L,  TOP_R,   DIM_Y1, "W_array")
dim(TOP_R,  ARRAY_R, DIM_Y1, "D_stair")
dim(GREEN_L, ORANGE_R, DIM_Y1, "periphery")

DIM_Y2 = 0.52
dim(SUB_L, SUB_R, DIM_Y2,
    "footprint = array + 2\u00b7staircase + periphery", ly_off=-0.24)

# ----------------------------------------------------------------------------
fig.subplots_adjust(left=0.005, right=0.995, top=0.995, bottom=0.005)
fig.savefig("/scratch/rxc5747/NeuroSim---3D-Ferroelectric/feram_floorplan.png", dpi=200,
            bbox_inches="tight", facecolor="white")
fig.savefig("/scratch/rxc5747/NeuroSim---3D-Ferroelectric/feram_floorplan.svg", bbox_inches="tight",
            facecolor="white")
print("saved")
