from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath

OUTPUT_PATH = Path(__file__).resolve().with_name("rhd_abcd_state_transition.png")

fig, ax = plt.subplots(figsize=(19, 10.5))
ax.set_xlim(0, 19)
ax.set_ylim(0, 10.5)
ax.axis("off")

navy = "#2C3E50"
text_dark = "#17202A"
muted = "#607286"
state_fill = "#EAF2F8"
input_fill = "#F7F9F9"
intervention_fill = "#FFF3CD"
intervention_edge = "#B7791F"
intervention_line = "#C27C0E"
service_fill = "#F1ECF8"
service_edge = "#6C4A8B"
death_fill = "#FDECEC"
death_edge = "#B83227"
recovery_fill = "#EAF7EE"
recovery_edge = "#2E7D32"
cycle_fill = "#EEF5F1"
cycle_edge = "#3F7C64"
white = "#FFFFFF"


def node(x, y, w, h, title, subtitle="", fill=state_fill, edge=navy,
         title_fs=13, sub_fs=9.1, z=5, linewidth=1.8):
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2), w, h,
        boxstyle="round,pad=0.04,rounding_size=0.14",
        facecolor=fill, edgecolor=edge, linewidth=linewidth, zorder=z,
    )
    patch.set_path_effects([
        pe.SimplePatchShadow(offset=(2, -2), shadow_rgbFace=(0, 0, 0), alpha=0.10),
        pe.Normal(),
    ])
    ax.add_patch(patch)
    ax.text(
        x, y + (0.17 if subtitle else 0), title,
        ha="center", va="center", fontsize=title_fs, fontweight="bold",
        color=text_dark, zorder=z + 1, linespacing=1.05,
    )
    if subtitle:
        ax.text(
            x, y - 0.27, subtitle,
            ha="center", va="center", fontsize=sub_fs, color=muted,
            zorder=z + 1, linespacing=1.13,
        )
    return patch


def arrow(start, end, label=None, color=navy, lw=2.0, rad=0.0,
          label_xy=None, linestyle="-", z=3, mutation=17, label_fs=9.2):
    patch = FancyArrowPatch(
        start, end, arrowstyle="-|>", mutation_scale=mutation,
        connectionstyle=f"arc3,rad={rad}", linewidth=lw, color=color,
        linestyle=linestyle, zorder=z,
    )
    ax.add_patch(patch)
    if label:
        lx, ly = label_xy if label_xy is not None else (
            (start[0] + end[0]) / 2, (start[1] + end[1]) / 2
        )
        ax.text(
            lx, ly, label, ha="center", va="center", fontsize=label_fs,
            color=color, zorder=z + 4, linespacing=1.12,
            bbox={"boxstyle": "round,pad=0.16", "facecolor": white,
                  "edgecolor": "none", "alpha": 0.95},
        )
    return patch


def curved_path(points, label=None, label_xy=None, color=navy, lw=1.8,
                linestyle="-", z=2, label_fs=8.9):
    path = MplPath(points, [MplPath.MOVETO, MplPath.CURVE4,
                           MplPath.CURVE4, MplPath.CURVE4])
    patch = FancyArrowPatch(
        path=path, arrowstyle="-|>", mutation_scale=16,
        linewidth=lw, color=color, linestyle=linestyle, zorder=z,
    )
    ax.add_patch(patch)
    if label and label_xy is not None:
        ax.text(
            label_xy[0], label_xy[1], label, ha="center", va="center",
            fontsize=label_fs, color=color, zorder=z + 4, linespacing=1.12,
            bbox={"boxstyle": "round,pad=0.15", "facecolor": white,
                  "edgecolor": "none", "alpha": 0.95},
        )
    return patch


def pill(x, y, text, fill, edge, width=1.85):
    patch = FancyBboxPatch(
        (x - width / 2, y - 0.23), width, 0.46,
        boxstyle="round,pad=0.03,rounding_size=0.10",
        facecolor=fill, edgecolor=edge, linewidth=1.2, zorder=8,
    )
    ax.add_patch(patch)
    ax.text(x, y, text, ha="center", va="center", fontsize=8.3,
            fontweight="bold", color=edge, zorder=9)


# Header
ax.text(9.5, 10.05, "Rheumatic Heart Disease A/B/C/D State-Transition Model",
        ha="center", va="center", fontsize=23, fontweight="bold", color=text_dark)
ax.text(
    9.5, 9.62,
    "Age–sex structured annual cohort model for Indonesia; incident RHD enters stage A and surgery is a service, not a health state",
    ha="center", va="center", fontsize=11.7, color=muted,
)

# Top callouts
node(
    2.55, 8.55, 4.25, 1.25, "Data-fed epidemiology",
    "GBD 2023 age–sex incidence,\nprevalence, RHD mortality,\nand other-cause mortality",
    fill=input_fill, edge=navy, title_fs=12.2, sub_fs=8.8,
)
node(
    9.45, 8.55, 5.10, 1.38, "Secondary-prevention care cascade",
    "Screening → diagnosis → optimal treatment (SAP)\n"
    "Treatment is capped by diagnosis and screening\n"
    "SAP reduces stage-specific RHD mortality in A, B, C, and D",
    fill=intervention_fill, edge=intervention_edge, title_fs=12.3, sub_fs=8.6,
)
node(
    16.15, 8.55, 3.75, 1.30, "Surgery service",
    "Applied to fractions of stages C and D\n"
    "Reduces C→D progression and D→RHD death\n"
    "No surgery or post-surgery state",
    fill=service_fill, edge=service_edge, title_fs=12.0, sub_fs=8.5,
)

# Health states
node(1.25, 5.75, 2.20, 1.25, "No RHD", "Population without\nprevalent disease",
     fill=recovery_fill, edge=recovery_edge, title_fs=13.0)
node(4.20, 5.75, 2.35, 1.36, "Stage A", "Early / mild RHD\nIncident cases enter here",
     fill=state_fill, edge=navy, title_fs=13.4)
node(7.30, 5.75, 2.35, 1.36, "Stage B", "Established RHD\nBidirectional with A and C",
     fill=state_fill, edge=navy, title_fs=13.4)
node(10.40, 5.75, 2.35, 1.36, "Stage C", "Advanced RHD\nMay require surgery",
     fill=state_fill, edge=navy, title_fs=13.4)
node(13.50, 5.75, 2.35, 1.36, "Stage D", "Severe / end-stage RHD\nHighest RHD mortality",
     fill=state_fill, edge=navy, title_fs=13.4)
node(17.10, 3.10, 2.85, 1.25, "Death", "Stage-specific RHD death\nor other-cause death",
     fill=death_fill, edge=death_edge, title_fs=13.1, sub_fs=8.9)

# Main progression
arrow((2.35, 5.75), (3.00, 5.75), "Incident RHD", color=navy, lw=2.2,
      label_xy=(2.67, 6.10), label_fs=8.8)
arrow((5.38, 5.75), (6.12, 5.75), "A → B", color=navy, lw=2.15,
      label_xy=(5.75, 6.08), label_fs=8.7)
arrow((8.48, 5.75), (9.22, 5.75), "B → C", color=navy, lw=2.15,
      label_xy=(8.85, 6.08), label_fs=8.7)
arrow((11.58, 5.75), (12.32, 5.75), "C → D", color=service_edge, lw=2.6,
      label_xy=(11.95, 6.08), label_fs=8.7)

# Regression / reverse transitions
curved_path([(3.85, 5.08), (3.60, 4.18), (2.15, 4.15), (1.82, 5.10)],
            "A → No RHD", (2.77, 4.35), recovery_edge, 1.8, label_fs=8.7)
curved_path([(6.95, 5.08), (6.70, 4.55), (5.75, 4.55), (5.48, 5.08)],
            "B → A", (6.22, 4.60), cycle_edge, 1.7, label_fs=8.6)
curved_path([(10.05, 5.08), (9.80, 4.30), (8.85, 4.30), (8.58, 5.08)],
            "C → B", (9.32, 4.35), cycle_edge, 1.7, label_fs=8.6)
curved_path([(13.15, 5.08), (12.90, 4.05), (11.95, 4.05), (11.68, 5.08)],
            "D → C\n(optional; default 0)", (12.42, 4.18), cycle_edge, 1.6,
            linestyle="--", label_fs=8.3)

# Prevalence seed
curved_path([(2.50, 7.90), (4.00, 7.35), (6.25, 7.08), (7.12, 6.47)],
            "Prevalent RHD seed split across A/B/C/D", (5.55, 7.25),
            navy, 1.55, linestyle="--", label_fs=8.6)

# Intervention overlays
for x in (4.20, 7.30, 10.40, 13.50):
    arrow((9.45, 7.82), (x, 6.45), color=intervention_line, lw=1.25,
          linestyle="--", z=2, mutation=12)
    pill(x, 6.75, "SAP lowers RHD death", intervention_fill, intervention_edge)

arrow((16.15, 7.88), (11.95, 6.20), color=service_edge, lw=1.55,
      linestyle="--", z=3, mutation=13)
arrow((16.15, 7.88), (14.55, 5.10), color=service_edge, lw=1.55,
      linestyle="--", z=3, mutation=13)

# Mortality pathways from A-D
curved_path([(4.55, 5.06), (6.30, 3.55), (12.40, 3.15), (15.72, 3.12)],
            "Stage A RHD death + other-cause death", (8.35, 3.62),
            death_edge, 1.55, label_fs=8.2)
curved_path([(7.65, 5.06), (9.10, 3.95), (13.45, 3.55), (15.72, 3.20)],
            "Stage B mortality", (10.70, 3.95), death_edge, 1.55, label_fs=8.3)
curved_path([(10.75, 5.06), (11.70, 4.40), (14.65, 3.80), (15.78, 3.30)],
            "Stage C mortality", (13.45, 4.30), death_edge, 1.65, label_fs=8.3)
curved_path([(13.85, 5.06), (14.60, 4.45), (15.35, 3.85), (15.92, 3.45)],
            "Stage D mortality\n(surgery lowers RHD component)", (15.00, 4.45),
            death_edge, 1.85, label_fs=8.3)

# Annual ageing
node(
    3.35, 2.10, 4.45, 1.20, "Annual ageing and carry-forward",
    "Survivors remain in their resulting state,\nthen age a → a + 1; age 100 is open",
    fill=cycle_fill, edge=cycle_edge, title_fs=11.8, sub_fs=8.8, linewidth=1.6,
)
for x in (1.25, 4.20, 7.30, 10.40, 13.50):
    curved_path([(x, 5.08), (x - 0.10, 3.65), (5.10, 2.90), (4.85, 2.40)],
                color=cycle_edge, lw=1.15, linestyle="--", z=1)
arrow((1.50, 2.62), (1.20, 5.06), "Next annual cycle", color=cycle_edge,
      lw=1.55, rad=-0.20, label_xy=(0.85, 3.72), linestyle="--",
      mutation=13, label_fs=8.5)

# Legend
legend_y = 0.78
ax.text(0.45, legend_y, "Legend", fontsize=10.2, fontweight="bold",
        color=text_dark, va="center")


def legend_swatch(x, fill, edge, label):
    swatch = FancyBboxPatch(
        (x, legend_y - 0.18), 0.42, 0.36,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor=fill, edgecolor=edge, linewidth=1.25,
    )
    ax.add_patch(swatch)
    ax.text(x + 0.55, legend_y, label, fontsize=9.1, color=text_dark, va="center")


legend_swatch(1.15, recovery_fill, recovery_edge, "No RHD / regression")
legend_swatch(4.05, state_fill, navy, "Living WHF RHD stage")
legend_swatch(7.15, intervention_fill, intervention_edge, "Secondary prevention")
legend_swatch(11.30, service_fill, service_edge, "Surgery service")
legend_swatch(14.15, death_fill, death_edge, "Absorbing death state")

ax.text(0.45, 0.25, "Health states: No RHD, A, B, C, D, Death",
        ha="left", va="center", fontsize=9.0, color=muted, fontweight="bold")
ax.text(18.55, 0.25,
        "Reference: cascade held at baseline | Scale-up: screening, diagnosis, and treatment ramp to targets",
        ha="right", va="center", fontsize=8.9, color=muted)

plt.tight_layout(pad=0.8)
fig.savefig(OUTPUT_PATH, dpi=300, bbox_inches="tight", facecolor="white")
plt.close(fig)
print(f"Saved PNG: {OUTPUT_PATH}")
