from pathlib import Path

import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch
from matplotlib.path import Path as MplPath


# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
OUTPUT_PATH = Path(__file__).resolve().with_name(
    "rhd_secondary_prevention_state_transition.png"
)


# -----------------------------------------------------------------------------
# Figure setup
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(18, 10))
ax.set_xlim(0, 18)
ax.set_ylim(0, 10)
ax.axis("off")


# -----------------------------------------------------------------------------
# Palette
# -----------------------------------------------------------------------------
navy = "#2C3E50"
text_dark = "#17202A"
muted = "#607286"
state_fill = "#EAF2F8"
input_fill = "#F7F9F9"

intervention_fill = "#FFF3CD"
intervention_edge = "#B7791F"
intervention_line = "#C27C0E"

purple_fill = "#F1ECF8"
purple_edge = "#6C4A8B"

red_fill = "#FDECEC"
red_edge = "#B83227"

cycle_fill = "#EEF5F1"
cycle_edge = "#3F7C64"

white = "#FFFFFF"


# -----------------------------------------------------------------------------
# Drawing helpers
# -----------------------------------------------------------------------------
def node(
    x,
    y,
    w,
    h,
    title,
    subtitle="",
    fill=state_fill,
    edge=navy,
    title_fs=13,
    sub_fs=9.5,
    z=5,
    linewidth=1.8,
):
    """Draw a rounded state, input, or intervention node."""
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.035,rounding_size=0.14",
        facecolor=fill,
        edgecolor=edge,
        linewidth=linewidth,
        zorder=z,
    )
    patch.set_path_effects(
        [
            pe.SimplePatchShadow(
                offset=(2, -2),
                shadow_rgbFace=(0, 0, 0),
                alpha=0.10,
            ),
            pe.Normal(),
        ]
    )
    ax.add_patch(patch)

    title_y = y + (0.17 if subtitle else 0)
    ax.text(
        x,
        title_y,
        title,
        ha="center",
        va="center",
        fontsize=title_fs,
        fontweight="bold",
        color=text_dark,
        zorder=z + 1,
        linespacing=1.05,
    )

    if subtitle:
        ax.text(
            x,
            y - 0.27,
            subtitle,
            ha="center",
            va="center",
            fontsize=sub_fs,
            color=muted,
            zorder=z + 1,
            linespacing=1.14,
        )
    return patch


def arrow(
    start,
    end,
    label=None,
    color=navy,
    lw=2.0,
    rad=0.0,
    label_xy=None,
    linestyle="-",
    z=3,
    mutation=17,
    label_fs=9.5,
):
    """Draw a straight or curved arrow."""
    patch = FancyArrowPatch(
        start,
        end,
        arrowstyle="-|>",
        mutation_scale=mutation,
        connectionstyle=f"arc3,rad={rad}",
        linewidth=lw,
        color=color,
        linestyle=linestyle,
        zorder=z,
    )
    ax.add_patch(patch)

    if label:
        lx, ly = (
            label_xy
            if label_xy is not None
            else ((start[0] + end[0]) / 2, (start[1] + end[1]) / 2)
        )
        ax.text(
            lx,
            ly,
            label,
            ha="center",
            va="center",
            fontsize=label_fs,
            color=color,
            zorder=z + 3,
            linespacing=1.13,
            bbox={
                "boxstyle": "round,pad=0.16",
                "facecolor": white,
                "edgecolor": "none",
                "alpha": 0.95,
            },
        )
    return patch


def curved_path(
    points,
    label=None,
    label_xy=None,
    color=red_edge,
    lw=1.9,
    linestyle="-",
    z=2,
    label_fs=9.2,
):
    """Draw a cubic Bézier arrow from four control points."""
    path = MplPath(
        points,
        [
            MplPath.MOVETO,
            MplPath.CURVE4,
            MplPath.CURVE4,
            MplPath.CURVE4,
        ],
    )
    patch = FancyArrowPatch(
        path=path,
        arrowstyle="-|>",
        mutation_scale=16,
        linewidth=lw,
        color=color,
        linestyle=linestyle,
        zorder=z,
    )
    ax.add_patch(patch)

    if label and label_xy is not None:
        ax.text(
            label_xy[0],
            label_xy[1],
            label,
            ha="center",
            va="center",
            fontsize=label_fs,
            color=color,
            zorder=z + 4,
            linespacing=1.12,
            bbox={
                "boxstyle": "round,pad=0.15",
                "facecolor": white,
                "edgecolor": "none",
                "alpha": 0.95,
            },
        )
    return patch


# -----------------------------------------------------------------------------
# Header
# -----------------------------------------------------------------------------
ax.text(
    9,
    9.63,
    "Rheumatic Heart Disease Secondary-Prevention Model",
    ha="center",
    va="center",
    fontsize=23,
    fontweight="bold",
    color=text_dark,
)
ax.text(
    9,
    9.22,
    (
        "Age–sex structured annual cohort model for Indonesia (2026–2100); "
        "ARF is not represented as a health state"
    ),
    ha="center",
    va="center",
    fontsize=11.8,
    color=muted,
)


# -----------------------------------------------------------------------------
# Top callouts
# -----------------------------------------------------------------------------
node(
    2.55,
    8.10,
    4.15,
    1.25,
    "Data-fed epidemiology",
    "GBD 2023 age–sex rates:\nincidence, prevalence, RHD deaths,\nand other-cause mortality",
    fill=input_fill,
    edge=navy,
    title_fs=12.3,
    sub_fs=8.8,
)

node(
    8.75,
    8.10,
    4.55,
    1.30,
    "Echo screening + SAP scale-up",
    "School-age screening (5–15 years)\nCoverage 5% → 40% by 2030 | 55% RRR\n"
    "pMS(t) = pMS,0 [1 − 0.55 × coverage(t)]",
    fill=intervention_fill,
    edge=intervention_edge,
    title_fs=12.4,
    sub_fs=8.7,
)

node(
    14.60,
    8.10,
    3.75,
    1.20,
    "Tertiary care held fixed",
    "HF management and surgery coverage\nremain at baseline in both scenarios",
    fill=purple_fill,
    edge=purple_edge,
    title_fs=11.8,
    sub_fs=8.9,
)


# -----------------------------------------------------------------------------
# Main states and inputs
# -----------------------------------------------------------------------------
node(
    1.45,
    5.65,
    2.50,
    1.25,
    "Population",
    "Single-year age × sex\nprojection",
    fill=input_fill,
    edge=navy,
    title_fs=13,
)

node(
    5.25,
    5.65,
    2.85,
    1.42,
    "Mild RHD",
    "Asymptomatic / subclinical\nIncident cases enter here\nSurvivors remain in state",
    fill=state_fill,
    edge=navy,
    title_fs=13.3,
    sub_fs=9.0,
)

node(
    10.10,
    5.65,
    2.95,
    1.42,
    "Severe RHD",
    "Heart failure / advanced disease\nNon-surgical survivors remain",
    fill=state_fill,
    edge=navy,
    title_fs=13.3,
    sub_fs=9.0,
)

node(
    14.70,
    5.65,
    2.85,
    1.42,
    "Post-surgery",
    "Valve-intervention survivors\nSurvivors remain in state",
    fill=purple_fill,
    edge=purple_edge,
    title_fs=13.0,
    sub_fs=9.0,
)

node(
    10.20,
    2.15,
    3.10,
    1.25,
    "Death",
    "RHD, operative, or\nother-cause mortality",
    fill=red_fill,
    edge=red_edge,
    title_fs=13.0,
    sub_fs=9.2,
)


# -----------------------------------------------------------------------------
# Data-fed entry mechanisms
# -----------------------------------------------------------------------------
arrow(
    (2.72, 5.65),
    (3.78, 5.65),
    "Annual incident mild RHD\n= incidence rate × population",
    color=navy,
    lw=2.15,
    label_xy=(3.25, 6.16),
    label_fs=9.1,
)

curved_path(
    [
        (2.50, 7.46),
        (3.35, 6.95),
        (3.75, 6.55),
        (4.42, 6.30),
    ],
    label="Prevalence seed in 2026\n96% mild",
    label_xy=(3.52, 7.02),
    color=navy,
    lw=1.65,
    linestyle="--",
    label_fs=8.8,
)

curved_path(
    [
        (2.82, 7.48),
        (5.50, 7.05),
        (7.80, 6.80),
        (9.08, 6.30),
    ],
    label="3% severe",
    label_xy=(6.95, 7.10),
    color=navy,
    lw=1.55,
    linestyle="--",
    label_fs=8.8,
)

curved_path(
    [
        (3.00, 7.52),
        (7.90, 7.35),
        (12.30, 7.10),
        (13.80, 6.30),
    ],
    label="1% post-surgery",
    label_xy=(11.55, 7.35),
    color=navy,
    lw=1.55,
    linestyle="--",
    label_fs=8.8,
)


# -----------------------------------------------------------------------------
# Main disease transitions
# -----------------------------------------------------------------------------
arrow(
    (6.70, 5.65),
    (8.60, 5.65),
    color=intervention_line,
    lw=2.9,
    mutation=19,
)
ax.text(
    7.65,
    6.05,
    "Mild → severe progression",
    ha="center",
    va="center",
    fontsize=9.7,
    color=intervention_edge,
    fontweight="bold",
    zorder=7,
    bbox={
        "boxstyle": "round,pad=0.14",
        "facecolor": white,
        "edgecolor": "none",
        "alpha": 0.96,
    },
)

arrow(
    (11.60, 5.65),
    (13.25, 5.65),
    "Surgery among eligible\nsevere cases",
    color=purple_edge,
    lw=2.15,
    label_xy=(12.42, 6.08),
    label_fs=9.0,
)


# -----------------------------------------------------------------------------
# Intervention pointers
# -----------------------------------------------------------------------------
arrow(
    (8.75, 7.43),
    (7.72, 6.05),
    color=intervention_line,
    lw=1.8,
    linestyle="--",
    mutation=14,
    z=4,
)
arrow(
    (14.60, 7.49),
    (12.50, 6.12),
    color=purple_edge,
    lw=1.6,
    linestyle="--",
    mutation=14,
    z=4,
)


# -----------------------------------------------------------------------------
# Mortality pathways
# -----------------------------------------------------------------------------
curved_path(
    [
        (5.30, 4.92),
        (5.90, 3.85),
        (7.55, 2.80),
        (8.70, 2.28),
    ],
    label="Other-cause mortality",
    label_xy=(6.85, 3.48),
    color=red_edge,
    lw=1.8,
)

curved_path(
    [
        (10.10, 4.92),
        (10.10, 4.05),
        (10.15, 3.12),
        (10.18, 2.78),
    ],
    label="RHD mortality reduced by HF management\n+ other-cause mortality",
    label_xy=(11.45, 3.62),
    color=red_edge,
    lw=1.9,
    label_fs=8.8,
)

curved_path(
    [
        (12.55, 5.42),
        (12.90, 4.25),
        (12.05, 3.08),
        (11.55, 2.60),
    ],
    label="Operative mortality",
    label_xy=(13.15, 4.18),
    color=red_edge,
    lw=1.7,
    linestyle="--",
    label_fs=8.8,
)

curved_path(
    [
        (14.70, 4.92),
        (14.55, 3.75),
        (13.05, 2.75),
        (11.72, 2.30),
    ],
    label="Residual RHD + other-cause mortality",
    label_xy=(14.45, 3.35),
    color=red_edge,
    lw=1.8,
    label_fs=8.8,
)


# -----------------------------------------------------------------------------
# Cohort ageing / next-cycle mechanism
# -----------------------------------------------------------------------------
node(
    2.90,
    2.05,
    4.15,
    1.15,
    "Annual ageing of survivors",
    "After each cycle: age a → a + 1\nAge 100 is an open terminal group",
    fill=cycle_fill,
    edge=cycle_edge,
    title_fs=11.8,
    sub_fs=8.9,
    linewidth=1.6,
)

curved_path(
    [
        (5.00, 4.92),
        (4.55, 4.00),
        (3.82, 3.10),
        (3.35, 2.60),
    ],
    color=cycle_edge,
    lw=1.6,
    linestyle="--",
)
curved_path(
    [
        (9.30, 4.95),
        (7.70, 3.95),
        (5.00, 2.98),
        (4.25, 2.40),
    ],
    color=cycle_edge,
    lw=1.5,
    linestyle="--",
)
curved_path(
    [
        (13.95, 4.95),
        (11.30, 3.70),
        (6.20, 2.80),
        (4.55, 2.23),
    ],
    color=cycle_edge,
    lw=1.5,
    linestyle="--",
)

arrow(
    (1.95, 2.62),
    (1.45, 5.00),
    "Next annual cycle",
    color=cycle_edge,
    lw=1.6,
    rad=-0.18,
    label_xy=(1.05, 3.72),
    linestyle="--",
    mutation=14,
    label_fs=8.8,
)


# -----------------------------------------------------------------------------
# Legend and footer
# -----------------------------------------------------------------------------
legend_y = 0.68

ax.text(
    0.45,
    legend_y,
    "Legend",
    fontsize=10.2,
    fontweight="bold",
    color=text_dark,
    va="center",
)

def legend_swatch(x, fill, edge, label):
    swatch = FancyBboxPatch(
        (x, legend_y - 0.18),
        0.42,
        0.36,
        boxstyle="round,pad=0.02,rounding_size=0.04",
        facecolor=fill,
        edgecolor=edge,
        linewidth=1.25,
    )
    ax.add_patch(swatch)
    ax.text(
        x + 0.55,
        legend_y,
        label,
        fontsize=9.2,
        color=text_dark,
        va="center",
    )

legend_swatch(1.20, input_fill, navy, "Data or population input")
legend_swatch(4.25, state_fill, navy, "Living RHD state")
legend_swatch(
    7.10,
    intervention_fill,
    intervention_edge,
    "Scaled secondary-prevention intervention",
)
legend_swatch(12.10, purple_fill, purple_edge, "Tertiary-care pathway")
legend_swatch(15.20, red_fill, red_edge, "Absorbing death state")

ax.text(
    17.55,
    0.20,
    "Reference: SAP coverage fixed at 5%  |  Scale-up: 5% → 40% (2026–2030)",
    ha="right",
    va="center",
    fontsize=9.0,
    color=muted,
)
ax.text(
    0.45,
    0.20,
    "No ARF state or ARF-to-RHD transition",
    ha="left",
    va="center",
    fontsize=9.0,
    color=muted,
    fontweight="bold",
)


# -----------------------------------------------------------------------------
# Save
# -----------------------------------------------------------------------------
plt.tight_layout(pad=0.8)
fig.savefig(
    OUTPUT_PATH,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)
plt.close(fig)
print(f"Saved PNG: {OUTPUT_PATH}")
