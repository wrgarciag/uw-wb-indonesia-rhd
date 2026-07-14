from pathlib import Path
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
from matplotlib.path import Path as MplPath
import matplotlib.patheffects as pe

# -----------------------------------------------------------------------------
# Output
# -----------------------------------------------------------------------------
output_path = Path("rhd_secondary_prevention_state_transition.png")

# -----------------------------------------------------------------------------
# Figure setup
# -----------------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(18, 10))
ax.set_xlim(0, 18)
ax.set_ylim(0, 10)
ax.axis("off")

# -----------------------------------------------------------------------------
# Colors
# -----------------------------------------------------------------------------
navy = "#2C3E50"
text_dark = "#17202A"
muted = "#607286"
state_fill = "#EAF2F8"
risk_fill = "#F7F9F9"

intervention_fill = "#FFF3CD"
intervention_edge = "#B7791F"
intervention_line = "#C27C0E"

green_fill = "#EAF7EE"
green_edge = "#2E7D32"

purple_fill = "#F1ECF8"
purple_edge = "#6C4A8B"

red_fill = "#FDECEC"
red_edge = "#B83227"

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
    sub_fs=9.7,
    z=5,
):
    """Draw a rounded state or intervention node."""
    patch = FancyBboxPatch(
        (x - w / 2, y - h / 2),
        w,
        h,
        boxstyle="round,pad=0.035,rounding_size=0.14",
        facecolor=fill,
        edgecolor=edge,
        linewidth=1.8,
        zorder=z,
    )
    patch.set_path_effects(
        [
            pe.SimplePatchShadow(
                offset=(2, -2),
                shadow_rgbFace=(0, 0, 0),
                alpha=0.11,
            ),
            pe.Normal(),
        ]
    )
    ax.add_patch(patch)

    ax.text(
        x,
        y + (0.16 if subtitle else 0),
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
            y - 0.26,
            subtitle,
            ha="center",
            va="center",
            fontsize=sub_fs,
            color=muted,
            zorder=z + 1,
            linespacing=1.15,
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
):
    """Draw a straight or curved arrow between two points."""
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
            else (
                (start[0] + end[0]) / 2,
                (start[1] + end[1]) / 2,
            )
        )
        ax.text(
            lx,
            ly,
            label,
            ha="center",
            va="center",
            fontsize=9.8,
            color=color,
            zorder=z + 3,
            linespacing=1.15,
            bbox={
                "boxstyle": "round,pad=0.16",
                "facecolor": white,
                "edgecolor": "none",
                "alpha": 0.94,
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
):
    """Draw a cubic Bezier arrow using four control points."""
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
            fontsize=9.3,
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
    9.6,
    "Rheumatic Heart Disease Secondary-Prevention Model",
    ha="center",
    va="center",
    fontsize=23,
    fontweight="bold",
    color=text_dark,
)

ax.text(
    9,
    9.18,
    (
        "Annual cohort state transitions; secondary prevention reduces "
        "ARF-to-RHD and mild-to-severe progression"
    ),
    ha="center",
    va="center",
    fontsize=12,
    color=muted,
)


# -----------------------------------------------------------------------------
# Intervention callouts
# -----------------------------------------------------------------------------
node(
    5.55,
    8.15,
    3.25,
    1.08,
    "1  SAP after ARF",
    "Coverage 5% → 40%  |  55% RRR\n"
    "p = p₀[1 − 0.55 × coverage]",
    fill=intervention_fill,
    edge=intervention_edge,
    title_fs=12.5,
    sub_fs=9.2,
)

node(
    9.25,
    8.15,
    4.05,
    1.18,
    "2  Echo screening + SAP",
    "Detect mild RHD, then prophylaxis\n"
    "Coverage 5% → 40%  |  55% RRR\n"
    "p = p₀[1 − 0.55 × coverage]",
    fill=intervention_fill,
    edge=intervention_edge,
    title_fs=12.5,
    sub_fs=8.9,
)

node(
    14.1,
    8.15,
    3.15,
    1.05,
    "Tertiary pathways held fixed",
    "HF management and surgery remain\n"
    "at baseline coverage",
    fill=purple_fill,
    edge=purple_edge,
    title_fs=11.7,
    sub_fs=9.0,
)


# -----------------------------------------------------------------------------
# Main health states
# -----------------------------------------------------------------------------
node(
    1.35,
    5.7,
    2.35,
    1.2,
    "Population at risk",
    "Exogenous annual inflows",
    fill=risk_fill,
)

node(
    4.55,
    5.7,
    2.55,
    1.3,
    "Acute rheumatic\nfever (ARF)",
    "New cases each year",
)

node(
    8.05,
    5.7,
    2.55,
    1.35,
    "Mild RHD",
    "Asymptomatic / subclinical\n"
    "Survivors remain in state",
)

node(
    11.55,
    5.7,
    2.7,
    1.35,
    "Severe RHD",
    "Heart failure / advanced disease\n"
    "Non-surgical survivors remain",
)

node(
    15.05,
    5.7,
    2.5,
    1.35,
    "Post-surgery",
    "Valve intervention survivors\n"
    "Survivors remain in state",
    fill=purple_fill,
    edge=purple_edge,
)


# -----------------------------------------------------------------------------
# Resolved and absorbing states
# -----------------------------------------------------------------------------
node(
    4.55,
    1.95,
    2.55,
    1.1,
    "Resolved / no RHD",
    "ARF remission",
    fill=green_fill,
    edge=green_edge,
)

node(
    11.8,
    1.85,
    2.85,
    1.25,
    "Death",
    "ARF, RHD, operative,\nor other-cause",
    fill=red_fill,
    edge=red_edge,
)


# -----------------------------------------------------------------------------
# Main disease pathway
# -----------------------------------------------------------------------------
arrow(
    (2.53, 5.7),
    (3.25, 5.7),
    "Incident ARF",
    label_xy=(2.9, 6.05),
)

arrow(
    (5.83, 5.7),
    (6.76, 5.7),
    color=intervention_line,
    lw=2.8,
)

ax.text(
    6.30,
    6.02,
    "ARF → RHD",
    ha="center",
    va="center",
    fontsize=9.5,
    color=intervention_edge,
    fontweight="bold",
    zorder=7,
    bbox={
        "boxstyle": "round,pad=0.13",
        "facecolor": white,
        "edgecolor": "none",
        "alpha": 0.96,
    },
)

arrow(
    (9.34, 5.7),
    (10.18, 5.7),
    color=intervention_line,
    lw=2.8,
)

ax.text(
    9.76,
    6.02,
    "Mild → severe RHD",
    ha="center",
    va="center",
    fontsize=9.3,
    color=intervention_edge,
    fontweight="bold",
    zorder=7,
    bbox={
        "boxstyle": "round,pad=0.13",
        "facecolor": white,
        "edgecolor": "none",
        "alpha": 0.96,
    },
)

arrow(
    (12.92, 5.7),
    (13.78, 5.7),
    "Surgery among eligible cases",
    color=purple_edge,
    lw=2.1,
    label_xy=(13.35, 6.12),
)


# -----------------------------------------------------------------------------
# Incident asymptomatic RHD inflow
# -----------------------------------------------------------------------------
curved_path(
    [
        (2.45, 5.35),
        (3.7, 4.15),
        (5.65, 4.15),
        (6.82, 5.35),
    ],
    label="Incident asymptomatic RHD",
    label_xy=(4.7, 4.28),
    color=navy,
    lw=1.9,
    z=2,
)


# -----------------------------------------------------------------------------
# Resolution
# -----------------------------------------------------------------------------
arrow(
    (4.55, 5.03),
    (4.55, 2.53),
    "Resolution",
    color=green_edge,
    label_xy=(4.05, 3.75),
)


# -----------------------------------------------------------------------------
# Mortality pathways
# -----------------------------------------------------------------------------
curved_path(
    [
        (5.15, 5.05),
        (6.6, 3.85),
        (8.9, 2.65),
        (10.45, 2.05),
    ],
    label="Acute ARF death",
    label_xy=(7.1, 3.48),
)

curved_path(
    [
        (8.2, 5.02),
        (8.55, 4.0),
        (9.35, 2.95),
        (10.55, 2.18),
    ],
    label="Other-cause death",
    label_xy=(9.05, 3.3),
)

curved_path(
    [
        (11.55, 5.02),
        (11.55, 4.15),
        (11.7, 3.0),
        (11.75, 2.49),
    ],
    label=(
        "RHD + other-cause death\n"
        "HF management lowers RHD death"
    ),
    label_xy=(12.75, 3.72),
)

curved_path(
    [
        (13.35, 5.48),
        (13.7, 4.35),
        (13.2, 3.0),
        (12.82, 2.42),
    ],
    label="Operative mortality",
    label_xy=(14.0, 4.25),
    linestyle="--",
    lw=1.8,
)

curved_path(
    [
        (15.05, 5.02),
        (15.0, 3.8),
        (14.2, 2.75),
        (13.15, 2.15),
    ],
    label="Residual RHD + other-cause death",
    label_xy=(15.05, 3.2),
)


# -----------------------------------------------------------------------------
# Intervention pointers
# -----------------------------------------------------------------------------
arrow(
    (5.55, 7.60),
    (6.15, 6.0),
    color=intervention_line,
    lw=1.8,
    linestyle="--",
    z=4,
    mutation=14,
)

arrow(
    (9.25, 7.55),
    (9.72, 6.0),
    color=intervention_line,
    lw=1.8,
    linestyle="--",
    z=4,
    mutation=14,
)

arrow(
    (14.1, 7.62),
    (13.38, 6.03),
    color=purple_edge,
    lw=1.6,
    linestyle="--",
    z=4,
    mutation=14,
)


# -----------------------------------------------------------------------------
# Legend and footer
# -----------------------------------------------------------------------------
legend_y = 0.65

ax.text(
    0.45,
    legend_y,
    "Legend",
    fontsize=10.3,
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
        fontsize=9.5,
        color=text_dark,
        va="center",
    )


legend_swatch(
    1.2,
    state_fill,
    navy,
    "Disease state",
)

legend_swatch(
    3.25,
    intervention_fill,
    intervention_edge,
    "Scaled secondary-prevention intervention",
)

legend_swatch(
    7.45,
    purple_fill,
    purple_edge,
    "Tertiary-care pathway held at baseline",
)

legend_swatch(
    11.55,
    red_fill,
    red_edge,
    "Absorbing death state",
)

ax.text(
    17.55,
    0.65,
    "Cycle length: 1 year  |  Living-state stocks carry forward",
    ha="right",
    va="center",
    fontsize=9.3,
    color=muted,
)

ax.text(
    9,
    0.18,
    (
        "Model logic from 06_secondary_prevention_model(1).R. "
        "ARF and asymptomatic RHD incidence are exogenous; "
        "secondary prevention changes progression probabilities."
    ),
    ha="center",
    va="center",
    fontsize=9.0,
    color=muted,
)


# -----------------------------------------------------------------------------
# Save output
# -----------------------------------------------------------------------------
plt.tight_layout(pad=0.8)

fig.savefig(
    output_path,
    dpi=300,
    bbox_inches="tight",
    facecolor="white",
)

plt.show()

print(f"Saved PNG: {output_path.resolve()}")
