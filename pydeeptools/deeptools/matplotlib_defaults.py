import matplotlib

# Backend must be set before importing pyplot
matplotlib.use("Agg")

matplotlib.rcParams.update({

    # --------------------------------------------------
    # Output compatibility
    # --------------------------------------------------
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",

    # --------------------------------------------------
    # Fonts
    # --------------------------------------------------
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans"],
    "font.size": 8.0,

    # Stabilize FreeType rendering
    "text.hinting": "none",

    # --------------------------------------------------
    # Figure defaults
    # --------------------------------------------------
    "figure.dpi": 100,
    "savefig.dpi": 100,
    "figure.figsize": [6.4, 4.8],

    "figure.facecolor": "white",
    "axes.facecolor": "white",
    "savefig.facecolor": "white",

    # Disable automatic layout changes
    "figure.autolayout": False,
    "figure.constrained_layout.use": False,

    # --------------------------------------------------
    # Axes / labels
    # --------------------------------------------------
    "axes.titlesize": 8.0,
    "axes.labelsize": 8.0,

    "xtick.labelsize": 8.0,
    "ytick.labelsize": 8.0,

    "axes.autolimit_mode": "data",
    "axes.xmargin": 0,
    "axes.ymargin": 0,

    # --------------------------------------------------
    # Lines and paths
    # --------------------------------------------------
    "lines.linewidth": 1.5,

    "path.simplify": False,
    "path.simplify_threshold": 0.0,

    # --------------------------------------------------
    # Images
    # --------------------------------------------------
    "image.interpolation": "nearest",
    "image.resample": False,

    # --------------------------------------------------
    # Legend
    # --------------------------------------------------
    "legend.fontsize": 8.0,
})