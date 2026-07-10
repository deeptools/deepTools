import matplotlib

# Backend should be set before importing pyplot
matplotlib.use("Agg")

# Central matplotlib configuration
matplotlib.rcParams.update({

    # Output compatibility
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",

    # Fonts
    "font.size": 8.0,

    # Figure defaults
    "figure.dpi": 100,
    "savefig.dpi": 200,

    # Layout
    #"figure.autolayout": False,

    # Lines
    #"lines.linewidth": 1.0,

    # Tick appearance
    #"xtick.direction": "out",
    #"ytick.direction": "out",

})


