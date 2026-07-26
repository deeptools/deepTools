import matplotlib
matplotlib.use("Agg", force=True)

# Central matplotlib configuration:
matplotlib.rcParams.update({

    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "font.size": 8.0,
    "figure.dpi": 100,
    "savefig.dpi": 200,

})
