width, height = 150, 150
colors = {
    "blue": "#000080",
    "ct": "#7570B3",
    "oa": "#D95F02",
    "ms": "#E6AB02",
    "at": "#1B9E77",
    "Spent media Ct": "#1B9E77",
    "Spent media Oa": "#E7298A",
    "H20": "gray",
    "st": "#de5959",
    "st_error": "#f1a2a2",
    "ecoli": "#50a250",
    "ecoli_error": "#a0d0a0",
}

colors_heatmap = [
    [0.0, "#7570B3"],  # deep blue
    [0.5, "white"],
    [1.0, "#D95F02"],  # deep red
]

colors_metabolites = {
    "Nucleotide related": "#1f77b4",
    "Carbohydrates": "#ff7f0e",
    "Fatty Acids": "#2ca02c",
    "Amino Acids": "#9467bd",
    "Organic Acids": "#8c564b",
    "Coenzymes": "#e377c2",
    "Others": "#7f7f7f",
}


def style_plot(
    fig,
    marker_size=3,
    top_margin=10,
    left_margin=10,
    right_margin=10,
    buttom_margin=10,
    font_size=14,
    line_thickness=1.5,
):
    """Style function for figures setting fot size and true black color."""
    fig.update_layout(
        {
            "plot_bgcolor": "#FFFFFF",
            "paper_bgcolor": "#FFFFFF",
        },
        font={"size": font_size, "color": "black"},
    )
    for d in fig["data"]:
        try:
            d["marker"]["size"] = marker_size
        except KeyError:
            pass
        try:
            d["line"]["width"] = line_thickness
        except KeyError:
            pass
        try:
            d["error_y"]["thickness"] = line_thickness / 2
        except KeyError:
            pass
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = font_size
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = font_size
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = font_size
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"

    fig.update_layout(
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=buttom_margin),
        hoverlabel=dict(font_size=font_size),
    )
    gridline_width = 0.2
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="gray",
        zeroline=False,
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="gray",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=0.5,
        showline=True,
        mirror=True,
        linecolor="black",
        linewidth=0.5,
        zeroline=False,
        tickcolor="black",
        tickwidth=0.5,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig
