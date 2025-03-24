width, height = 500, 400
colors = {
    "st": "#FF1C1C",
    "ecoli": "#33CC33",
}
font_size = 14
coolwarm_colorscale = [
    [0.0, "rgb(58, 76, 192)"],  # Dark blue
    [0.5, "rgb(221, 221, 221)"],  # White
    [1.0, "rgb(180, 4, 38)"],  # Dark red
]


def style_plot(
    fig,
    marker_size=3,
    top_margin=20,
    left_margin=40,
    right_margin=0,
    buttom_margin=50,
    font_size=font_size,
    line_thickness=3,
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
        # d["error_y"]["thickness"] = line_thickness
    for a in fig["layout"]["annotations"]:
        a["font"]["size"] = font_size
        a["font"]["color"] = "black"
    fig["layout"]["title"]["font"]["size"] = font_size
    fig["layout"]["title"]["font"]["color"] = "black"
    fig["layout"]["legend"]["title"]["font"]["size"] = font_size
    fig["layout"]["legend"]["title"]["font"]["color"] = "black"
    # fig["layout"]["contour"]["colorbar"]["titlefont"]["size"] = font_size

    fig.update_layout(
        margin=dict(l=left_margin, r=right_margin, t=top_margin, b=buttom_margin),
        hoverlabel=dict(font_size=font_size),
    )
    gridline_width = 0.5
    fig.update_yaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.update_xaxes(
        title_standoff=0,
        gridcolor="black",
        zerolinecolor="black",
        gridwidth=gridline_width,
        zerolinewidth=gridline_width,
    )
    fig.for_each_xaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    fig.for_each_yaxis(
        lambda axis: axis.title.update(font=dict(size=font_size, color="black"))
    )
    return fig
