import scipy
from scipy.spatial.distance import pdist, squareform
import numpy as np
import colorsys

import matplotlib
import matplotlib.pyplot as plt
from matplotlib import rc, rcParams
import matplotlib.tri as tri
import matplotlib.patches as patches
import matplotlib.colors as mc
from matplotlib.pyplot import cm


def log_interp1d(xx, yy, kind="linear", **kwargs):
    """Return an interpolating function using Scipy's interp1d for log-spaced data

    Args:
        xx (ndarray): x-axis data
        yy (ndarray): y-axis data
        kind (str, optional): Interpolating method to use. Defaults to 'linear'.

    Returns:
        scipy.interpolate._interpolate.interp1d: interpolating function.
    """

    logx = np.log10(xx)
    logy = np.log10(yy)
    lin_interp = scipy.interpolate.interp1d(logx, logy, kind=kind, **kwargs)
    log_interp = lambda zz: np.power(10.0, lin_interp(np.log10(zz)))
    return log_interp


###########################
fsize = 12
fsize_annotate = 10
std_figsize = (9, 3.75)
std_axes_form = [0.085, 0.14, 0.9, 0.81]


# standard figure creation
def std_fig(ax_form=std_axes_form, figsize=std_figsize, rasterized=True):
    rcparams = {
        "axes.labelsize": fsize,
        "xtick.labelsize": fsize,
        "ytick.labelsize": fsize,
        "figure.figsize": std_figsize,
        "legend.frameon": False,
        "legend.loc": "best",
    }
    plt.rcParams["text.latex.preamble"] = r"\usepackage{amsmath}\usepackage{amssymb}"
    rc("font", **{"family": "serif", "serif": ["Computer Modern"]})
    rc("text", usetex=True)
    rcParams.update(rcparams)
    matplotlib.rcParams["hatch.linewidth"] = 0.3
    fig = plt.figure(figsize=figsize)
    ax = fig.add_axes(ax_form, rasterized=rasterized)
    # ax.patch.set_alpha(0.0)
    return fig, ax


# standard saving function
def std_savefig(fig, path, dpi=500, **kwargs):
    fig.savefig(path, dpi=dpi, **kwargs, bbox_inches="tight")
    if ".pdf" in path:
        fig.savefig(path.replace(".pdf", ".png"), dpi=dpi, **kwargs, bbox_inches="tight")
        fig.savefig(path.replace(".pdf", "_white.png"), dpi=dpi, facecolor="white", **kwargs, bbox_inches="tight")


def std_plot_limits(
    case,
    figsize=std_figsize,
    ax_form=std_axes_form,
    skip_ids=[],
    xrange=(1e-3, 1e2),
    yrange=(1e-10, 1e-1),
    title=None,
    new_labelpos={},
    new_color={},
    new_dash={},
    new_rotation={},
    grid=False,
    annotate=True,
    color_fill=False,
    color_contour=False,
    color_only=[],
    npoints_interp=int(1e5),
    colormap="tab20",
    colormap_range=(0, 1),
    force_label_color=None,
    force_contour_color=None,
    grey_fill_color="lightgrey",
    label_fontsize=6,
    linewidth=0.5,
    alpha=1,
):
    """std_plot_limits Plotting the mixing limits in the standard format

    Parameters
    ----------
    case : HNLimits.Limits()
        The class containing all the limits, data, and additional information on this particular case.
    skip_ids : list, optional
        A list of strings containing the ids of the limits to be skipped, by default []
    xrange : tuple, optional
        the x axis range, by default (1e-3, 1e2)
    yrange : tuple, optional
        the y axis range, by default (1e-10,1e-1)
    title : str, optional
        the title of the plot to be added to the axis, by default None
    new_labelpos : dict, optional
        dictionary with the ids of the limits and the corresponding desired new labels, by default {}
    new_color : dict, optional
        dictionary with the ids of the limits and the correspoinding desired new colors, by default {}
    new_rotation : dict, optional
        dictionary with the ids of the limits and the correspoinding desired new rotations of the label, by default {}
    grid : bool, optional
        whether to add a grid to the plot, by default False
    annotate: bool, optional
        whether to add annotations to the plot, by default True
    color_fill : list, optional
        list with the ids of the limits to be colorfilled, by default {}
    color_only: list, optional
        color only the limits in this list, by default []
    color_contour: list, optional
        color only the limits in this list, by default []
    npoints_interp : int, optional
        number of points to use when drawing the curves, by default int(1e5)
    suffix : str, optional
        suffix of the plot name and file name, by default ''
    colormap : str, optional
        what colormap to iterate over to color the limits, by default 'tab20'
    colormap_range : tuple, optional
        the range of the colormap, by default (0, 1)
    force_label_color : str, optional
        force the label color to be this color, by default None
    force_contour_color : str, optional
        force the contour color to be this color, by default None
    linewidth : float, optional
        the linewidth of the limits, by default 0.5
    alpha : float, optional
        the alpha of the filled regions, by default 1

    Returns
    -------
    fig, ax
        The figure and axis objects.
    """

    fig, ax = std_fig(figsize=figsize, ax_form=ax_form)

    x = np.geomspace(xrange[0], xrange[1], int(npoints_interp))

    labelpos_dic = {}
    for id, limit in case.limits.iterrows():
        if limit.interp_func is not None:
            ilabel, ival = np.nanargmin(limit.interp_func(x)), np.nanmin(limit.interp_func(x))
            if ilabel is None or ival is None:
                print("Could not find label position for ", id)
                labelpos_dic[id] = (None, None)
            else:
                labelpos_dic[id] = (x[ilabel] * 0.9, ival / 2.5)

    ##########
    # find the average value of U, and determine the zorder based on that
    df_order = case.limits.ualpha4.apply(lambda x: np.min(x) if x is not None else 0)

    df_order_in_x = case.limits.m4.apply(lambda x: np.mean(np.log10(x)[~np.isnan(x)]) if x is not None else 1)
    df_order_in_x = (df_order_in_x - np.min(df_order_in_x)) / np.max(df_order_in_x - np.min(df_order_in_x))
    df_order_in_x = df_order_in_x * (colormap_range[1] - colormap_range[0]) + colormap_range[0]

    color_dic = dict(zip(case.limits.index, getattr(cm, colormap)(df_order_in_x)))  # a list of RGB tuples
    dash_dic = dict(zip(case.limits.index, (1 + len(color_dic.keys())) * [(1, (1, 0))]))
    rot_dic = dict(zip(case.limits.index, (1 + len(color_dic.keys())) * [0]))
    if force_label_color is not None:
        label_color = force_label_color
    if force_contour_color is not None:
        contour_color = force_contour_color

    labelpos_dic.update(new_labelpos)
    color_dic.update(new_color)
    dash_dic.update(new_dash)
    rot_dic.update(new_rotation)

    for id, limit in case.limits.iterrows():
        if (id not in skip_ids) & (limit.interp_func is not None):

            limit_zorder = df_order[f"{id}"] / df_order.max()

            if len(color_only) > 0:
                if id in color_only:
                    if force_contour_color is None:
                        contour_color = lighten_color(color_dic[id], 1) if color_contour else "None"
                    fill_color = contour_color if color_fill else grey_fill_color
                    if force_label_color is None:
                        label_color = color_dic[id]
                # else:
                #     contour_color = 'black'
                #     fill_color = lighten_color('lightgrey', 0.3)
                else:
                    if force_contour_color is None:
                        contour_color = "grey" if color_contour else "None"
                    fill_color = grey_fill_color
                    if force_label_color is None:
                        label_color = "black"
            else:
                # contour_color = lighten_color(color_dic[id], 1)
                if force_contour_color is None:
                    contour_color = color_dic[id] if color_contour else "None"
                fill_color = color_dic[id] if color_fill else grey_fill_color
                if force_label_color is None:
                    label_color = color_dic[id]

            # independently of what happened -- color the regions accordingly to color_fill or not at all.

            if limit.year is None:
                dash = (1, (4, 2))
                contour_color = color_dic[id] if color_contour else "None"
            else:
                dash = dash_dic[id]

            if annotate:
                label = rf"\noindent {limit.plot_label}".replace("(", r" \noindent {\tiny \textsc{(").replace(")", r")} }").replace(r"\\", r"\vspace{-2ex} \\ ")
                _ = ax.annotate(label, xy=labelpos_dic[id], xycoords="data", color=label_color, zorder=5, fontsize=label_fontsize, rotation=rot_dic[id])
            # t.set_bbox(dict(facecolor=background_grey, alpha=0.2, edgecolor='None'))

            if (limit.file_top == limit.file_bottom) and limit.m4 is not None and limit.ualpha4 is not None:

                x_ordered, y_ordered = get_ordered_closed_region((limit.m4, limit.ualpha4), logx=False, logy=True)
                ax.fill(x_ordered, y_ordered, edgecolor=contour_color, facecolor="None", linestyle=dash, zorder=3, lw=linewidth)

                # Filling
                if ("cosmo" in id) or (limit.year is None):
                    continue
                else:
                    ax.fill(x_ordered, y_ordered, facecolor=fill_color, edgecolor="None", alpha=alpha, zorder=limit_zorder)

            else:
                ax.plot(x, limit.interp_func(x), color=contour_color, linestyle=dash, zorder=3, lw=linewidth)
                ax.plot(x, limit.interp_func_top(x), color=contour_color, linestyle=dash, zorder=3, lw=linewidth)

                # Filling
                if ("cosmo" in id) or (limit.year is None):
                    # continue
                    ax.fill_between(x, limit.interp_func(x), limit.interp_func_top(x), facecolor=fill_color, edgecolor="None", alpha=0.1, zorder=limit_zorder)
                    ax.fill_between(
                        x, limit.interp_func(x), limit.interp_func_top(x), facecolor="None", edgecolor="black", alpha=1, linestyle=dash, zorder=limit_zorder + 1
                    )

                else:
                    if limit.interp_func_top(x) is not None and limit.interp_func(x) is not None and limit_zorder is not None:
                        ax.fill_between(
                            x, limit.interp_func(x), limit.interp_func_top(x), facecolor=fill_color, edgecolor="None", alpha=alpha, zorder=limit_zorder
                        )
                        # ax.fill_between(
                        #     x,
                        #     limit.interp_func(x),
                        #     limit.interp_func_top(x),
                        #     facecolor="None",
                        #     edgecolor=fill_color,
                        #     alpha=1,
                        #     zorder=limit_zorder + 1,
                        # )

    ax.set_yscale("log")
    ax.set_xscale("log")

    ax.set_ylabel(rf"{case.latexflavor}", fontsize=12)
    ax.set_xlabel(r"$M_N$ (GeV)", fontsize=12)

    ax.tick_params(which="both", bottom=True, top=True, left=True, right=True)
    ax.tick_params(labelbottom=True, labeltop=False, labelleft=True, labelright=False)

    major = np.array([1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 1])
    minor = np.array([2, 3, 4, 5, 6, 7, 8, 9])
    minor = np.array([m * minor for m in major]).flatten()[:-9]
    ax.set_yticks(major)
    ax.set_yticks(minor, minor=True)

    ax.tick_params(which="major", axis="x", length=6, direction="in")
    ax.tick_params(which="minor", axis="x", length=3, direction="in")
    ax.tick_params(which="major", axis="y", length=6, direction="in")
    ax.tick_params(which="minor", axis="y", length=3, direction="in")

    if grid:
        ax.grid(axis="y", which="both", dashes=(6, 1), alpha=0.25, c="black", lw=0.2)
        ax.grid(axis="x", which="both", dashes=(6, 1), alpha=0.25, c="black", lw=0.2)

    ax.set_ylim(*yrange)
    ax.set_xlim(*xrange)
    ax.set_title(title)

    return fig, ax


# https://stackoverflow.com/questions/37765197/darken-or-lighten-a-color-in-matplotlib
def lighten_color(color, amount=0.5):
    """
    Lightens the given color by multiplying (1-luminosity) by the given amount.
    Input can be matplotlib color string, hex string, or RGB tuple.

    Examples:
    >> lighten_color('g', 0.3)
    >> lighten_color('#F034A3', 0.6)
    >> lighten_color((.3,.55,.1), 0.5)
    """
    try:
        c = mc.cnames[color]
    except:
        c = color
    c = colorsys.rgb_to_hls(*mc.to_rgb(c))
    return colorsys.hls_to_rgb(c[0], 1 - amount * (1 - c[1]), c[2])


def double_axes_fig(height=0.5, gap=0.1, axis_base=[0.14, 0.1, 0.80, 0.18], figsize=std_figsize, split_y=False, split_x=False, rasterized=True):
    fig = plt.figure(figsize=figsize)

    if split_y and not split_x:
        axis_base = [0.14, 0.1, 0.80, 0.4 - gap / 2]
        axis_appended = [0.14, 0.5 + gap / 2, 0.80, 0.4 - gap / 2]

    elif not split_y and split_x:
        axis_appended = [0.14, 0.1, 0.4 - gap / 2, 0.8]
        axis_base = [0.14 + 0.4 + gap / 2, 0.1, 0.4 - gap / 2, 0.8]

    else:
        axis_base[-1] = height
        axis_appended = axis_base + np.array([0, height + gap, 0, 1 - 2 * height - gap - axis_base[1] - 0.05])

    ax1 = fig.add_axes(axis_appended, rasterized=rasterized)
    ax2 = fig.add_axes(axis_base, rasterized=rasterized)
    ax1.patch.set_alpha(0.0)
    ax2.patch.set_alpha(0.0)

    return fig, [ax1, ax2]


def data_plot(ax, X, Y, xerr, yerr, zorder=1, label="data"):
    return ax.errorbar(
        X,
        Y,
        yerr=yerr,
        xerr=xerr,
        marker="o",
        markeredgewidth=0.5,
        capsize=1.0,
        markerfacecolor="black",
        markeredgecolor="black",
        ms=2,
        lw=0.0,
        elinewidth=0.8,
        color="black",
        label=label,
        zorder=zorder,
    )


def step_plot(ax, x, y, lw=1, color="red", label="signal", where="post", dashes=(3, 0), zorder=1):
    return ax.step(np.append(x, np.max(x) + x[-1]), np.append(y, 0.0), where=where, lw=lw, dashes=dashes, color=color, label=label, zorder=zorder)


def get_ordered_closed_region(points, logx=False, logy=False):
    x, y = points
    # check for nans
    if np.isnan(points).sum() > 0:
        raise ValueError("NaN's were found in input data. Cannot order the contour.")

    # check for repeated x-entries --
    # this is an error because
    x, mask_diff = np.unique(x, return_index=True)
    y = y[mask_diff]

    if logy:
        if (y == 0).any():
            raise ValueError("y values cannot contain any zeros in log mode.")
        sy = 1  # np.sign(y)
        ssy = 1  # (np.abs(y) < 1) * (-1) + (np.abs(y) > 1) * (1)
        y = ssy * np.log10(y * sy)
    if logx:
        if (x == 0).any():
            raise ValueError("x values cannot contain any zeros in log mode.")
        sx = 1  # np.sign(x)
        ssx = 1  # (x < 1) * (-1) + (x > 1) * (1)
        x = ssx * np.log10(x * sx)

    xmin, ymin = np.min(x), np.min(y)
    x, y = x - xmin, y - ymin

    points = np.array([x, y]).T
    # points_s     = (points - points.mean(0))
    # angles       = np.angle((points_s[:,0] + 1j*points_s[:,1]))
    # points_sort  = points_s[angles.argsort()]
    # points_sort += points.mean(0)

    # if np.isnan(points_sort).sum()>0:
    #     raise ValueError("NaN's were found in sorted points. Cannot order the contour.")
    # # print(points.mean(0))
    # # return points_sort
    # tck, u = splprep(points_sort.T, u=None, s=0.0, per=0, k=1)
    # # u_new = np.linspace(u.min(), u.max(), len(points[:,0]))
    # x_new, y_new = splev(u, tck, der=0)
    # # x_new, y_new = splev(u_new, tck, der=0)
    dist_matrix = squareform(pdist(points))

    # Set diagonal to a large number to avoid self-loop
    np.fill_diagonal(dist_matrix, np.inf)

    # Start from the first point
    current_point = 0
    path = [current_point]

    # Find the nearest neighbor of each point
    while len(path) < len(points):
        # Find the nearest point that is not already in the path
        nearest = np.argmin(dist_matrix[current_point])
        # Add the nearest point to the path
        path.append(nearest)
        # Update the current point
        current_point = nearest
        # Mark the visited point so it's not revisited
        dist_matrix[:, current_point] = np.inf

    # Return the ordered path indices and the corresponding points
    x_new, y_new = points[path].T
    x_new, y_new = x_new + xmin, y_new + ymin

    if logx:
        x_new = sx * 10 ** (ssx * x_new)
    if logy:
        y_new = sy * 10 ** (ssy * y_new)

    return x_new, y_new


################################################
def interp_grid(x, y, z, fine_gridx=False, fine_gridy=False, logx=False, logy=False, method="interpolate", smear_stddev=False):
    # default
    if not fine_gridx:
        fine_gridx = 100
    if not fine_gridy:
        fine_gridy = 100

    # log scale x
    if logx:
        xi = np.logspace(np.min(np.log10(x)), np.max(np.log10(x)), fine_gridx)
    else:
        xi = np.linspace(np.min(x), np.max(x), fine_gridx)

    # log scale y
    if logy:
        y = -np.log(y)
        yi = np.logspace(np.min(np.log10(y)), np.max(np.log10(y)), fine_gridy)

    else:
        yi = np.linspace(np.min(y), np.max(y), fine_gridy)

    Xi, Yi = np.meshgrid(xi, yi)
    if logy:
        Yi = np.exp(-Yi)

    # triangulation
    if method == "triangulation":
        triang = tri.Triangulation(x, y)
        interpolator = tri.LinearTriInterpolator(triang, z)
        Zi = interpolator(Xi, Yi)

    elif method == "interpolate":
        Zi = scipy.interpolate.griddata((x, y), z, (xi[None, :], yi[:, None]), method="linear", rescale=True)
    else:
        print(f"Method {method} not implemented.")

    # gaussian smear -- not recommended
    if smear_stddev:
        Zi = scipy.ndimage.filters.gaussian_filter(Zi, smear_stddev, mode="nearest", order=0, cval=0)

    return Xi, Yi, Zi


def my_histogram(ax, df, var, color="black", label=r"new", density=True, ls="-", var_range=(0, 1), units=1, hatch=""):
    out = ax.hist(
        df[var, ""] * units,
        bins=30,
        range=var_range,
        weights=df["weight",],
        facecolor=color,
        ls=ls,
        edgecolor=color,
        histtype="step",
        density=density,
        lw=1,
        zorder=10,
    )

    out = ax.hist(
        df[var, ""] * units,
        label=label,
        bins=30,
        range=var_range,
        weights=df["weight",],
        facecolor="None",
        hatch=hatch,
        edgecolor=color,
        density=density,
        lw=0.0,
        alpha=1,
    )


def error_histogram(ax, df, var, color="black", label=r"new", density=True, ls="-", var_range=(0, 1), units=1, hatch="", cumulative=False):
    w = df["weight",]
    x = df[var, ""] * units
    # if cumulative:
    #     var_range = (np.min(df[var, '']), np.max(df[var, '']))

    prediction, bin_edges = np.histogram(
        x,
        range=var_range,
        bins=50,
        weights=w,
    )

    errors2 = np.histogram(
        x,
        range=var_range,
        bins=50,
        weights=w**2,
    )[0]
    if cumulative == "cum_sum":
        prediction = np.cumsum(prediction)
        errors2 = np.cumsum(errors2)
    elif cumulative == "cum_sum_prior_to":
        prediction = np.cumsum(prediction[::-1])[::-1]
        errors2 = np.cumsum(errors2[::-1])[::-1]

    area = np.sum(prediction)
    prediction /= area
    errors = np.sqrt(errors2) / area

    # plotting
    ax.plot(bin_edges, np.append(prediction, [0]), ds="steps-post", label=label, color=color, lw=0.8)

    for edge_left, edge_right, pred, err in zip(bin_edges[:-1], bin_edges[1:], prediction, errors):
        ax.add_patch(
            patches.Rectangle(
                (edge_left, pred - err),
                edge_right - edge_left,
                2 * err,
                color=color,
                hatch=hatch,
                fill=False,
                linewidth=0,
                alpha=0.5,
            )
        )
