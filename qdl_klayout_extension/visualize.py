"""
A module that contains functions that use matplotlib to plot the coordinates of the KLayout patterns.

Attributes
----------
default_fill_kwargs: Dict[str, Any]
    Contains the default keyword arguments for the matplotlib.pyplot.fill function.
    These include the color, alpha, and hatch parameters.
default_fill_kwargs: Dict[str, Any]
    Contains the default keyword arguments for the matplotlib.pyplot.plot function.
    These include the color, and marker parameters.
"""

from matplotlib import pyplot as plt

from pint import Quantity as Qty
from qdl_klayout_extension.core.coordinates import CoordinatesList

default_fill_kwargs = {'color': 'crimson', 'alpha': 0.3, 'hatch': 'XX'}
default_edge_kwargs = {'color': 'crimson', 'marker': 'o'}


def get_plot_kwargs(edge_kwargs, fill_kwargs):
    """
    Returns keyword arguments for matplotlib.pyplot.plot and matplotlib.pyplot.fill functions.

    Parameters
    ----------
    edge_kwargs : dict or None, optional
        Keyword arguments for matplotlib.pyplot.plot function, by default None.
    fill_kwargs : dict or None, optional
        Keyword arguments for matplotlib.pyplot.fill function, by default None.

    Returns
    -------
    tuple
        A tuple containing the keyword arguments for matplotlib.pyplot.plot and matplotlib.pyplot.fill functions.
    """
    edge_kwargs = {} if edge_kwargs is None else edge_kwargs
    fill_kwargs = {} if fill_kwargs is None else edge_kwargs

    for key in default_edge_kwargs.keys():
        if key not in edge_kwargs.keys():
            edge_kwargs[key] = default_edge_kwargs[key]

    for key in default_fill_kwargs.keys():
        if key not in fill_kwargs.keys():
            fill_kwargs[key] = default_fill_kwargs[key]

    return edge_kwargs, fill_kwargs


def plot_coords(coords: CoordinatesList, edge_kwargs: dict | None = None,
                fill_kwargs: dict | None = None, units='uu'):
    """
    Plot the coordinates of a shape using matplotlib.pyplot.

    Parameters
    ----------
    coords : CoordinatesList
        The list of coordinates of the shape.
    edge_kwargs : dict or None, optional
        Keyword arguments for the edge of the shape, by default None.
    fill_kwargs : dict or None, optional
        Keyword arguments for the fill of the shape, by default None.
    units : str, optional
        The units of the coordinates, can be either 'uu' or 'dbu', by default 'uu'.

    Raises
    ------
    ValueError
        If the input units are not 'uu' or 'dbu'.
    """
    if units not in ['uu', 'dbu']:
        raise ValueError(f"Only acceptable units are {['uu', 'dbu']}.")

    edge_kwargs, fill_kwargs = get_plot_kwargs(edge_kwargs, fill_kwargs)

    x = list(getattr(coords, f'x_{units}'))
    y = list(getattr(coords, f'y_{units}'))
    x = Qty.from_list(x + [x[0]]).m if isinstance(x[0], Qty) else x + [x[0]]
    y = Qty.from_list(y + [y[0]]).m if isinstance(y[0], Qty) else y + [y[0]]

    plt.plot(x, y, **edge_kwargs)
    plt.fill(x, y, **fill_kwargs)

    ax: plt.Axes = plt.gca()
    ax.set_aspect('equal')


def plot_pattern(pattern_array: CoordinatesList, coords: CoordinatesList, edge_kwargs: dict | None = None,
                 fill_kwargs: dict | None = None, units='uu'):
    """
    Plot the repeated pattern of coordinates of a shape using matplotlib.pyplot.

    Parameters
    ----------
    pattern_array: CoordinatesList
        The list of coordinates of a pattern, as defined by Pattern
    coords : CoordinatesList
        The list of coordinates of the shape to be repeatedly plotted.
    edge_kwargs : dict or None, optional
        Keyword arguments for the edge of the shape, by default None.
    fill_kwargs : dict or None, optional
        Keyword arguments for the fill of the shape, by default None.
    units : str, optional
        The units of the coordinates, can be either 'uu' or 'dbu', by default 'uu'.

    Raises
    ------
    ValueError
        If the input units are not 'uu' or 'dbu'.
    """

    for pa_coord in pattern_array:
        trans_coords = coords.get_translated(pa_coord)
        plot_coords(trans_coords, edge_kwargs, fill_kwargs, units)

# TODO: define new functions to plot shapes (simple polyogns and paths) and patterns rather than just coords.
