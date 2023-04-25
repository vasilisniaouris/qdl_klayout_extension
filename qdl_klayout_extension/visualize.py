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
import warnings
from typing import Tuple, Dict, List, Pattern

import numpy as np
from matplotlib.lines import Line2D
from matplotlib import pyplot as plt

from pint import Quantity as Qty
from qdl_klayout_extension.core.coordinates import CoordinatesList
from qdl_klayout_extension.core.shapes import SimplePath, SimplePolygon, Shape

default_fill_kwargs = {'color': 'crimson', 'alpha': 0.3, 'hatch': 'XX'}
default_edge_kwargs = {'color': 'crimson', 'marker': 'o'}
plot_styles = ["polygon", "path"]


class Line2DDataUnits(Line2D):
    """
    A class that inherits from matplotlib.lines.Line2D were linewidth is instead given in data units.
    We use this class to easily plot SimplePath-type shapes.

    With this class definition there may be an issue with the size of the legend's linewidth.

    Found on `Stackoverflow <https://stackoverflow.com/questions/19394505/expand-the-line-with-specified-width-in-data-unit/42972469#42972469>`_.
    """
    def __init__(self, *args, **kwargs):
        _lw_data = kwargs.pop("linewidth", 1)
        super().__init__(*args, **kwargs)
        self._lw_data = _lw_data

    def _get_lw(self):
        if self.axes is not None:
            ppd = 72./self.axes.figure.dpi
            trans = self.axes.transData.transform
            return ((trans((1, self._lw_data))-trans((0, 0)))*ppd)[1]
        else:
            return 1

    def _set_lw(self, lw):
        self._lw_data = lw

    _linewidth = property(_get_lw, _set_lw)


def get_plot_kwargs(edge_kwargs, fill_kwargs) -> Tuple[Dict, Dict]:
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
    Tuple[Dict, Dict]
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


def get_xy_from_coords(coords: CoordinatesList, units='uu', closed_shape: bool = True) -> Tuple[np.ndarray, np.ndarray]:
    """
    Get the x and y coordinates from a CoordinatesList.

    Parameters
    ----------
    coords : CoordinatesList
        The list of coordinates of the shape.
    units : str, optional
        The units of the coordinates, can be either 'uu' or 'dbu', by default 'uu'.
    closed_shape: bool
        If the shape is closed (e.g. a polygon), add the start point to the end to close the shape.

    Returns
    -------
    Tuple[np.ndarray, np.ndarray]
        A tuple containing the x and y coordinates of the shape in the given unit.

    Raises
    ------
    ValueError
        If the input units are not 'uu' or 'dbu'.
    """

    if units not in ['uu', 'dbu']:
        raise ValueError(f"Only acceptable units are {['uu', 'dbu']}.")

    x = list(getattr(coords, f'x_{units}'))
    y = list(getattr(coords, f'y_{units}'))
    if closed_shape:
        x = Qty.from_list(x + [x[0]]).m if isinstance(x[0], Qty) else x + [x[0]]
        y = Qty.from_list(y + [y[0]]).m if isinstance(y[0], Qty) else y + [y[0]]
    else:
        x = Qty.from_list(x).m if isinstance(x[0], Qty) else x
        y = Qty.from_list(y).m if isinstance(y[0], Qty) else y

    return x, y


def plot_polygon(polygon: SimplePolygon, edge_kwargs: dict | None = None,
                 fill_kwargs: dict | None = None, units='uu'):
    """
    Plot the coordinates of a shape using matplotlib.pyplot.

    Parameters
    ----------
    polygon : SimplePolygon
        The polygon object to be plotted.
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
    edge_kwargs, fill_kwargs = get_plot_kwargs(edge_kwargs, fill_kwargs)
    x, y = get_xy_from_coords(polygon.coords, units)

    plt.plot(x, y, **edge_kwargs)
    plt.fill(x, y, **fill_kwargs)

    ax: plt.Axes = plt.gca()
    ax.set_aspect('equal')


def plot_path(path: SimplePath, edge_kwargs: dict | None = None, fill_kwargs: dict | None = None, units='uu'):
    """
    Plot a SimplePath using matplotlib.pyplot. Takes into account the end-cap style and size.

    Parameters
    ----------
    path : SimplePath
        The path object to be plotted.
    edge_kwargs : dict or None, optional
        Keyword arguments for the edge of the path, by default None.
    fill_kwargs : dict or None, optional
        Keyword arguments for the fill of the path, by default None.
    units : str, optional
        The units of the coordinates, can be either 'uu' or 'dbu', by default 'uu'.

    Raises
    ------
    ValueError
        If the input units are not 'uu' or 'dbu'.
    """
    edge_kwargs, fill_kwargs = get_plot_kwargs(edge_kwargs, fill_kwargs)
    fill_kwargs.pop('hatch')
    x, y = get_xy_from_coords(path.coords, units, closed_shape=False)
    width = getattr(path, f'width_{units}')
    width = width.m if isinstance(width, Qty) else width

    if path.rounded_end_cap:
        # If the end caps are round, an automatic half width will be added on both ends by matplotlib
        cap_style = 'round'
    else:
        # If the end caps are flat and have non-zero length, we need to add new points to the beginning and the end of
        #  the path coordinates.
        cap_style = 'butt'
        c_init, c_final = path.get_end_cap_coords()
        if path.initial_end_cap_uu > 0:
            np.insert(x, 0, c_init.x_uu)
            np.insert(y, 0, c_init.y_uu)
        if path.final_end_cap_uu > 0:
            np.append(x, c_final.x_uu)
            np.append(y, c_final.y_uu)

    line_fill = Line2DDataUnits(x, y, linewidth=width, **fill_kwargs, solid_capstyle=cap_style,
                                solid_joinstyle='miter')

    ax: plt.Axes = plt.gca()
    ax.add_line(line_fill)
    ax.autoscale(None)  # to update the plot limits, if autoscale is on
    ax.set_aspect('equal')


def plot_shape(shape: Shape | np.ndarray[Shape] | List[List[Shape]], edge_kwargs: dict | None = None,
               fill_kwargs: dict | None = None, units='uu'):
    """
    Plot the coordinates of a shape using matplotlib.pyplot. Different classes of Shapes will be plotted differently.

    Parameters
    ----------
    shape : shape: Shape | np.ndarray[Shape] | List[List[Shape]]
        A shape to be plotted.
        If an array, it must be a 2D array of the same size as the pattern.
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

    if isinstance(shape, SimplePath):
        plot_path(shape, edge_kwargs, fill_kwargs, units)
    elif isinstance(shape, SimplePolygon):
        plot_polygon(shape, edge_kwargs, fill_kwargs, units)
    else:
        warnings.warn(f"Shape {shape}, of type {type(shape)} is not recognized. Nothing will be plotted.")


def plot_pattern(pattern: Pattern, shape: Shape | np.ndarray[Shape] | List[List[Shape]],
                 edge_kwargs: dict | None = None,
                 fill_kwargs: dict | None = None, units='uu'):
    """
    Plot the repeated pattern of coordinates of a shape using matplotlib.pyplot.

    Parameters
    ----------
    pattern: Pattern
        The pattern that will provide the pattern array to be plotted.
    shape : shape: Shape | np.ndarray[Shape] | List[List[Shape]]
        A shape to be repeatedly plotted or an array of shapes to be plotted once each.
        If an array, it must be a 2D array of the same size as the pattern.
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

    pattern_array = pattern.get_pattern_array()

    if isinstance(shape, Shape):
        shapes = [[shape] * pattern.n_y] * pattern.n_x
    else:
        if np.shape(shape) != np.shape(pattern_array):
            raise ValueError(f'The pattern({np.shape(pattern_array)}) and shape({np.shape(shape)}) '
                             f'arrays are not of the same size.')
        shapes = shape

    for i in range(pattern.n_x):
        pattern_column = pattern_array[i]
        shapes_column = shapes[i]
        for j in range(pattern.n_y):
            pattern_coords = pattern_column[j]
            sp = shapes_column[j]
            translated_shape = sp.get_translated(pattern_coords)
            plot_shape(translated_shape, edge_kwargs, fill_kwargs, units)
