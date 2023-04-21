"""
A module containing various miscellaneous methods used throughout this package.
Primarily contains methods for basic geometric calculations, as well as unit conversions between
user and database units (See `constants` for more).
"""

import numpy as np
from pint import Quantity as Qty
from typing import List, Sequence, Tuple

from qdl_klayout_extension.constants import num, multi_num, multi_qty, multi_num_ext, DBU_UM, user_units, num_ext
from qdl_klayout_extension.errors import XYLengthError


def v_in_uu(value: multi_num, physical_type='length') -> multi_qty:
    """
    Convert a scalar value in user-defined units to a Quantity object in the corresponding user-defined units.

    Parameters
    ----------
    value : float | int | Sequence[float | int]
        The scalar value to convert.
    physical_type : str, optional
        The physical type of the scalar value, which must be present in `user_units`. Default is `'length'`.

    Returns
    -------
    pint.Quantity | Sequence[pint.Quantity]
        A Quantity object with the converted value and corresponding user units.

    Raises
    ------
    ValueError
        If `physical_type` is not a valid key in `user_units`.
    """

    if physical_type in user_units.keys():
        return Qty(value, user_units[physical_type])
    else:
        raise ValueError(f"Physical type {physical_type} is not in the list of the accepted types "
                         f"{list(user_units.keys())}")


def qty_in_uu(qty: multi_qty) -> multi_qty:
    """
    Convert a Quantity object to a Quantity object in user-defined units.

    Parameters
    ----------
    qty : pint.Quantity | Sequence[pint.Quantity]
        The quantity to convert.

    Returns
    -------
    pint.Quantity | Sequence[pint.Quantity]
        A Quantity object with the converted value and corresponding unit.

    Raises
    ------
    ValueError
        If the physical type of the input Quantity object is not a valid key in `user_units`.
    """

    physical_type = ', '.join([dim[1:-1] for dim in list(qty.units.dimensionality.keys())])
    if physical_type in user_units.keys():
        return qty.to(user_units[physical_type])
    else:
        raise ValueError(f"Physical type {physical_type} of quantity {qty} is not in the list of the accepted types "
                         f"{list(user_units.keys())}")


def ulu2dbu(v: multi_num_ext) -> int | List[int]:
    """
    Convert a scalar value or Quantity object or a list of those in user-defined length units to database units.

    Parameters
    ----------
    v : float | int | pint.Quantity | Sequence[float | int | pint.Quantity]
        The scalar value or Quantity object or list of those to convert to database units.

    Returns
    -------
    int | List[int]
        The converted quantity in database units, rounded to the nearest integer.
    """

    if not isinstance(v, Qty):
        v = v_in_uu(v)

    rv: Qty = v.to('um') / DBU_UM
    if isinstance(rv.m, float | int):
        return int(round(rv.m))
    else:
        return [int(round(element)) for element in rv.magnitude]


def uau2rad(v: multi_num_ext) -> num | Sequence[num]:
    """
    Convert a scalar value or Quantity object or a list of those in user-defined angle units to radians.

    Parameters
    ----------
    v : float | int | pint.Quantity | Sequence[float | int | pint.Quantity]
        The scalar value or Quantity object to convert to database units.

    Returns
    -------
    float | int | Sequence[float | int]
        The converted quantity in database units, rounded to the nearest integer.
    """

    if not isinstance(v, Qty):
        v = v_in_uu(v, 'angle')

    return v.to('rad').m


def find_centroid(x: Sequence[num_ext] | Qty, y: Sequence[num_ext] | Qty) -> \
        Tuple[num_ext, num_ext]:
    """
    Calculate the centroid of a set of points with coordinates (x, y).

    Parameters
    ----------
    x : Sequence[int | float | pint.Quantity]
        x-coordinates of the points
    y : Sequence[int | float | pint.Quantity]
        y-coordinates of the points

    Returns
    -------
    Tuple[int | float | pint.Quantity, int | float | pint.Quantity]
        (xc, yc) - the x and y coordinates of the centroid

    Raises
    ------
    XYLengthError
        If the length of x and y sequences do not match.
    """

    if len(x) != len(y):
        raise XYLengthError(x, y)

    xc = sum(x) / len(x)
    yc = sum(y) / len(y)

    return xc, yc


def cartesian_to_polar(x: num_ext, y: num_ext) -> Tuple[num_ext, num_ext]:
    """
    Converts Cartesian coordinates to polar coordinates.

    Parameters
    ----------
    x : int | float | pint.Quantity
        x-coordinate of the point.
    y : int | float | pint.Quantity
        y-coordinate of the point.

    Returns
    -------
    Tuple[int | float | pint.Quantity, int | float | pint.Quantity]
        The radius and angle in polar coordinates.
    """
    if isinstance(x, List):
        x = Qty.from_list(x) if isinstance(x[0], Qty) else np.array(x)
    if isinstance(y, List):
        y = Qty.from_list(y) if isinstance(y[0], Qty) else np.array(y)

    radius, angle = np.sqrt(x ** 2 + y ** 2), np.arctan2(y, x)
    return radius, angle


def cartesian_to_polar_from_reference_point(x: Sequence, y: Sequence,
                                            xref: num_ext | None = None, yref: num_ext | None = None):
    """
    Converts multiple Cartesian coordinates to polar coordinates, centered at the centroid.

    Parameters
    ----------
    x : Sequence[int | float | pint.Quantity]
        x-coordinate of the points.
    y : Sequence[int | float | pint.Quantity]
        y-coordinate of the points.
    xref
    yref

    Returns
    -------
    Tuple[List[int | float | pint.Quantity], List[int | float | pint.Quantity]]
        The radii and angles in polar coordinates, centered at the centroid.

    Raises
    ------
    XYLengthError
        If the length of x and y sequences do not match.
    """
    if xref is None:
        xref, yref = find_centroid(x, y)

    x_centered = [xp - xref for xp in x]
    y_centered = [yp - yref for yp in y]
    radius, angle = cartesian_to_polar(x_centered, y_centered)

    return radius, angle


def sort_rotationally(x: Sequence[num_ext] | Qty, y: Sequence[num_ext] | Qty,
                      counterclockwise=False, xref: num_ext | None = None, yref: num_ext | None = None) -> \
        Tuple[Sequence[num_ext], Sequence[num_ext]]:
    """
    Sorts coordinates `x` and `y` in a counterclockwise or clockwise direction about their centroid.

    Parameters
    ----------
    x : Sequence[int | float | pint.Quantity]
        x-coordinate of the points.
    y : Sequence[int | float | pint.Quantity]
        y-coordinate of the points.
    counterclockwise:  bool, optional
        If True, sorts coordinates in counterclockwise direction. Defaults to False.
    xref: int | float | pint.Quantity
        The x-coordinate of the reference point for the converting to polar coordinates.
    yref: int | float | pint.Quantity
        The y-coordinate of the reference point for the converting to polar coordinates.

    Returns
    -------
    Tuple[List[num_ext], List[num_ext]]
        The sorted x-coordinates and y-coordinates.

    Raises
    ------
    XYLengthError
        If the length of x and y sequences do not match.
    """

    radius, angle = cartesian_to_polar_from_reference_point(x, y, xref, yref)

    sorting_list = [(xp, yp, a, r) for xp, yp, a, r in zip(x, y, angle, radius)]
    sorted_coords = sorted(sorting_list, key=lambda tpl: (tpl[2], tpl[3]), reverse=counterclockwise)
    x_sorted, y_sorted, _, _ = map(list, zip(*sorted_coords))

    x_sorted = Qty.from_list(x_sorted) if isinstance(x_sorted[0], Qty) else np.array(x_sorted)
    y_sorted = Qty.from_list(y_sorted) if isinstance(y_sorted[0], Qty) else np.array(y_sorted)

    return x_sorted, y_sorted


def line_intersection(c11: Tuple[num_ext, num_ext], c12: Tuple[num_ext, num_ext],
                      c21: Tuple[num_ext, num_ext], c22: Tuple[num_ext, num_ext]) -> \
        Tuple[num_ext | None, num_ext | None]:
    """
    Calculate the intersection point between two line segments.

    More info:  https://en.wikipedia.org/wiki/Line%E2%80%93line_intersection#Given_two_points_on_each_line_segment

    Parameters
    ----------
    c11 : Tuple[int | float | Quantity, int | float | Quantity]
        The first point of the first line segment.
    c12 : Tuple[int | float | Quantity, int | float | Quantity]
        The second point of the first line segment.
    c21 : Tuple[int | float | Quantity, int | float | Quantity]
        The first point of the second line segment.
    c22 : Tuple[int | float | Quantity, int | float | Quantity]
        The second point of the second line segment.

    Returns
    -------
    Tuple[int | float | Quantity | None, int | float | Quantity | None]
        If the line segments intersect, returns the intersection point as a tuple of num_ext values.
        Otherwise, returns (None, None).
    """

    x1, y1 = c11
    x2, y2 = c12
    x3, y3 = c21
    x4, y4 = c22

    t_num = (x1 - x3) * (y3 - y4) - (y1 - y3) * (x3 - x4)
    u_num = (x1 - x3) * (y1 - y2) - (y1 - y3) * (x1 - x2)
    denominator = (x1 - x2) * (y3 - y4) - (y1 - y2) * (x3 - x4)
    t = t_num / denominator
    u = u_num / denominator

    if 0 <= t <= 1 and 0 <= u <= 1:
        x = x1 + t * (x2 - x1)
        y = y1 + t * (y2 - y1)
        return x, y

    return None, None


