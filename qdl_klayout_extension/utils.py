"""
A module containing various miscellaneous methods used throughout this package.
Primarily contains methods for basic geometric calculations, as well as unit conversions between
user and database units (See `constants` for more).
"""

import numpy as np

from collections import deque
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
    if physical_type == '' and qty.units in ['degree', 'radian']:
        physical_type = 'angle'
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


def dbu2ulu(v: int | List[int]) -> int | List[int]:
    """
    Convert a scalar value or Quantity object or a list of those in database units to user-defined length units.

    Parameters
    ----------
    v : int | List[int]
        The scalar value or Quantity object or list of those to convert to user-defined length units.

    Returns
    -------
    pint.Quantity | Sequence[pint.Quantity]
        The converted quantity in user-defined length units.
    """

    qty = Qty(v, 'um') * DBU_UM
    return qty.to(user_units['length'])


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


def cartesian_to_polar_from_reference_point(x: Sequence | num_ext, y: Sequence | num_ext,
                                            xref: num_ext | None = None, yref: num_ext | None = None):
    """
    Converts multiple Cartesian coordinates to polar coordinates, centered at the centroid.

    Parameters
    ----------
    x: Sequence[int | float | pint.Quantity] | int | float | pint.Quantity
        x-coordinate of the points.
    y: Sequence[int | float | pint.Quantity] | int | float | pint.Quantity
        y-coordinate of the points.
    xref: int | float | pint.Quantity
        The x-coordinate of the reference point for the converting to polar coordinates.
    yref: int | float | pint.Quantity
        The y-coordinate of the reference point for the converting to polar coordinates.

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

    x_centered = [xp - xref for xp in x] if isinstance(x, Sequence) else x - xref
    y_centered = [yp - yref for yp in y] if isinstance(y, Sequence) else y - yref
    radius, angle = cartesian_to_polar(x_centered, y_centered)

    return radius, angle


def sort_rotationally(x: Sequence[num_ext] | Qty, y: Sequence[num_ext] | Qty,
                      counterclockwise=False, xref: num_ext | None = None, yref: num_ext | None = None,
                      start_angle: num_ext | None = None) -> \
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
    start_angle: int | float | pint.Quantity | None
        The angle to start sorting coordinates. Defaults to None.

    Returns
    -------
    Tuple[List[num_ext], List[num_ext]]
        The sorted x-coordinates and y-coordinates.

    Raises
    ------
    XYLengthError
        If the length of x and y sequences do not match.
    """

    radii, angles = cartesian_to_polar_from_reference_point(x, y, xref, yref)
    angles = shift_angles(Qty(angles, 'rad'), start_angle)

    sorting_list = [(xp, yp, a, r) for xp, yp, a, r in zip(x, y, angles, radii)]
    sorted_coords = sorted(sorting_list, key=lambda tpl: (tpl[2], tpl[3]), reverse=(not counterclockwise))
    x_sorted, y_sorted, _, _ = map(list, zip(*sorted_coords))

    x_sorted = Qty.from_list(x_sorted) if isinstance(x_sorted[0], Qty) else np.array(x_sorted)
    y_sorted = Qty.from_list(y_sorted) if isinstance(y_sorted[0], Qty) else np.array(y_sorted)

    return x_sorted, y_sorted


def shift_angles(angles: multi_num_ext, start_angle: num_ext) -> Qty:
    """
    Shifts the given angles so that they start at the given start_angle.

    Parameters
    ----------
    angles : Sequence[int | float | pint.Quantity]
        The angles of the points.
    start_angle: int | float | pint.Quantity
        The angle to start sorting coordinates.

    Returns
    -------
    Sequence[int | float | pint.Quantity]
        The shifted angles.
    """

    angles = Qty.from_list(angles).to('rad').m if isinstance(angles[0], Qty) else v_in_uu(angles, 'angle').to('rad').m

    if start_angle is None:
        start_angle = angles[0]

    start_angle = start_angle.to('rad').m if isinstance(start_angle, Qty) else v_in_uu(start_angle, 'angle').to('rad').m
    start_angle = start_angle % (2 * np.pi)

    angles[angles < 0] += 2 * np.pi
    angles -= start_angle
    angles[angles < 0] += 2 * np.pi

    return Qty(angles, 'rad')


def is_sorted_rotationally(x: Sequence[num_ext] | Qty, y: Sequence[num_ext] | Qty,
                           xref: num_ext | None = None, yref: num_ext | None = None,
                           start_angle: num_ext | None = None):
    """
    Checks if the given x and y coordinates are sorted in a counterclockwise or clockwise direction.
    To check this, we employ a method where we rotate the list of angles by one element at a time, and
    try to see if these angles monotonically increase or decrease. If none of the rotations provide a monotonically
    ordered list, then the list is not sorted rotationally.

    Providing a starting angle helps make this process faster, since the for loop for the list rotation happen
    fewer times.

    Parameters
    ----------
    x : Sequence[int | float | pint.Quantity]
        x-coordinate of the points.
    y : Sequence[int | float | pint.Quantity]
        y-coordinate of the points.
    xref: int | float | pint.Quantity
        The x-coordinate of the reference point for the converting to polar coordinates.
    yref: int | float | pint.Quantity
        The y-coordinate of the reference point for the converting to polar coordinates.
    start_angle: int | float | pint.Quantity | None
        The angle to start sorting coordinates. Defaults to None.
    Returns
    -------
    int
        -1 if the radii and angles are sorted in a clockwise direction.
        0 if the radii and angles are not sorted.
        1 if the radii and angles are sorted in a counterclockwise direction.
    """

    radii, angles = cartesian_to_polar_from_reference_point(x, y, xref, yref)

    if start_angle is not None:
        angles = shift_angles(angles, start_angle)

    angles = np.around(angles, 13)

    # We find all the possible iterations of angles and check if they are sorted.
    # This alleviates us from the pain of dealing with 0-360 extrema.
    original_angles = angles
    for i in range(len(original_angles)):
        shifted_angles = deque(original_angles)
        shifted_angles.rotate(-i)
        angles = type(original_angles)(shifted_angles)

        if all((angles[i], radii[i]) >= (angles[i-1], radii[i-1]) for i in range(1, len(angles))):
            return 1  # counterclockwise
        elif all((angles[i], radii[i]) <= (angles[i-1], radii[i-1]) for i in range(1, len(angles))):
            return -1  # clockwise

    return 0  # if nothing works, the list is not sorted.


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


def line_segment_length(c1: Tuple[num_ext, num_ext], c2: Tuple[num_ext, num_ext]):
    """
    Calculate the length of a line segment.

    Parameters
    ----------
    c1 : Tuple[int | float | Quantity, int | float | Quantity]
        The first point of the line segment.
    c2 : Tuple[int | float | Quantity, int | float | Quantity]
        The second point of the line segment.

    Returns
    -------
    float | Quantity
        The length of the line segment.
    """
    x1, y1 = c1
    x2, y2 = c2

    return np.sqrt((x2 - x1) ** 2 + (y2 - y1) ** 2)


def point_line_segment_distance(c_point: Tuple[num_ext, num_ext], c1_line: Tuple[num_ext, num_ext],
                                c2_line: Tuple[num_ext, num_ext]):
    """
    Calculate the distance from a point to a line segment.

    More info in https://en.wikipedia.org/wiki/Distance_from_a_point_to_a_line#Line_defined_by_two_points

    Parameters
    ----------
    c_point : Tuple[int | float | Quantity, int | float | Quantity]
        The point of interest
    c1_line : Tuple[int | float | Quantity, int | float | Quantity]
        The first point of the line segment.
    c2_line : Tuple[int | float | Quantity, int | float | Quantity]
        The second point of the line segment.

    Returns
    -------
    float | Quantity
        The distance from the point to the line segment.

    """
    x0, y0 = c_point
    x1, y1 = c1_line
    x2, y2 = c2_line

    x21 = x2 - x1
    y21 = y2 - y1
    x10 = x1 - x0
    y10 = y1 - y0

    line_length = line_segment_length(c1_line, c2_line)
    numerator = np.abs(x21 * y10 - x10 * y21)
    return numerator / line_length

