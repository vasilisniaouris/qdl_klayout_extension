"""
This module contains methods that retrieve coordinates of basic shapes, e.g. rectangle, where the origin (0, 0) is
the same as the centroid of the shape.
"""
import numpy as np

from pint import Quantity as Qty

from qdl_klayout_extension.constants import num_ext
from qdl_klayout_extension.core.coordinates import CoordinatesList
from qdl_klayout_extension.utils import qty_in_uu, v_in_uu


def rectangle_coords(width: num_ext, length: num_ext) -> CoordinatesList:
    """
    Return the coordinates of a rectangle centered at (0, 0) with the given width and length in clockwise order.

    Parameters
    ----------
    width : int | float | pint.Quantity
        The width of the rectangle.
    length : int | float | pint.Quantity
        The length of the rectangle.

    Returns
    -------
    CoordinatesList
        A `CoordinatesList` object containing the x- and y-coordinates of the rectangle vertices in
        clockwise order.
    """
    half_width = width / 2
    half_length = length / 2

    x_shape = [- half_width, half_width, half_width, - half_width]
    y_shape = [- half_length, - half_length, half_length, half_length]

    coords = CoordinatesList(x_shape, y_shape)
    coords.sort_clockwise()
    return coords


def circle_coords(radius: num_ext, num_points: int) -> CoordinatesList:
    """
    Return the coordinates of a circle centered at (0, 0) with the given radius.

    Parameters
    ----------
    radius : int | float | pint.Quantity
        The radius of the circle.
    num_points : int
        The number of points to use to approximate the circle.

    Returns
    -------
    CoordinatesList
        A `CoordinatesList` object containing the x- and y-coordinates of the circle vertices.
    """
    angles = np.linspace(-np.pi, np.pi, num_points)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    coords = CoordinatesList(x_shape, y_shape)
    coords.sort_clockwise()
    return coords


# TODO: Fix sorting issue with start and end angles.
def arc_coords(radius: num_ext, start_angle: num_ext, end_angle: num_ext, num_points: int) -> CoordinatesList:
    """
    Return the coordinates of an arc centered at (0, 0) with the given radius, start and end angles, and number of
    points. Angle 0 deg is parallel to the y-axis. Rotation is counterclockwise

    Parameters
    ----------
    radius : int | float | pint.Quantity
        The radius of the arc.
    start_angle : int | float | pint.Quantity
        The start angle of the arc.
    end_angle : int | float | pint.Quantity
        The end angle of the arc.
    num_points : int
        The number of points to use to approximate the arc.

    Returns
    -------
    CoordinatesList
        A `CoordinatesList` object containing the x- and y-coordinates of the arc vertices.
    """
    start_angle = qty_in_uu(start_angle) if isinstance(start_angle, Qty) else v_in_uu(start_angle, 'angle')
    end_angle = qty_in_uu(end_angle) if isinstance(end_angle, Qty) else v_in_uu(end_angle, 'angle')

    angles = np.linspace(start_angle.to('rad').m, end_angle.to('rad').m, num_points)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    coords = CoordinatesList(x_shape, y_shape)
    coords.sort_clockwise()
    return coords

