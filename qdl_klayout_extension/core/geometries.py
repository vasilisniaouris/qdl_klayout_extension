"""
This module contains methods that retrieve coordinates of basic shapes, e.g. rectangle, where the origin (0, 0) is
the same as the centroid of the shape.
"""
import numpy as np

from qdl_klayout_extension.constants import num_ext
from qdl_klayout_extension.core.coordinates import CoordinatesList


def rectangle_coords(width: num_ext, length: num_ext) -> CoordinatesList:
    """
    Return the coordinates of a rectangle centered at (0, 0) with the given width and length in counterclockwise order.

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
        counterclockwise order.
    """
    half_width = width / 2
    half_length = length / 2

    x_shape = [- half_width, half_width, half_width, - half_width]
    y_shape = [- half_length, - half_length, half_length, half_length]

    return CoordinatesList(x_shape, y_shape)


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
    angles = np.linspace(0, 2 * np.pi, num_points)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    return CoordinatesList(x_shape, y_shape)


def arc_coords(radius: num_ext, start_angle: num_ext, end_angle: num_ext, num_points: int) -> CoordinatesList:
    """
    Return the coordinates of an arc centered at (0, 0) with the given radius, start and end angles, and number of
    points.

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
    angles = np.linspace(start_angle, end_angle, num_points)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    return CoordinatesList(x_shape, y_shape)

