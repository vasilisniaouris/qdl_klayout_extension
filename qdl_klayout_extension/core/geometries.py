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
    coords.sort_counterclockwise()
    return coords


def circle_coords(radius: num_ext, num_points: int) -> CoordinatesList:
    """
    Return the coordinates of a circle centered at (0, 0) with the given radius.
    In order to define a circle, we make sure that there is a small overlap between start and end,
    so that the circle does not stay accidentally open.

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
    angles = np.linspace(-np.pi, np.pi, num_points + 1)
    angles = np.append(angles, angles[1] + 2 * np.pi)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    coords = CoordinatesList(x_shape, y_shape)

    return coords


def arc_coords(radius: num_ext, start_angle: num_ext, angle_range: num_ext, num_points: int) -> CoordinatesList:
    """
    Return the coordinates of an arc centered at (0, 0) with the given radius, start angle, angle range, and number of
    points. Angle 0 deg is parallel to the positive x-axis. Rotation is counterclockwise.

    Parameters
    ----------
    radius : int | float | pint.Quantity
        The radius of the arc.
    start_angle : int | float | pint.Quantity
        The start angle of the arc.
    angle_range : int | float | pint.Quantity
        The angle range of the arc
    num_points : int
        The number of points to use to approximate the arc.

    Returns
    -------
    CoordinatesList
        A `CoordinatesList` object containing the x- and y-coordinates of the arc vertices.
    """

    start_angle = qty_in_uu(start_angle) if isinstance(start_angle, Qty) else v_in_uu(start_angle, 'angle')
    angle_range = qty_in_uu(angle_range) if isinstance(angle_range, Qty) else v_in_uu(angle_range, 'angle')
    start_angle = start_angle.to('rad')
    angle_range = angle_range.to('rad')
    end_angle = start_angle + angle_range

    if angle_range.m >= 2 * np.pi:  # If it is more than a full circle, returns full circle.
        return circle_coords(radius, num_points)

    angles = np.linspace(start_angle.m, end_angle.m, num_points)
    x_shape = radius * np.cos(angles)
    y_shape = radius * np.sin(angles)

    coords = CoordinatesList(x_shape, y_shape)
    return coords

