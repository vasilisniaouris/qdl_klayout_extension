"""
This module contains methods that retrieve coordinates of basic shapes, e.g. rectangle, where the origin (0, 0) is
the same as the centroid of the shape.
"""

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


