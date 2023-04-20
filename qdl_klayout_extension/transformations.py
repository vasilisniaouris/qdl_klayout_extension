"""
A module containing different methods that can implement various types of coordinate transformations.

More info on https://en.wikipedia.org/wiki/Transformation_matrix.
"""


import math
from typing import Iterable
from pint import Quantity as Qty

import numpy as np

from qdl_klayout_extension.constants import num_ext, num, multi_num_ext
from qdl_klayout_extension.utils import uau2rad


def rotation_matrix(angle: num_ext, counterclockwise=True):
    """
    Returns a rotation matrix that rotates a coordinate around the origin (0, 0) by a given angle.

    Parameters
    ----------
    angle : int | float
        The angle in degrees to rotate the coordinate.
    counterclockwise : bool, optional
        If True (default), rotates the coordinate counterclockwise. If False, rotates clockwise.

    Returns
    -------
    np.ndarray
        A 2x2 numpy array representing the rotation matrix.
    """
    angle = uau2rad(angle) if counterclockwise else -uau2rad(angle)
    c = math.cos(angle)
    s = math.sin(angle)

    rot_mat = np.array([[c, - s], [s, c]])
    return rot_mat


def stretching_matrix(x_stretch: num, y_stretch: num):
    """
    Returns a stretching matrix that enlarges or compresses a shape in the x and y axes with
    the shape's centroid at (0, 0).

    Parameters
    ----------
    x_stretch : int | float
        The stretch factor in the x-axis.
    y_stretch : int | float
        The stretch factor in the y-axis.

    Returns
    -------
    np.ndarray
        A 2x2 numpy array representing the stretching matrix.
    """
    return np.array([[x_stretch, 1], [1, y_stretch]])


def squeezing_matrix(squeeze_value: num):
    """
    Returns a squeezing matrix that "squeezes" a shape by enlarging it in one direction and
    compressing it in the other, with the shape's centroid at (0, 0). Squeezing conserves the shape's area.

    Parameters
    ----------
    squeeze_value : int | float
        The factor by which to squeeze the shape. A value larger than 1 enlarges the shape in the x direction.
        A value smaller than 1 enlarges the shape in the y direction.

    Returns
    -------
    np.ndarray
        A 2x2 numpy array representing the squeezing matrix.
    """
    return np.array([[squeeze_value, 1], [1, 1/squeeze_value]])


def reflection_matrix(lx: num, ly: num):
    """
    Returns a reflection matrix that reflects a shape on a line defined by (0, 0) and the given (lx, ly) coordinates.

    Parameters
    ----------
    lx : int | float
        The x-coordinate of a point on the reflection line.
    ly : int | float
        The y-coordinate of a point on the reflection line.

    Returns
    -------
    np.ndarray
        A 2x2 numpy array representing the reflection matrix.
    """
    lx2 = lx ** 2
    ly2 = ly ** 2
    l2 = lx2 + ly2
    two_lxly = 2 * lx * ly
    square_diff_xy = lx2 - ly2

    ref_mat = np.array([[square_diff_xy, two_lxly], [two_lxly, - square_diff_xy]]) / l2
    return ref_mat


def translate(x: multi_num_ext, y: multi_num_ext, xt: multi_num_ext, yt: multi_num_ext):
    """
    Translate points in Cartesian coordinates by specified amounts in xt and yt.

    Parameters
    ----------
    x : array_like or Quantity
        The x-coordinates of the points to be translated.
    y : array_like or Quantity
        The y-coordinates of the points to be translated.
    xt : array_like or Quantity
        The amount to translate each point in x.
    yt : array_like or Quantity
        The amount to translate each point in y.

    Returns
    -------
    Tuple
        A tuple containing the new x-coordinates and y-coordinates of the translated points.
    """

    x = np.array(x) if not isinstance(x, Qty) else x
    y = np.array(y) if not isinstance(y, Qty) else y
    xt = np.array(xt) if not isinstance(xt, Qty) else xt
    yt = np.array(yt) if not isinstance(yt, Qty) else yt

    return x + xt, y + yt


def transform(x: multi_num_ext, y: multi_num_ext, trans_matrix: np.ndarray):
    """
    Transform points in Cartesian coordinates using a transformation matrix.

    Parameters
    ----------
    x : array_like or Quantity
        The x-coordinates of the points to be transformed.
    y : array_like or Quantity
        The y-coordinates of the points to be transformed.
    trans_matrix : array_like
        The 2x2 transformation matrix to apply to the points.

    Returns
    -------
    Tuple
        A tuple containing the new x-coordinates and y-coordinates of the transformed points.
    """

    if not isinstance(x, Iterable):
        x = [x]
    if not isinstance(y, Iterable):
        y = [y]

    units = None
    if isinstance(x, Qty):
        units = x.units
        x = x.magnitude

    if isinstance(y, Qty):
        units = y.units
        y = y.magnitude

    coords = np.array([x, y])

    xr, yr = np.squeeze(trans_matrix @ coords)
    return (xr, yr) if units is None else (xr * units, yr * units)


def rotate(x: multi_num_ext, y: multi_num_ext, angle: num_ext, counterclockwise=True):
    """
    Rotate points in Cartesian coordinates around the origin (0, 0).

    Parameters
    ----------
    x : array_like | Quantity
        The x-coordinates of the points to be rotated.
    y : array_like | Quantity
        The y-coordinates of the points to be rotated.
    angle : int | float | Quantity
        The angle to rotate the points, in degrees. If a Quantity is passed in, it must have dimensions of angle.
    counterclockwise : bool, optional
        Whether to rotate the points counterclockwise (default) or clockwise.

    Returns
    -------
    Tuple
        A tuple containing the new x-coordinates and y-coordinates of the rotated points.
    """

    rot_mat = rotation_matrix(angle, counterclockwise)
    x_trans, y_trans = transform(x, y, rot_mat)
    return x_trans, y_trans


def stretch(x: multi_num_ext, y: multi_num_ext, x_stretch: num, y_stretch: num):
    """
    Enlarges or Compresses a shape in the x and y axes with the shape's centroid at (0, 0).

    Parameters
    ----------
    x : array_like or Quantity
        The x-coordinates of the points to be stretched.
    y : array_like or Quantity
        The y-coordinates of the points to be stretched.
    x_stretch : float or Quantity
        The amount to stretch each point along the x-axis. If a Quantity is passed in,
         it must have dimensions of length.
    y_stretch : float or Quantity
        The amount to stretch each point along the y-axis. If a Quantity is passed in,
        it must have dimensions of length.

    Returns
    -------
    Tuple
        A tuple containing the new x-coordinates and y-coordinates of the stretched points.
    """

    stretch_mat = stretching_matrix(x_stretch, y_stretch)
    x_trans, y_trans = transform(x, y, stretch_mat)
    return x_trans, y_trans


def squeeze(x: multi_num_ext, y: multi_num_ext, squeeze_value: num):
    """
    Squeezes a shape by enlarging it in one direction and compressing it in the other, with the shape's
    centroid at (0, 0). Squeezing conserves the shape's area.

    Parameters:
    -----------
    x : multi_num_ext
        x-coordinates of the points to be squeezed.
    y : multi_num_ext
        y-coordinates of the points to be squeezed.
    squeeze_value : num
        Amount of squeezing.

    Returns:
    --------
    Tuple[multi_num_ext, multi_num_ext]
        The squeezed x and y coordinates.
    """

    squeeze_mat = squeezing_matrix(squeeze_value)
    x_trans, y_trans = transform(x, y, squeeze_mat)
    return x_trans, y_trans


def reflect(x: multi_num_ext, y: multi_num_ext, lx: multi_num_ext, ly: multi_num_ext):
    """
    Reflects the coordinates about a line passing through a point (lx, ly) and the origin (0, 0).

    Parameters:
    -----------
    x : multi_num_ext
        x-coordinates of the points to be reflected.
    y : multi_num_ext
        y-coordinates of the points to be reflected.
    lx : multi_num_ext
        x-coordinate of the point through which the line passes.
    ly : multi_num_ext
        y-coordinate of the point through which the line passes.

    Returns:
    --------
    Tuple[multi_num_ext, multi_num_ext]
        The reflected x and y coordinates.
    """
    ref_mat = reflection_matrix(lx, ly)
    x_trans, y_trans = transform(x, y, ref_mat)
    return x_trans, y_trans


def reflect_x(x: multi_num_ext, y: multi_num_ext):
    """
    Reflects the coordinates about the x-axis.

    Parameters:
    -----------
    x : multi_num_ext
        x-coordinates of the points to be reflected.
    y : multi_num_ext
        y-coordinates of the points to be reflected.

    Returns:
    --------
    Tuple[multi_num_ext, multi_num_ext]
        The reflected x and y coordinates.
    """

    ref_mat = reflection_matrix(0, 1)
    x_trans, y_trans = transform(x, y, ref_mat)
    return x_trans, y_trans


def reflect_y(x: multi_num_ext, y: multi_num_ext):
    """
    Reflects the coordinates about the y-axis.

    Parameters:
    -----------
    x : multi_num_ext
        x-coordinates of the points to be reflected.
    y : multi_num_ext
        y-coordinates of the points to be reflected.

    Returns:
    --------
    Tuple[multi_num_ext, multi_num_ext]
        The reflected x and y coordinates.
    """

    ref_mat = reflection_matrix(1, 0)
    x_trans, y_trans = transform(x, y, ref_mat)
    return x_trans, y_trans

