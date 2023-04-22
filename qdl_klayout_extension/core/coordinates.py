"""
This module contains classes that help with the creation and manipulation of 2D coordinates.
"""

import warnings
from typing import List, Iterable, Sequence, Tuple

import numpy as np
import pya
from multipledispatch import dispatch
from pint import Quantity as Qty
from qdl_klayout_extension.constants import num_ext, multi_num_ext, num
from qdl_klayout_extension.errors import XYLengthError, OperationTypeError
from qdl_klayout_extension.transformations import translate, rotate, transform, stretch, squeeze, reflect, reflect_x, \
    reflect_y
from qdl_klayout_extension.utils import v_in_uu, qty_in_uu, ulu2dbu, find_centroid, sort_rotationally, \
    line_intersection, line_segment_length, point_line_segment_distance


def get_coords(coords: Tuple | "Coordinates" | Sequence | "CoordinatesList") -> Tuple[Qty | None, Qty | None]:
    """
    Extracts x and y coordinates from a given set of coordinates or coordinates list.

    Parameters
    ----------
    coords : Tuple | Coordinates | Sequence | CoordinatesList
        The coordinates to extract from. It can be a tuple of x and y values,
        a single Coordinates object, a sequence of Coordinates objects,
        or a CoordinatesList object.

    Returns
    -------
    x : pint.Quantity | None
        A single value or an array with the x values of the given coordinates. If the input was
        not recognized, or if there were no valid x values in the input,
        returns None.
    y : numpy.ndarray or None
        A single value or an array with the y values of the given coordinates. If the input was
        not recognized, or if there were no valid y values in the input,
        returns None.
    """
    if hasattr(coords, "in_uu"):
        return coords.x_uu, coords.y_uu
    elif isinstance(coords, Tuple):
        return coords
    elif isinstance(coords, Sequence):
        if hasattr(coords[0], "in_uu"):
            list_of_tuples = [(coord.x_uu, coord.y_uu) for coord in coords]
            x, y = list(map(tuple, zip(*list_of_tuples)))
            return x, y
        elif isinstance(coords[0], Tuple):
            x, y = list(map(tuple, zip(*coords)))
            return x, y

    return None, None


class Coordinates:
    """
    A class representing coordinates in user-defined and database units.
    The class provides getter methods for accessing the x and y coordinates in both user-defined and database units,
    a tuple of the coordinates in each unit as well as a KLayout Point and Vector object representation.
    It also provides methods for getting a translated, transformed, or rotated version of the point.
    Lastly, it implements basic operators such as "+", "-" and "==" between coordinates, as well as "*" , "/", "//" and
    "%" between a coordinate and a numeric.
    """

    def __init__(self, x: num_ext, y: num_ext):
        """
        Parameters
        ----------
        x : float | int | pint.Quantity
            The x-coordinate value in user-defined units.
        y : float | int | pint.Quantity
            The y-coordinate value in user-defined units.

        """
        self._x_uu: Qty = qty_in_uu(x) if isinstance(x, Qty) else v_in_uu(x)
        self._y_uu: Qty = qty_in_uu(y) if isinstance(y, Qty) else v_in_uu(y)
        self._klayout_point = None if np.isnan(self.x_uu) or np.isnan(self.y_uu) else pya.Point(self.x_dbu, self.y_dbu)
        self._klayout_vector = None if np.isnan(self.x_uu) or np.isnan(self.y_uu) \
            else pya.Vector(self.x_dbu, self.y_dbu)

    @property
    def x_uu(self) -> Qty:
        """ pint.Quantity: The x-coordinate value in user-defined units as a Quantity object."""
        return self._x_uu

    @property
    def y_uu(self) -> Qty:
        """ pint.Quantity: The y-coordinate value in user-defined units as a Quantity object."""
        return self._y_uu

    @property
    def x_dbu(self) -> int:
        """ pint.Quantity: The x-coordinate value in database units."""
        return ulu2dbu(self._x_uu)

    @property
    def y_dbu(self) -> int:
        """ pint.Quantity: The y-coordinate value in database units."""
        return ulu2dbu(self._y_uu)

    @property
    def in_uu(self) -> Tuple[Qty, Qty]:
        """ Coordinates of the object in user units as a tuple. """
        return self._x_uu, self._y_uu

    @property
    def in_dbu(self) -> Tuple[int, int]:
        """ Coordinates of the object in database units as a tuple. """
        return self.x_dbu, self.y_dbu

    @property
    def klayout_point(self) -> pya.Point:
        """ KLayout Point object representing the coordinates of the object. """
        return self._klayout_point

    @property
    def klayout_vector(self) -> pya.Vector:
        """ KLayout Vector object representing the coordinates of the object. """
        return self._klayout_vector

    @dispatch(object)
    def get_translated(self, coords: Tuple | "Coordinates"):
        """
        Returns a new Coordinates object that is translated by the distance (x, y) passed in as a coordinate argument.

        Parameters
        ----------
        coords : tuple or Coordinates
            Tuple of x and y distance or a Coordinates object.

        Returns
        -------
        Coordinates
            A new Coordinates object that is translated by the given distance in (x, y).
        """
        x_coords, y_coords = get_coords(coords)
        x, y = translate(self.x_uu, self.y_uu, x_coords, y_coords)
        return Coordinates(x, y)

    @dispatch(object, object)
    def get_translated(self, x_coords: num_ext, y_coords: num_ext):
        """
        Returns a new Coordinates object that is translated by the distance in x, y passed in as arguments.

        Parameters
        ----------
        x_coords : int | float | pint.Quantity
            The x distance to translate the Coordinates object.
        y_coords : int | float | pint.Quantity
            The y distance to translate the Coordinates object.

        Returns
        -------
        Coordinates
            A new Coordinates object that is translated by the given distance in (x, y).
        """
        x, y = translate(self.x_uu, self.y_uu, x_coords, y_coords)
        return Coordinates(x, y)

    def get_transformed(self, trans_matrix: np.ndarray):
        """
        Returns a new Coordinates object that is transformed by the given transformation matrix.

        Parameters
        ----------
        trans_matrix : numpy.ndarray
            A 2x2 numpy array representing the transformation matrix.

        Returns
        -------
        Coordinates
            A new Coordinates object that is transformed by the given transformation matrix.
        """
        x, y = transform(self.x_uu, self.y_uu, trans_matrix)
        return Coordinates(x, y)

    def get_rotated(self, angle: num_ext, counterclockwise=True):
        """
        Returns a new Coordinates object that is rotated by the given angle.

        Parameters
        ----------
        angle : int | float | pint.Quantity
            The angle (in degrees) by which the Coordinates object is to be rotated.
        counterclockwise : bool, optional
            Whether the rotation should be counterclockwise (default) or clockwise.

        Returns
        -------
        Coordinates
            A new Coordinates object that is rotated by the given angle.
        """
        x, y = rotate(self.x_uu, self.y_uu, angle, counterclockwise)
        return Coordinates(x, y)

    def get_stretched(self, x_stretch: num, y_stretch: num):
        """
        Returns a new Coordinates object that is stretched along the x and/or y axes by the given factors.

        Parameters
        ----------
        x_stretch : int | float
            The factor by which the Coordinates object is to be stretched along the x-axis.
        y_stretch : int | float
            The factor by which the Coordinates object is to be stretched along the y-axis.

        Returns
        -------
        Coordinates
            A new Coordinates object that is stretched along the x and/or y axes by the given factors.
        """
        x, y = stretch(self.x_uu, self.y_uu, x_stretch, y_stretch)
        return Coordinates(x, y)

    def get_squeezed(self, squeeze_value: num):
        """
        Returns a new Coordinates object that is squeezed  by the given factor. Squeezing conserves the shape area.

        Parameters
        ----------
        squeeze_value : int | float
            The factor by which the Coordinates object is to be squeezed or expanded.

        Returns
        -------
        Coordinates
            A new Coordinates object that is squeezed or expanded by the given factor.
        """
        x, y = squeeze(self.x_uu, self.y_uu, squeeze_value)
        return Coordinates(x, y)

    @dispatch(object)
    def get_reflected(self, coords: Tuple | "Coordinates"):
        """
        Returns a new Coordinates object that is reflected about the line passing through (o, 0)
        and the given coordinates.

        Parameters
        ----------
        coords : tuple | Coordinates
            Tuple of x and y coordinates or a Coordinates object representing the line about which the Coordinates
            object is to be reflected.

        Returns
        -------
        Coordinates
            A new Coordinates object that is reflected about the line passing through (0, 0) and the given coordinates.
        """
        x_coords, y_coords = get_coords(coords)
        x, y = reflect(self.x_uu, self.y_uu, x_coords, y_coords)
        return Coordinates(x, y)

    @dispatch(object, object)
    def get_reflected(self, x_coords: num_ext, y_coords: num_ext):
        """
        Returns a new Coordinates object that is reflected about the line passing through (0, 0)
        and the given coordinates (x_coords, y_coords).

        Parameters
        ----------
        x_coords : int | float | pint.Quantity
            The x coordinates of the defining the reflection line (0, 0), (x_coords, y_coords).
        y_coords : int | float | pint.Quantity
            The y coordinates of the defining the reflection line (0, 0), (x_coords, y_coords).

        Returns
        -------
        Coordinates
            A new Coordinates object that is reflected about the line passing through (0, 0) and the given coordinates.
        """
        x, y = reflect(self.x_uu, self.y_uu, x_coords, y_coords)
        return Coordinates(x, y)

    def get_reflected_x(self):
        """
        Returns a new Coordinates object that is reflected about the x-axis.

        Returns
        -------
        Coordinates
            A new Coordinates object that is reflected about the x-axis.
        """
        x, y = reflect_x(self.x_uu, self.y_uu)
        return Coordinates(x, y)

    def get_reflected_y(self):
        """
        Returns a new Coordinates object that is reflected about the y-axis.

        Returns
        -------
        Coordinates
            A new Coordinates object that is reflected about the y-axis.
        """
        x, y = reflect_y(self.x_uu, self.y_uu)
        return Coordinates(x, y)

    def __repr__(self):
        return repr(self.in_uu)

    def __add__(self, other: "Coordinates" | Tuple[num_ext]):
        x, y = get_coords(other)
        if x is not None and y is not None:
            return Coordinates(self.x_uu + x, self.y_uu + y)
        else:
            raise OperationTypeError(self, other, '+')

    def __sub__(self, other: "Coordinates" | Tuple[num_ext]):
        x, y = get_coords(other)
        if x is not None and y is not None:
            return Coordinates(self.x_uu - x, self.y_uu - y)
        else:
            raise OperationTypeError(self, other, '-')

    def __pos__(self):
        return Coordinates(+self.x_uu, +self.y_uu)

    def __neg__(self):
        return Coordinates(-self.x_uu, -self.y_uu)

    def __mul__(self, other):
        if isinstance(other, num):
            return Coordinates(self.x_uu * other, self.y_uu * other)
        else:
            raise OperationTypeError(self, other, '*')

    def __truediv__(self, other):
        if isinstance(other, num):
            return Coordinates(self.x_uu / other, self.y_uu / other)
        else:
            raise OperationTypeError(self, other, '/')

    def __floordiv__(self, other):
        if isinstance(other, num):
            return Coordinates(self.x_uu // other, self.y_uu // other)
        else:
            raise OperationTypeError(self, other, '//')

    def __mod__(self, other):
        if isinstance(other, num):
            return Coordinates(self.x_uu % other, self.y_uu % other)
        else:
            raise OperationTypeError(self, other, '%')

    def __eq__(self, other):
        x, y = get_coords(other)
        if x is not None and y is not None:
            x_condition = (np.isnan(self.x_uu) and np.isnan(x)) or self.x_uu == x
            y_condition = (np.isnan(self.y_uu) and np.isnan(y)) or self.y_uu == y
            return x_condition and y_condition
        else:
            raise OperationTypeError(self, other, '==')


class Line:
    """
    This class represents a line in 2D Cartesian space. It takes two points as input, which define the line.
    It calculates the slope and y-intercept of the line using these points, as well as properties to access
    the two points, slope and y-intercept of the line. Lastly, it also provides a method to calculate the
    intersection point between two lines.
    """
    def __init__(self, point_1: Tuple | Coordinates, point_2: Tuple | Coordinates):
        """
        Parameters
        ----------
        point_1 : Tuple | Coordinates
            The first point that defines the line.
        point_2 : Tuple | Coordinates
            The second point that defines the line.
        """

        if isinstance(point_1, Tuple):
            point_1 = Coordinates(*point_1)
        if isinstance(point_2, Tuple):
            point_2 = Coordinates(*point_2)

        self._point_1 = point_1
        self._point_2 = point_2

        self.__post_init__()

    def __post_init__(self):
        x1, y1 = self.point_1.x_uu, self.point_1.y_uu
        x2, y2 = self.point_2.x_uu, self.point_2.y_uu
        self._slope = (y2 - y1) / (x2 - x1)
        self._intercept = y1 - self._slope * x1
        self._length = line_segment_length((x1, y1), (x2, y2))

    @property
    def point_1(self):
        """ The first point that defines the line. """
        return self._point_1

    @property
    def point_2(self):
        """ The second point that defines the line. """
        return self._point_2

    @property
    def slope(self):
        """ The slope of the line. """
        return self._slope

    @property
    def intercept(self):
        """ The intercept of the line. """
        return self._intercept

    @property
    def length(self):
        """ The length of the line (distance between the two points). """
        return self._length

    def intersection(self, other: "Line"):
        """
        Returns the intersection point of the line with another line.

        Parameters
        ----------
        other : Line
            The other line to intersect with.

        Returns
        -------
        Coordinates | None
            The intersection point of the two lines as a Coordinates object, or `None` if they don't intersect.
        """

        x, y = line_intersection(self.point_1.in_uu, self.point_2.in_uu, other.point_1.in_uu, other.point_2.in_uu)

        if x is not None:
            return Coordinates(x, y)

        return None

    def distance_from_point(self, point: Coordinates | Tuple[num_ext | num_ext]):
        point_line_segment_distance()

    def __repr__(self):
        return f"Line({self.point_1}, {self.point_2})"


class CoordinatesList:
    """
    The CoordinatesList class is used for managing lists of Coordinates, and it can be initialized in two ways:

    - With a list of Coordinates
    - With X and Y coordinates

    The class provides getter methods for accessing the x and y coordinates in both user-defined and database units,
    a list of tuples of the coordinates in each unit as well as lists of KLayout Point and Vector objects.
    It also provides methods for getting a translated, transformed, or rotated version of the point.
    Additionally, it includes methods for computing various geometric properties of a polygon represented by the list
    of coordinates, such as the centroid, coordinates between two adjacent points, lines and length of lines between
    adjacent points, angles between edges and the y-axis, and angles between adjacent edges.

    Parameters
    ----------
    coords : Sequence[Coordinates]
        List of Coordinates to be managed.
    ref_point : Coordinates, optional
        Reference point for the CoordinatesList. It is used for rotational sorting purposes and defaults to the
        coordinates' centroid point.
    x : Sequence[num_ext] | Qty
        X-coordinates of the CoordinatesList, by default None.
    y : Sequence[num_ext] | Qty
        Y-coordinates of the CoordinatesList, by default None.
    """

    @dispatch(object, ref_point=Coordinates)
    def __init__(self, coords: Sequence[Coordinates], ref_point: Coordinates = Coordinates(np.nan, np.nan)):
        """
        Initialize CoordinatesList with a list of Coordinates.

        Parameters
        ----------
        coords : Sequence[Coordinates]
            List of Coordinates to be managed.
        ref_point : Coordinates, optional
            Reference point for the CoordinatesList. It is used for rotational sorting purposes and defaults to the
            coordinates' centroid point.
        """

        self._coords: List[Coordinates] = list(coords)

        self._x_uu = Qty.from_list([coord.x_uu for coord in self._coords])
        self._y_uu = Qty.from_list([coord.y_uu for coord in self._coords])

        self._ref_point = ref_point
        self.__post_init__()

    @dispatch(object, object, ref_point=Coordinates)
    def __init__(self, x: Sequence[num_ext] | Qty, y: Sequence[num_ext] | Qty,
                 ref_point: Coordinates = Coordinates(np.nan, np.nan)):
        """
        Initialize CoordinatesList with X and Y coordinates.

        Parameters
        ----------
        x : Sequence[num_ext] | Qty
            X-coordinates of the CoordinatesList.
        y : Sequence[num_ext] | Qty
            Y-coordinates of the CoordinatesList.
        ref_point : Coordinates, optional
            Reference point for the CoordinatesList. It is used for rotational sorting purposes and defaults to the
            coordinates' centroid point.
        """

        x, y = self._are_inputs_valid(x, y)

        self._x_uu = qty_in_uu(x) if isinstance(x, Qty) else v_in_uu(x)
        self._y_uu = qty_in_uu(y) if isinstance(y, Qty) else v_in_uu(y)

        self._coords: List[Coordinates] = [Coordinates(self._x_uu[i], self._y_uu[i]) for i in range(len(x))]

        self._ref_point = ref_point
        self.__post_init__()

    def __post_init__(self):
        self._ref_point = self.get_centroid() if self._ref_point == Coordinates(np.nan, np.nan) else self._ref_point
        self._is_sorted = self._check_sorted()
        self._klayout_points = [pya.Point(x, y) for x, y in self.in_dbu]
        self._klayout_vectors = [pya.Vector(x, y) for x, y in self.in_dbu]

    @staticmethod
    def _are_inputs_valid(x, y) -> Tuple[Sequence, Sequence]:
        """
        Check if input x and y are valid.

        Parameters
        ----------
        x : Sequence[num_ext] | Qty
            X-coordinates of the CoordinatesList.
        y : Sequence[num_ext] | Qty
            Y-coordinates of the CoordinatesList.

        Returns
        -------
        Tuple[Sequence, Sequence]
            Validated x and y coordinates.

        Raises
        ------
        XYLengthError
            If the lengths between the x and y arrays do not match.
        """
        if not isinstance(x, Iterable):
            x, y = [x], [y]

        if len(x) != len(y):
            raise XYLengthError(x, y)

        return x, y

    def _check_sorted(self):
        """
        Check if the CoordinatesList is sorted.

        Returns
        -------
        int
            -1 if sorted clockwise, 1 if sorted counterclockwise, and 0 if not sorted.
        """

        x_sorted, y_sorted = sort_rotationally(self.x_uu, self.y_uu, xref=self.ref_point.x_uu, yref=self.ref_point.y_uu)

        if all(self._x_uu == x_sorted):
            return -1

        x_sorted, y_sorted = sort_rotationally(self.x_uu, self.y_uu, counterclockwise=True, xref=self.ref_point.x_uu,
                                               yref=self.ref_point.y_uu)
        if all(self._x_uu == x_sorted):
            return 1

        return 0

    @property
    def coords_list(self) -> List[Coordinates]:
        """ List of Coordinates managed by the CoordinatesList. """
        return self._coords

    @property
    def x_uu(self) -> Qty:
        """ The x coordinates of the CoordinatesList in user units. """
        return self._x_uu

    @property
    def y_uu(self) -> Qty:
        """ The y coordinates of the CoordinatesList in user units. """
        return self._y_uu

    @property
    def x_dbu(self) -> List[int]:
        """ The x coordinates of the CoordinatesList in database units. """
        return ulu2dbu(self._x_uu)

    @property
    def y_dbu(self) -> List[int]:
        """ The y coordinates of the CoordinatesList in database units. """
        return ulu2dbu(self._y_uu)

    @property
    def in_uu(self) -> List[Tuple[Qty, Qty]]:
        """ List of Coordinates of the CoordinatesList in user units as a tuple with two pint.Quantity elements. """
        return [coord.in_uu for coord in self._coords]

    @property
    def in_dbu(self) -> List[Tuple[int, int]]:
        """ List of Coordinates of the CoordinatesList in database units as a tuple with two int elements. """
        return [coord.in_dbu for coord in self._coords]

    @property
    def ref_point(self):
        """ Reference point for the CoordinatesList. The reference point is use for rotational sorting purposes. """
        return self._ref_point

    @ref_point.setter
    def ref_point(self, ref_point: Coordinates):
        self._ref_point = ref_point

    @property
    def klayout_points(self) -> List[pya.Point]:
        """ A list of KLayout Point objects representing the CoordinatesList. """
        return self._klayout_points

    @property
    def klayout_vectors(self) -> List[pya.Vector]:
        """ A list of KLayout Vector objects representing the CoordinatesList. """
        return self._klayout_vectors

    @property
    def is_sorted(self):
        """ Returns 0 if not sorted, -1 if clockwise, 1 if counterclockwise. """
        return self._is_sorted

    def sort_clockwise(self, ref_point: Coordinates | None = None):
        """
        Sorts the coordinates of the polygon in clockwise order with respect to the reference point.

        Parameters
        ----------
        ref_point : Coordinates or None, optional
            The reference point to sort the polygon around. If None, uses the current reference point.
        """

        xref = self.ref_point.x_uu if ref_point is None else ref_point.x_uu
        yref = self.ref_point.y_uu if ref_point is None else ref_point.y_uu
        x_sorted, y_sorted = sort_rotationally(self.x_uu, self.y_uu, xref=xref, yref=yref)
        self.__init__(x_sorted, y_sorted, ref_point=Coordinates(xref, yref))

    def sort_counterclockwise(self, ref_point: Coordinates | None = None):
        """
        Sorts the coordinates of the polygon in counterclockwise order with respect to the reference point.

        Parameters
        ----------
        ref_point : Coordinates or None, optional
            The reference point to sort the polygon around. If None, uses the current reference point.
        """
        xref = self.ref_point.x_uu if ref_point is None else ref_point.x_uu
        yref = self.ref_point.y_uu if ref_point is None else ref_point.y_uu

        x_sorted, y_sorted = sort_rotationally(self.x_uu, self.y_uu, counterclockwise=True, xref=xref, yref=yref)
        self.__init__(x_sorted, y_sorted, ref_point=Coordinates(xref, yref))

    def get_edge_centers(self) -> "CoordinatesList":
        """
        Calculates the coordinates of the centers of the edges of the polygon.

        Returns
        -------
        CoordinatesList
            A list of the coordinates of the centers of the edges.
        """
        if not self._is_sorted:
            warnings.warn('Edge centers may not be correct since the coordinate list is not sorted.')

        x = list(self._x_uu) + [self._x_uu[0]]
        y = list(self._y_uu) + [self._y_uu[0]]

        x_edge_centers = Qty.from_list([(x[i] + x[i + 1]) / 2 for i in range(len(self))])
        y_edge_centers = Qty.from_list([(y[i] + y[i + 1]) / 2 for i in range(len(self))])

        return CoordinatesList(x_edge_centers, y_edge_centers)

    def get_centroid(self) -> Coordinates:
        """
        Calculates the centroid of the polygon.

        Returns
        -------
        Coordinates
            The centroid of the polygon.
        """
        xc, yc = find_centroid(self._x_uu, self._y_uu)
        return Coordinates(xc, yc)

    def get_edge_lines(self) -> List[Line]:
        """
        Returns a list of the edge lines of the polygon.

        Returns
        -------
        List[Line]
            A list of the edge lines of the polygon.
        """
        if not self._is_sorted:
            warnings.warn('Edge lines may not be correct since the coordinate list is not sorted.')

        coords = list(self._coords) + [self._coords[0]]

        lines = [Line(coords[i], coords[i + 1]) for i in range(len(self))]
        return lines

    def get_edge_lengths(self) -> Qty:
        """
        Calculates the lengths of the edges of the polygon.

        Returns
        -------
        pint.Quantity[np.ndarray]
            An array containing the lengths of the edges.
        """
        lines = self.get_edge_lines()

        lengths = Qty.from_list([line.length for line in lines])
        return lengths

    def get_edge_angles(self) -> Qty:
        """
        Calculates the angles between the edges of the polygon and the y-axis.

        Returns
        -------
        pint.Quantity[np.ndarray]
            An array containing the angles between the edges and the y-axis.
        """
        if not self._is_sorted:
            warnings.warn('Edge angles may not be correct since the coordinate list is not sorted.')

        x = list(self._x_uu) + [self._x_uu[0]]
        y = list(self._y_uu) + [self._y_uu[0]]

        angles = Qty.from_list([np.arctan2(y[i + 1] - y[i], x[i + 1] - x[i]) % Qty(360, 'deg')
                                for i in range(len(self))])
        return angles

    def get_point_angles(self) -> Qty:
        """
        Calculates the angles between the adjacent edges of the polygon.

        Returns
        -------
        pint.Quantity[np.ndarray]
            An array containing the angles between the adjacent edges.
        """
        if not self._is_sorted:
            warnings.warn('Point angles may not be correct since the coordinate list is not sorted.')

        edge_angles = self.get_edge_angles()
        edge_angles = [edge_angles[-1]] + list(edge_angles)

        angles = Qty.from_list([(edge_angles[i + 1] - edge_angles[i]) % Qty(360, 'deg') for i in range(len(self))])
        return angles

    @dispatch(object)
    def get_translated(self, coords: Tuple | Coordinates | Sequence | "CoordinatesList"):
        """
        Translate the CoordinatesList by given values in Coordinates.

        Parameters
        ----------
        coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The distance to translate the CoordinatesList.

        Returns
        -------
        CoordinatesList
            The translated coordinates.
        """

        x_coords, y_coords = get_coords(coords)
        x, y = translate(self.x_uu, self.y_uu, x_coords, y_coords)
        ref_point = self.ref_point.get_translated(coords)
        return CoordinatesList(x, y, ref_point=ref_point)

    @dispatch(object, object)
    def get_translated(self, x_coords: multi_num_ext, y_coords: multi_num_ext):
        """
        Translate the CoordinatesList by given x and y values.

        Parameters
        ----------
        x_coords : int | float | pint.Quantity | Sequence[int | float | pint.Quantity]
            The x distance to translate the CoordinatesList in user units.
            If a list, it must be the same size as the CoordinateList to perform element-wise translation
        y_coords : int | float | pint.Quantity | Sequence[int | float | pint.Quantity]
            The y distance to translate the CoordinatesList in user units.
            If a list, it must be the same size as the CoordinateList to perform element-wise translation

        Returns
        -------
        CoordinatesList
            The translated coordinates.
        """

        x, y = translate(self.x_uu, self.y_uu, x_coords, y_coords)
        ref_point = self.ref_point.get_translated(x_coords, y_coords)
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_transformed(self, trans_matrix: np.ndarray):
        """
        Transform the coordinates using the given transformation matrix.

        Parameters
        ----------
        trans_matrix : np.ndarray
            The transformation matrix.

        Returns
        -------
        CoordinatesList
            The transformed coordinates.
        """

        x, y = transform(self.x_uu, self.y_uu, trans_matrix)
        ref_point = self.ref_point.get_transformed(trans_matrix)
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_rotated(self, angle: num_ext, counterclockwise=True):
        """

        Rotate the coordinates by a given angle around the origin.

        Parameters
        ----------
        angle : int | float | pint.Quantity
            The angle to rotate by in degrees.
        counterclockwise : bool, optional
            If True (default), rotate counterclockwise. Otherwise, rotate clockwise.

        Returns
        -------
        CoordinatesList
            The rotated coordinates.
        """
        x, y = rotate(self.x_uu, self.y_uu, angle, counterclockwise)
        ref_point = self.ref_point.get_rotated(angle, counterclockwise)
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_stretched(self, x_stretch: num, y_stretch: num):
        """
        Stretch (enlarge or compress) the coordinates by given x and y factors.

        Parameters
        ----------
        x_stretch : int | float
            The factor to stretch the x-coordinates by.
        y_stretch : int | float
            The factor to stretch the y-coordinates by.

        Returns
        -------
        CoordinatesList
            The stretched coordinates.
        """

        x, y = stretch(self.x_uu, self.y_uu, x_stretch, y_stretch)
        ref_point = self.ref_point.get_stretched(x_stretch, y_stretch)
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_squeezed(self, squeeze_value: num):
        """
        Squeeze the coordinates by given factor. Squeezing conserves the shape's area.

        Parameters
        ----------
        squeeze_value : int | float
            The factor to squeeze the coordinates by.

        Returns
        -------
        CoordinatesList
            The squeezed coordinates.
        """

        x, y = squeeze(self.x_uu, self.y_uu, squeeze_value)
        ref_point = self.ref_point.get_squeezed(squeeze_value)
        return CoordinatesList(x, y, ref_point=ref_point)

    @dispatch(object)
    def get_reflected(self, coords: Tuple | Coordinates | Sequence | "CoordinatesList"):
        """
        Reflect the coordinates across a line defined by the given coordinates (lx, ly) and the origin (0, 0).

        Parameters
        ----------
        coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The coordinates that define the reflection line.

        Returns
        -------
        CoordinatesList
            The reflected coordinates.
        """

        x_coords, y_coords = get_coords(coords)
        x, y = reflect(self.x_uu, self.y_uu, x_coords, y_coords)
        ref_point = self.ref_point.get_reflected(coords)
        return CoordinatesList(x, y, ref_point=ref_point)

    @dispatch(object, object)
    def get_reflected(self, x_coords: multi_num_ext, y_coords: multi_num_ext):
        """
        Reflect the coordinates across a line defined by the given coordinates (lx, ly) and the origin (0, 0).

        Parameters
        ----------
        x_coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The x coordinates in user units that define the reflection line.
        y_coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The y coordinates in user units that define the reflection line.

        Returns
        -------
        CoordinatesList
            The reflected coordinates.
        """
        x, y = reflect(self.x_uu, self.y_uu, x_coords, y_coords)
        ref_point = self.ref_point.get_reflected(x_coords, y_coords)
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_reflected_x(self):
        """
        Reflect the coordinates across the x-axis.

        Returns
        -------
        CoordinatesList
            The reflected coordinates.
        """

        x, y = reflect_x(self.x_uu, self.y_uu)
        ref_point = self.ref_point.get_reflected_x()
        return CoordinatesList(x, y, ref_point=ref_point)

    def get_reflected_y(self):
        """
        Reflect the coordinates across the y-axis.

        Returns
        -------
        CoordinatesList
            The reflected coordinates.
        """
        x, y = reflect_y(self.x_uu, self.y_uu)
        ref_point = self.ref_point.get_reflected_y()
        return CoordinatesList(x, y, ref_point=ref_point)

    def __iter__(self):
        return iter(self._coords)

    def __len__(self):
        return len(self._coords)

    def __repr__(self):
        return repr(self._coords)

    def __getitem__(self, item):
        return self._coords[item]

