"""
This module contains classes that help create various KLayout.Polygon shapes easily.
"""


from typing import Tuple, Sequence, List

import numpy as np
import pya
from pint import Quantity as Qty
from matplotlib import path as mpl_path

from qdl_klayout_extension.constants import num_ext
from qdl_klayout_extension.core.coordinates import CoordinatesList, Coordinates, Line
from qdl_klayout_extension.core.geometries import rectangle_coords
from qdl_klayout_extension.errors import ShapePointCountError, ShapePointAngleError, ShapeEdgeLengthError
from qdl_klayout_extension.utils import qty_in_uu, v_in_uu

_deviation_precision = 0.01


def _deviation_is_above_precision_limit(value1: num_ext, value2: num_ext):
    """
    Check if the deviation between two values is above the precision limit given by _deviation_precision.

    Parameters
    ----------
    value1 : int | float | pint.Quantity
        First value.
    value2 : int | float | pint.Quantity
        Second value.

    Returns
    -------
    bool
        True if the deviation is above the precision limit, False otherwise.
    """

    allowed_deviation = (np.abs(value1) + np.abs(value2)) / 2 * _deviation_precision
    diff = np.abs(value1 - value2)
    return diff > allowed_deviation


class SimplePolygon:
    """
    The SimplePolygon class is a class for representing and manipulating simple polygons in two-dimensional space.
    It defines a polygon as a list of coordinates that represent the vertices of the polygon in order.
    It provides various definitions of various geometric properties of the polygon, such as edge centers, edge angles to
    the y-axis, edge lines and adjacent edge angles.
    Additionally, it provides methods for checking whether a point or set of points are inside or outside the polygon.

    The SimplePolygon class serves as a foundation for more specialized polygon classes, such as the SimpleRectangle
    and SimpleCross classes, which inherit from it and provide additional functionality specific to those shapes.
    """
    _size: int | None = None

    def __init__(self, coords: CoordinatesList):
        """
        Parameters:
        -----------
        coords: CoordinatesList
            A list of coordinates representing the vertices of the rectangle.
        """

        if not coords.is_sorted:
            coords.sort_clockwise()

        self._coords: CoordinatesList = coords
        self.__post_init__()
        self._check_all()

    def __post_init__(self):
        self._klayout_simple_polygon: pya.SimplePolygon = pya.SimplePolygon(self._coords.klayout_points)
        self._centroid: Coordinates = self._coords.get_centroid()
        self._ref_point: Coordinates = self._coords.ref_point

        self._edge_centers: CoordinatesList = self._coords.get_edge_centers()
        self._edge_lines: List[Line] = self._coords.get_edge_lines()
        self._edge_lengths = self._coords.get_edge_lengths()

        self._edge_angles = self._coords.get_edge_angles()
        self._point_angles = self._coords.get_point_angles()

        if self._size is None:
            self._size = len(self._coords)

    def _check_all(self):
        """ Helper method to call all check methods. """
        self._check_coord_size()
        self._check_edge_lengths()
        self._check_point_angles()

    def _check_coord_size(self):
        """ Helper method to check if the number of coordinates is valid. """
        pass

    def _check_edge_lengths(self):
        """ Helper method to check if the edge lengths are valid. """
        pass

    def _check_point_angles(self):
        """ Helper method to check if the point angles are valid. """
        pass

    @property
    def coords(self) -> CoordinatesList:
        """ List of coordinate tuples of the polygon vertices. """
        return self._coords

    @property
    def klayout_simple_polygon(self) -> pya.SimplePolygon:
        """ KLayout SimplePolygon object representing the polygon vertices. """
        return self._klayout_simple_polygon

    @property
    def centroid(self) -> Coordinates:
        """ Coordinates of the centroid of the polygon. """
        return self._centroid

    @property
    def ref_point(self) -> Coordinates:
        """ Coordinates of the reference point of the polygon. """
        return self._ref_point

    @property
    def edge_centers(self) -> CoordinatesList:
        """ List of coordinate tuples of the centers of each polygon edge."""
        return self._edge_centers

    @property
    def edge_lines(self):
        """ List of Line objects representing each polygon edge. """
        return self._edge_lines

    @property
    def edge_lengths(self):
        """ List of lengths of each polygon edge in user-defined units. """
        return self._edge_lengths

    @property
    def edge_angles(self):
        """ List of angles between each polygon edge and the y-axis in user-defined units. """
        return self._edge_angles

    @property
    def point_angles(self):
        """ List of angles between adjacent edges in the polygon in user-defined units. """
        return self._point_angles

    def contains_point(self, point: Coordinates):
        """
        Check if a point is contained within the polygon.

        Parameters
        ----------
        point : Coordinates
            The point to check if it's contained in the polygon.

        Returns
        -------
        bool
            True if the point is contained within the polygon, False otherwise.
        """
        polygon_coordinates = [(coord[0].m, coord[1].m) for coord in self.coords.in_uu]
        polygon_path = mpl_path.Path(polygon_coordinates)

        return polygon_path.contains_point((point.x_uu.m, point.y_uu.m))

    def contains_points(self, points: CoordinatesList | Sequence[Coordinates]):
        """
        Check if a list of points are contained within the polygon.

        Parameters
        ----------
        points : CoordinatesList | Sequence[Coordinates]
            The list of points to check if they are contained in the polygon.

        Returns
        -------
        List[bool]
            A list of booleans indicating if each point is contained within the polygon.
        """
        return [self.contains_point(point) for point in points]


def coords_for_simple_polygon_merge(polygon_1: SimplePolygon, polygon_2: SimplePolygon,
                                    ref_point: Coordinates | None = None) -> CoordinatesList:

    """
    Merge two SimplePolygon objects by returning a list of Coordinates representing the merged polygon.

    Parameters
    ----------
    polygon_1 : SimplePolygon
        The first SimplePolygon object to merge.
    polygon_2 : SimplePolygon
        The second SimplePolygon object to merge.
    ref_point : Coordinates, optional
        The reference point to use for rotation sorting merged polygon coordinates, by default is the centroid.

    Returns
    -------
    CoordinatesList
        A CoordinatesList object representing the merged polygon.
    """

    intersections = [l2.intersection(l1) for l1 in polygon_1.edge_lines for l2 in polygon_2.edge_lines
                     if l2.intersection(l1) is not None]

    p2_in_p1 = polygon_1.contains_points(polygon_2.coords)
    p1_in_p2 = polygon_2.contains_points(polygon_1.coords)

    polygon_1_coords = [coord for coord, contained in zip(polygon_1.coords.coords_list, p1_in_p2) if not contained]
    polygon_2_coords = [coord for coord, contained in zip(polygon_2.coords.coords_list, p2_in_p1) if not contained]

    coords_list = polygon_1_coords + polygon_2_coords + intersections

    coords = CoordinatesList([cc for cc in coords_list])
    if ref_point is not None:
        coords.ref_point = ref_point
        coords.sort_clockwise()

    return coords


class SimpleRectangle(SimplePolygon):
    """
    The SimpleRectangle class is a subclass of the SimplePolygon class that represents a rectangle shape.
    It inherits all the attributes and methods of the SimplePolygon class and adds specific methods to check if the
    provided coordinates define a valid rectangle shape.

    The class provides two alternative constructors for creating a SimpleRectangle object. The first constructor
    from_centroid() defines the rectangle location in 2D space by its centroid coordinates. The second
    constructor from_edge() defines the rectangle location in 2D space by one of the four edge's coordinates.
    """
    _size = 4

    def __init__(self, coords: CoordinatesList):
        """
        Parameters:
        -----------
        coords: CoordinatesList
            A list of coordinates representing the vertices of the rectangle.
        """
        super().__init__(coords)

    def _check_coord_size(self):
        """
        Checks if the number of coordinates is equal to 4.

        Raises:
        -------
        ShapePointCountError:
            If the number of coordinates is not equal to 4.
        """
        if not len(self._coords) == self._size:
            raise ShapePointCountError('simple rectangle', self._size, len(self._coords))

    def _check_edge_lengths(self):
        """
        Checks if the lengths of the opposite edges are equal.

        Raises:
        -------
        ShapeEdgeLengthError:
            If the opposite edges of the rectangle are not equal.
        """

        if _deviation_is_above_precision_limit(self._edge_lengths[0], self._edge_lengths[2]) or \
                _deviation_is_above_precision_limit(self._edge_lengths[1], self._edge_lengths[3]):
            raise ShapeEdgeLengthError('simple rectangle', self._edge_lengths)

    def _check_point_angles(self):
        """
        Checks if the angles of the vertices are 90 degrees.

        Raises:
        -------
        ShapePointAngleError:
            If the angles of the vertices are not 90 degrees.
        """

        if any([_deviation_is_above_precision_limit(pa, Qty(90, 'deg')) for pa in self._point_angles]):
            raise ShapePointAngleError('simple rectangle', self._point_angles.to('deg'))

    @classmethod
    def from_centroid(cls, width: num_ext, length: num_ext, angle: num_ext = 0,
                      centroid: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0)):
        """
        Constructs a SimpleRectangle instance given the width and length of the rectangle, its angle of rotation,
        and the coordinates of the rectangle's centroid.

        Parameters:
        -----------
        width: int | float | pint.Quantity
            The width of the rectangle.
        length: int | float | pint.Quantity
            The length of the rectangle.
        angle: int | float | pint.Quantity
            The angle in degrees to rotate the rectangle counterclockwise around its centroid.
        centroid: Coordinates | Tuple[num_ext, num_ext]
            The centroid of the rectangle.

        Returns:
        --------
        SimpleRectangle
            A SimpleRectangle object constructed from the provided parameters.
        """
        if isinstance(centroid, Tuple):
            centroid = Coordinates(*centroid)

        coords00 = rectangle_coords(width, length)
        coords_rotated00 = coords00.get_rotated(angle)
        coords = coords_rotated00.get_translated(centroid)

        return cls(coords)

    @classmethod
    def from_edge(cls, width: num_ext, length: num_ext, angle: num_ext = 0,
                  edge: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0), edge_type='bottom'):
        """
        Constructs a SimpleRectangle instance given the width and length of the rectangle, its angle of rotation,
        and the coordinates of a point on the edge of the rectangle.

        Parameters
        ----------
        width : int | float | pint.Quantity
            The width of the rectangle.
        length : int | float | pint.Quantity
            The length of the rectangle.
        angle : int | float | pint.Quantity, optional
            The angle in degrees by which the rectangle is rotated around its centroid. The default is 0.
        edge : Coordinates or Tuple[int | float | pint.Quantity, int | float | pint.Quantity], optional
            A point on an edge of the rectangle. The default is Coordinates(0, 0).
        edge_type : {'bottom', 'top', 'left', 'right'}, optional
            The type of the edge that the input edge point is located on. The default is 'bottom'.

        Returns
        -------
        SimpleRectangle
            A SimpleRectangle instance with the specified properties.

        Raises
        ------
        ValueError
            If `edge_type` is not one of {'bottom', 'top', 'left', 'right'}.
        """
        if isinstance(edge, Tuple):
            edge = Coordinates(*edge)

        acceptable_edge_types = ['bottom', 'top', 'left', 'right']
        if edge_type not in acceptable_edge_types:
            raise ValueError(f"Edge type cannot be {edge_type}. Choose between {acceptable_edge_types} instead.")

        edge_idxs = {'bottom': 0, 'right': 1, 'top': 2, 'left': 3}

        coords00 = rectangle_coords(width, length)
        coords_rotated00 = coords00.get_rotated(angle)
        edge_coord: Coordinates = coords_rotated00.get_edge_centers()[edge_idxs[edge_type]]
        coords = coords_rotated00.get_translated(edge - edge_coord)

        return cls(coords)


class SimpleCross(SimplePolygon):
    """
    The SimpleCross class is a subclass of the SimplePolygon class that defines a simple (potentially asymmetric) cross
    shape, which is made up of two rectangles.
    The from_rectangles() method constructs a SimpleCross instance from the parameters that define the two rectangles.
    The reference point of this shape is the interaction of the two overlapping rectangles.
    """

    def __init__(self, coords: CoordinatesList):
        """
        Parameters:
        -----------
        coords: CoordinatesList
            A list of coordinates representing the vertices of the rectangle.
        """
        super().__init__(coords)

    @classmethod
    def from_rectangles(cls, width: num_ext, length_1: num_ext, length_2: num_ext,
                        angle: num_ext = 0, ref_point: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0),
                        shift_1: num_ext = 0, shift_2: num_ext = 0):
        """
        Constructs a SimpleCross instance from the dimensions of the two rectangles that define it.

        Parameters
        ----------
        width : float or pint.Quantity
            The width of the cross for both rectangles
        length_1 : float or pint.Quantity
            The length of the first (vertical) rectangle.
        length_2 : float or pint.Quantity
            The length of the second (horizontal) rectangle.
        angle : float or pint.Quantity, optional
            Angle in degrees defining the rotation of the cross from the cross's intersection point. Default is 0.
        ref_point : tuple of float or Coordinates, optional
            The reference point that determines the position of the cross's intersection. Default is (0, 0).
        shift_1 : float or pint.Quantity, optional
            The shift in the vertical direction of the first rectangle. Default is 0.
        shift_2 : float or pint.Quantity, optional
            The shift in the horizontal direction of the second rectangle. Default is 0.

        Returns
        -------
        SimpleCross
            A SimpleCross instance.

        Raises
        ------
        ValueError
            If the value of shift_1 or shift_2 is too large, and the two rectangles of the cross do not overlap.

        """

        if shift_1 > length_1 / 2. + width / 2.:
            raise ValueError(f'shift_1 {shift_1} is too large, and the two rectangles of the cross do not overlap.')
        if shift_2 > length_2 / 2. + width / 2.:
            raise ValueError(f'shift_2 {shift_2} is too large, and the two rectangles of the cross do not overlap.')

        if isinstance(ref_point, Tuple):
            ref_point = Coordinates(*ref_point)

        angle = qty_in_uu(angle) if isinstance(angle, Qty) else v_in_uu(angle, 'angle')

        centroid_1 = Coordinates(0, shift_1).get_rotated(angle)
        centroid_2 = Coordinates(shift_2, 0).get_rotated(angle)

        rect_1 = SimpleRectangle.from_centroid(width, length_1, angle, centroid_1)
        rect_2 = SimpleRectangle.from_centroid(width, length_2, angle + Qty(-90, 'deg'), centroid_2)

        cross_coords = coords_for_simple_polygon_merge(rect_1, rect_2, Coordinates(0, 0))
        cross_coords = cross_coords.get_translated(ref_point)

        return cls(cross_coords)
