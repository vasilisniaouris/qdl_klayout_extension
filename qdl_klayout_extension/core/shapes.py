"""
This module contains classes that help create various KLayout.Polygon shapes easily.
"""
import warnings
from typing import Tuple, Sequence, List

import numpy as np
import pya
from multipledispatch import dispatch
from pint import Quantity as Qty
from matplotlib import path as mpl_path

from qdl_klayout_extension.constants import num_ext, multi_num_ext, num
from qdl_klayout_extension.core.coordinates import CoordinatesList, Coordinates, Line
from qdl_klayout_extension.core.geometries import rectangle_coords, circle_coords, arc_coords
from qdl_klayout_extension.errors import ShapePointCountError, ShapePointAngleError, ShapeEdgeLengthError
from qdl_klayout_extension.utils import qty_in_uu, v_in_uu, ulu2dbu, dbu2ulu

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


class Shape:
    """
    Base class for all KLayout shapes, such as KLayout.SimplePolygon and KLayout.Path.
    It is used for representing and manipulating shapes in two-dimensional space.

    It is defined by a list of coordinates.

    It provides definitions of various geometric properties of the shape coordinates, such as edge centers, edge angles
    to the y-axis, edge lines and adjacent edge angles.
    It also provides methods for transforming the shape coordinates (and not other initializing parameters!),
    such as translation, rotation and scaling.
    Additionally, it provides the contain_point method, great for checking whether a point or set of points are inside
    or outside the shape. This method needs to be implemented by the subclass.
    """

    _size: int | None = None

    def __init__(self, coords: CoordinatesList, **kwargs):
        """
        Parameters
        ----------
        coords : CoordinatesList
            List of coordinates.
        kwargs
            Additional keyword arguments for the subclasses.
        """

        # if not coords.is_sorted:  # Might not work well with all types of shapes.
        #     coords.sort_clockwise()

        self._coords: CoordinatesList = coords
        self._set_re_init_kwargs(kwargs)
        self.__post_init__()
        self._check_all()

    def __post_init__(self):
        self._centroid: Coordinates = self._coords.get_centroid()
        self._ref_point: Coordinates = self._coords.ref_point

        self._edge_centers: CoordinatesList = self._coords.get_edge_centers()
        self._edge_lines: List[Line] = self._coords.get_edge_lines()
        self._edge_lengths = self._coords.get_edge_lengths()

        self._edge_angles = self._coords.get_edge_angles()
        self._point_angles = self._coords.get_point_angles()

        if self._size is None:
            self._size = len(self._coords)

    def _set_re_init_kwargs(self, kwargs):
        """
        Helper method to set the re-initialization keyword arguments.
        These are used in case the user wants to transform the shape.
        In that case the transformation should change the coordinate properties, but not any other initial properties.
        """
        self._re_init_kwargs = kwargs

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
    def klayout_object(self):
        """ KLayout object representing the shape. Defaults to None and must be redefined when Shape is subclassed. """
        return None

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

    @dispatch(object)
    def get_translated(self, coords: Tuple | Coordinates | Sequence | "CoordinatesList"):
        """
        Translate the Shape coordinates by given values in Coordinates.

        Parameters
        ----------
        coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The distance to translate the CoordinatesList.

        Returns
        -------
        Shape
            The translated shape.
        """
        coords = self.coords.get_translated(coords)
        return self.__class__(coords, **self._re_init_kwargs)

    @dispatch(object, object)
    def get_translated(self, x_coords: multi_num_ext, y_coords: multi_num_ext):
        """
        Translate the Shape coordinates by given x and y values.

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
        Shape
            The translated Shape.
        """
        coords = self.coords.get_translated(x_coords, y_coords)
        return self.__class__(coords, **self._re_init_kwargs)

    def get_transformed(self, trans_matrix: np.ndarray):
        """
        Transform the shape coordinates using the given transformation matrix.

        Parameters
        ----------
        trans_matrix : np.ndarray
            The transformation matrix.

        Returns
        -------
        Shape
            The transformed shape.
        """
        coords = self.coords.get_transformed(trans_matrix)
        return self.__class__(coords, **self._re_init_kwargs)

    def get_rotated(self, angle: num_ext, counterclockwise=True):
        """
        Rotate the shape coordinates by a given angle around the origin.

        Parameters
        ----------
        angle : int | float | pint.Quantity
            The angle to rotate by in degrees.
        counterclockwise : bool, optional
            If True (default), rotate counterclockwise. Otherwise, rotate clockwise.

        Returns
        -------
        Shape
            The rotated shape.
        """
        coords = self.coords.get_rotated(angle, counterclockwise)
        return self.__class__(coords, **self._re_init_kwargs)

    def get_stretched(self, x_stretch: num, y_stretch: num):
        """
        Stretch (enlarge or compress) the shape coordinates by given x and y factors.

        Parameters
        ----------
        x_stretch : int | float
            The factor to stretch the x-coordinates by.
        y_stretch : int | float
            The factor to stretch the y-coordinates by.

        Returns
        -------
        Shape
            The stretched shape coordinates.
        """
        coords = self.coords.get_stretched(x_stretch, y_stretch)
        return self.__class__(coords, **self._re_init_kwargs)

    def get_squeezed(self, squeeze_value: num):
        """
        Squeeze the shape coordinates by given factor. Squeezing conserves the shape's area.

        Parameters
        ----------
        squeeze_value : int | float
            The factor to squeeze the coordinates by.

        Returns
        -------
        Shape
            The squeezed shape.
        """
        coords = self.coords.get_squeezed(squeeze_value)
        return self.__class__(coords, **self._re_init_kwargs)

    @dispatch(object)
    def get_reflected(self, coords: Tuple | Coordinates | Sequence | "CoordinatesList"):
        """
        Reflect the shape coordinates across a line defined by the given coordinates (lx, ly) and the origin (0, 0).

        Parameters
        ----------
        coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The coordinates that define the reflection line.

        Returns
        -------
        Shape
            The reflected shape.
        """
        coords = self.coords.get_reflected(coords)
        return self.__class__(coords, **self._re_init_kwargs)

    @dispatch(object, object)
    def get_reflected(self, x_coords: multi_num_ext, y_coords: multi_num_ext):
        """
        Reflect the shape coordinates across a line defined by the given coordinates (lx, ly) and the origin (0, 0).

        Parameters
        ----------
        x_coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The x coordinates in user units that define the reflection line.
        y_coords : Tuple | Coordinates | Sequence | "CoordinatesList"
            The y coordinates in user units that define the reflection line.

        Returns
        -------
        Shape
            The reflected shape.
        """
        coords = self.coords.get_reflected(x_coords, y_coords)
        return self.__class__(coords, **self._re_init_kwargs)

    def get_reflected_x(self):
        """
        Reflect the shape coordinates across the x-axis.

        Returns
        -------
        Shape
            The reflected shape.
        """
        coords = self.coords.get_reflected_x()
        return self.__class__(coords, **self._re_init_kwargs)

    def get_reflected_y(self):
        """
        Reflect the shape coordinates across the y-axis.

        Returns
        -------
        Shape
            The reflected shape.
        """
        coords = self.coords.get_reflected_x()
        return self.__class__(coords, **self._re_init_kwargs)

    def contains_point(self, point: Coordinates):
        """
        Check if a point is contained within the shape. Must be overloaded when subclassed.

        Parameters
        ----------
        point : Coordinates
            The point to check if it's contained in the polygon.

        Returns
        -------
        bool
            True if the point is contained within the shape, False otherwise.
        """
        raise NotImplementedError('contains_point must be defined by the subclass.')

    def contains_points(self, points: CoordinatesList | Sequence[Coordinates]):
        """
        Check if a list of points are contained within the shape.

        Parameters
        ----------
        points : CoordinatesList | Sequence[Coordinates]
            The list of points to check if they are contained in the polygon.

        Returns
        -------
        List[bool]
            A list of booleans indicating if each point is contained within the shape.
        """
        return [self.contains_point(point) for point in points]


class SimplePolygon(Shape):
    """
    The SimplePolygon class is a subclass of the Shape class.
    It is used for representing and manipulating simple polygons in two-dimensional space.
    It defines a polygon as a list of coordinates that represent the vertices of the polygon in order.
    The polygon coordinates must be sorted appropriate. Most of the time, using the sort_clockwise() method of the
    coord class with a properly specified ref_point will ensure that the polygon will be visualized properly.

    The SimplePolygon class serves as a foundation for more specialized polygon classes, such as the SimpleRectangle
    and SimpleCross classes, which inherit from it and provide additional functionality specific to those shapes.
    """

    def __post_init__(self):
        self._klayout_simple_polygon: pya.SimplePolygon = pya.SimplePolygon(self._coords.klayout_points)
        super().__post_init__()

    @property
    def klayout_object(self) -> pya.SimplePolygon:
        """ KLayout SimplePolygon object representing the polygon vertices. """
        return self._klayout_simple_polygon

    @property
    def klayout_simple_polygon(self) -> pya.SimplePolygon:
        """ KLayout SimplePolygon object representing the polygon vertices. """
        return self._klayout_simple_polygon

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


class SimpleCircle(SimplePolygon):
    """
    The SimpleCircle class is a subclass of the SimplePolygon class that defines a simple circle.
    The from_center() method constructs a SimpleCircle instance from the parameters that define the circle.
    The reference point of this shape is the center of the circle.

    Warning
    -------
    Τhe SimpleCircle.from_center() is designed to provide you with overlapping coordinates, to alleviate edge issues.
    Since the is_sorted() method of the CoordinatesList class can only see if 1 period is sorted, it can not tell if
    there are more than one period in a given coordinate list.
    For similar reasons behind sort_rotationally(), trying to re-sort a SimpleCircle will yield a non-overlapped circle!
    """

    @classmethod
    def from_center(cls, radius: num_ext, num_points: int = 100,
                    center: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0)):
        """
        Constructs a SimpleCircle instance from the radius of the circle and the number of points.

        Parameters:
        -----------
        radius: int | float | pint.Quantity
            The radius of the circle.
        num_points: int
            The number of points used to construct the circle. Default is 100.
        center: Coordinates or Tuple[int | float | pint.Quantity, int | float | pint.Quantity], optional
            The center of the circle. Default is Coordinates(0, 0).

        Returns:
        --------
        SimpleCircle
            A SimpleCircle object constructed from the provided parameters.
        """
        if isinstance(center, Tuple):
            center = Coordinates(*center)

        coords = circle_coords(radius, num_points)
        coords = coords.get_translated(center)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # To catch the sorting warning
            return cls(coords)


class SimplePath(Shape):
    """
    The SimplePath class is a subclass of the Shape class that defines a path with a given width.
    It is used for representing and manipulating paths in two-dimensional space.
    It defines a path as a list of coordinates that represent the vertices of the path in order.
    The path coordinates must be sorted appropriate. Most of the time, using the sort_clockwise() method of the
    coord class with a properly specified ref_point will ensure that the path will be visualized properly.
    SimplePath transformations only transform the path coordinates, not the width or end caps.

    The SimplePath class serves as a foundation for more specialized polygon classes, such as the RingPath and ArcPath
    classes, which inherit from it and provide additional functionality specific to those shapes.
    """

    def __init__(self, coords: CoordinatesList, width: num_ext, initial_end_cap: num_ext = 0,
                 final_end_cap: num_ext = 0, rounded_end_caps: bool = False):
        """
        Parameters:
        -----------
        coords: CoordinatesList
            A list of coordinates representing the vertices of the path.
        width: int | float | pint.Quantity
            The width of the path.
        initial_end_cap: int | float | pint.Quantity, optional
            The extent of the initial end cap. Default is 0.
        final_end_cap: int | float | pint.Quantity, optional
            The extent of the final end cap. Default is 0.
        rounded_end_caps: bool, optional
            The boolean value that indicates whether the end caps are rounded or not. Default is False.
        """

        self._width_uu: Qty = qty_in_uu(width) if isinstance(width, Qty) else v_in_uu(width)
        self._rounded_end_caps: bool = rounded_end_caps
        if self._rounded_end_caps:
            end_cap_width = dbu2ulu(ulu2dbu(self._width_uu / 2))
            # The caps must be half the width, so, to avoid digitization issues, we slightly change the width definition
            self._width_uu = end_cap_width * 2
            self._initial_end_cap_uu = end_cap_width
            self._final_end_cap_uu = end_cap_width
        else:
            self._initial_end_cap_uu = qty_in_uu(initial_end_cap) if isinstance(initial_end_cap, Qty) \
                else v_in_uu(initial_end_cap)
            self._final_end_cap_uu = qty_in_uu(final_end_cap) if isinstance(final_end_cap, Qty) \
                else v_in_uu(final_end_cap)

        kwargs = {'width': self.width_uu, 'initial_end_cap': self.initial_end_cap_uu,
                  'final_end_cap': self.final_end_cap_uu, 'rounded_end_caps': self.rounded_end_cap}
        super().__init__(coords, **kwargs)

    def __post_init__(self):
        self._klayout_path: pya.Path = pya.Path(self.coords.klayout_points, self.width_dbu, self.initial_end_cap_dbu,
                                                self.final_end_cap_dbu, self.rounded_end_cap)
        super().__post_init__()

    @property
    def width_uu(self) -> Qty:
        """ The width of the path in user units. """
        return self._width_uu

    @property
    def initial_end_cap_uu(self) -> Qty:
        """ The extent of the initial end cap in user units. """
        return self._initial_end_cap_uu

    @property
    def final_end_cap_uu(self) -> Qty:
        """ The extent of the final end cap in user units. """
        return self._final_end_cap_uu

    @property
    def width_dbu(self) -> Qty:
        """ The width of the path in database units. """
        return ulu2dbu(self._width_uu)

    @property
    def initial_end_cap_dbu(self) -> Qty:
        """ The extent of the initial end cap in database units. """
        return ulu2dbu(self._initial_end_cap_uu)

    @property
    def final_end_cap_dbu(self) -> Qty:
        """ The extent of the final end cap in database units. """
        return ulu2dbu(self._final_end_cap_uu)

    @property
    def rounded_end_cap(self) -> bool:
        """ The boolean value that indicates whether the end caps are rounded or not. """
        return self._rounded_end_caps

    @property
    def klayout_object(self) -> pya.Path:
        """ KLayout Path object representing the path vertices. """
        return self._klayout_path

    @property
    def klayout_path(self) -> pya.Path:
        """ KLayout Path object representing the path vertices. """
        return self._klayout_path

    def get_end_cap_coords(self) -> Tuple[Coordinates, Coordinates]:
        """
        Get the coordinates of the initial and final end caps.
        The coordinates of the initial and final end caps are calculated by finding the first and last coordinates
        of the path, and then adding the width to each of them.

        Returns
        -------
        Tuple[Coordinates, Coordinates]
            The coordinates of the initial and final end caps.
        """

        x = self.coords.x_uu
        y = self.coords.y_uu
        width_init = self.initial_end_cap_uu
        width_final = self.final_end_cap_uu

        dx_init = x[1] - x[0]
        dy_init = y[1] - y[0]
        length_init = np.sqrt(dx_init ** 2 + dy_init ** 2)
        x_init = x[0] - width_init * dx_init / length_init
        y_init = y[0] - width_init * dy_init / length_init

        dx_final = x[-1] - x[-2]
        dy_final = y[-1] - y[-2]
        length_final = np.sqrt(dx_final ** 2 + dy_final ** 2)
        x_final = x[-1] + width_final * dx_final / length_final
        y_final = y[-1] + width_final * dy_final / length_final

        return Coordinates(x_init, y_init), Coordinates(x_final, y_final)

    def contains_point(self, point: Coordinates):
        """
        Check if a point is contained within the path.

        Parameters
        ----------
        point : Coordinates
            The point to check if it's contained in the path.

        Returns
        -------
        bool
            True if the point is contained within the path, False otherwise.

        Notes
        -----
        In this method, we first find the distance of the point from all line segments in the path.
        If one of these distances is less than the half-width, then the point may be contained in the shape.

        """
        edge_lines = self.edge_lines
        line_point_distances = Qty.from_list([line.distance_from_point(point) for line in edge_lines])

        if not any(line_point_distances <= self.width_uu / 2):
            return False

        # TODO: check if the point is within the edge_line limits.
        # TODO: figure out corner effects -> Probably have to define a get_polygon_coords() method, where you first
        #  find the extended line segments by changing the line intercept, and then taking the infinite line
        #  intersections between lines as the new points.
        return True


class RingPath(SimplePath):
    """
    The RingPath class is a subclass of the SimplePath class that defines a ring with a given width.
    The from_center() method constructs a SimplePath instance from the parameters that define the circle.
    The reference point of this shape is the center of the circle.

    Warning
    -------
    The ring path will yield a non-sorted warning for the edge and line calculations. Ignore it.
    The is_sorted() method of the CoordinatesList class can only see if 1 period is sorted,
    and
    Τhe RingPath.from_center() is designed to provide you with overlapping coordinates, to alleviate edge issues.
    Since the is_sorted() method of the CoordinatesList class can only see if 1 period is sorted, it can not tell if
    there are more than one period in a given coordinate list.
    For similar reasons behind sort_rotationally(), trying to re-sort a RingPath will yield a non-overlapped ring!
    """

    @classmethod
    def from_center(cls, radius: num_ext, width: num_ext, num_points: int = 100,
                    center: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0)):
        """
        Constructs a RingPath instance from the radius of the ring and the number of points.

        Parameters:
        -----------
        radius: int | float | pint.Quantity
            The radius of the ring.
        width: int | float | pint.Quantity
            The width of the ring.
        num_points: int
            The number of points used to construct the ring. Default is 100.
        center: Coordinates or Tuple[int | float | pint.Quantity, int | float | pint.Quantity], optional
            The center of the ring. Default is Coordinates(0, 0).

        Returns:
        --------
        RingPath
            A RingPath object constructed from the provided parameters.
        """
        if isinstance(center, Tuple):
            center = Coordinates(*center)

        coords = circle_coords(radius, num_points)
        coords = coords.get_translated(center)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)  # To catch the sorting warning
            return cls(coords, width)


class ArcPath(SimplePath):
    """
    The ArcPath class is a subclass of the SimplePath class that defines an arc with a given width.
    The from_center() method constructs a SimplePath instance from the parameters that define the arc.
    The reference point of this shape is the center of the arc.
    """

    @classmethod
    def from_center(cls, radius: num_ext, width: num_ext, num_points: int = 100,
                    center: Coordinates | Tuple[num_ext, num_ext] = Coordinates(0, 0),
                    start_angle: num_ext = 0, angle_range: num_ext = Qty(90, 'deg')):
        """
        Constructs an ArcPath instance from the radius of the arc and the number of points.
        If the angle range is greater than a full circle, it returns a RingPath instead.

        Parameters:
        -----------
        radius: int | float | pint.Quantity
            The radius of the arc.
        width: int | float | pint.Quantity
            The width of the arc.
        num_points: int
            The number of points used to construct the arc. Default is 100.
        center: Coordinates or Tuple[int | float | pint.Quantity, int | float | pint.Quantity], optional
            The center of the arc. Default is Coordinates(0, 0).
        start_angle: int | float | pint.Quantity, optional
            The starting angle of the arc. Default is 0.
        angle_range: int | float | pint.Quantity, optional
            The  angle range of the arc. Default is 90 degrees.

        Returns:
        --------
        ArcPath | RingPath
            An ArcPath object constructed from the provided parameters.
            A RingPath object is constructed if the angle range is greater than a full circle.
        """

        ar = qty_in_uu(angle_range) if isinstance(angle_range, Qty) else v_in_uu(angle_range, 'angle')
        if ar.to('deg').m >= 360:
            return RingPath.from_center(radius, width, num_points, center)

        if isinstance(center, Tuple):
            center = Coordinates(*center)

        coords = arc_coords(radius, start_angle, angle_range, num_points)
        coords = coords.get_translated(center)

        return cls(coords, width)
