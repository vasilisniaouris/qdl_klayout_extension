"""
This module contains classes that help create pattern arrays easily.
"""


from typing import Tuple, List

import pya

from qdl_klayout_extension.constants import num_ext
from qdl_klayout_extension.core.coordinates import Coordinates, CoordinatesList


class Pattern:
    """
    A pattern class for generating an array of coordinates in a rectangular pattern.
    """

    def __init__(self, dist_x: num_ext, dist_y: num_ext, n_x: int, n_y: int, angle: num_ext = 0,
                 global_displacement: Coordinates | Tuple[num_ext, num_ext] = (0, 0)):
        """
        Parameters
        ----------
        dist_x : int | float | pint.Quantity
            The distance between each pattern element in the x-direction in user units.
        dist_y : int | float | pint.Quantity
            The distance between each pattern element in the y-direction in user units.
        n_x : int
            The number of pattern elements in the x-direction.
        n_y : int
            The number of pattern elements in the y-direction.
        angle : int | float | pint.Quantity, optional
            The angle in user units by which to rotate the pattern array. Default is 0.
        global_displacement : Coordinates | Tuple[int | float | pint.Quantity, int | float | pint.Quantity], optional
            A tuple of x and y global displacement values, in the form of a Coordinates object or a tuple of two
            int | float | pint.Quantity values in user units.
            Default is (0, 0).
        """

        if isinstance(global_displacement, Tuple):
            global_displacement = Coordinates(*global_displacement)

        self._dist_x = dist_x
        self._dist_y = dist_y
        self._n_x = n_x
        self._n_y = n_y
        self._angle = angle
        self._global_displacement = global_displacement

        self.__post_init__()

    def __post_init__(self):
        self._separation_vx = Coordinates(self.dist_x, 0).get_rotated(self.angle)
        self._separation_vy = Coordinates(0, self.dist_y).get_rotated(self.angle)
        self._klayout_vx = self._separation_vx.klayout_vector
        self._klayout_vy = self._separation_vy.klayout_vector

        self._global_displacement_trans = pya.Trans(self.global_displacement.x_dbu, self.global_displacement.y_dbu)

    def get_pattern_array(self) -> List[List[Coordinates]]:
        """
        The coordinates array of the pattern.

        Returns
        -------
        List[List[Coordinates]]
            A list of lists of Coordinates objects. Each element holds a column and each column holds a
            Coordinates object.
        """

        coords = []
        for i in range(self.n_x):
            column = []
            for j in range(self.n_y):
                x = self.separation_vx.x_uu * i + self.separation_vy.x_uu * j + self.global_displacement.x_uu
                y = self.separation_vx.y_uu * i + self.separation_vy.y_uu * j + self.global_displacement.y_uu
                column.append(Coordinates(x, y))
            coords.append(column)

        return coords

    @property
    def dist_x(self) -> num_ext:
        """ The x distance between each element in the pattern in user units. """
        return self._dist_x

    @property
    def dist_y(self):
        """ The x distance between each element in the pattern in user units. """
        return self._dist_y

    @property
    def n_x(self):
        """ The number of pattern elements in the x-direction. """
        return self._n_x

    @property
    def n_y(self):
        """ The number of pattern elements in the y-direction. """
        return self._n_y

    @property
    def angle(self):
        """ The angle in user units by which to rotate the pattern array. """
        return self._angle

    @property
    def global_displacement(self):
        """ Global displacement of the pattern. """
        return self._global_displacement

    @property
    def separation_vx(self):
        """ Separation vector along the x-axis of the pattern. Pattern array rotation has been applied. """
        return self._separation_vx

    @property
    def separation_vy(self):
        """ Separation vector along the y-axis of the pattern. Pattern array rotation has been applied. """
        return self._separation_vy

    @property
    def klayout_vx(self):
        """ Separation vector along the x-axis of the pattern as a KLayout Vector. """
        return self._klayout_vx

    @property
    def klayout_vy(self):
        """ Separation vector along the y-axis of the pattern as a KLayout Vector. """
        return self._klayout_vy

    @property
    def global_displacement_trans(self):
        """ Global displacement of the pattern in KLayout transformation matrix format. """
        return self._global_displacement_trans

    def get_cell_array_instance(self, cell: pya.Cell):
        """
        Create a new cell instance array.

        Parameters
        ----------
        cell: pya.Cell
            The cell to create a new instance array for, with the given pattern array.

        Returns
        -------
        pya.CellInstArray
            The new cell instance array.
        """
        return pya.CellInstArray(cell.cell_index(), self.global_displacement_trans, self.klayout_vx,
                                 self.klayout_vy, self.n_x, self.n_y)

    def __repr__(self):
        return f'Pattern(Δx={self.dist_x}, Δy={self.dist_y}, n_x={self.n_x}, n_y={self.n_y}, ' \
               f'angle={self.angle}, global_displacement={self.global_displacement})'

