"""
This module contains the KLayout Layout class that helps with automatic database unit assignment.
"""

from typing import Sequence
from types import MethodType
import pya

from qdl_klayout_extension.constants import DBU_UM
from qdl_klayout_extension.core.shapes import Shape


def easy_insert(cell: pya.Cell, layer: int, elements: Sequence[Shape | pya.SimplePolygon | pya.Path]):
    """
    Insert the given elements into the specified layer of the cell.

    Parameters
    ----------
    cell : pya.Cell
        The cell to insert the elements into.
    layer : int
        The layer to insert the elements into.
    elements : Sequence[Shape | pya.SimplePolygon | pya.Path]
        The elements to insert into the cell.
    """
    for element in elements:
        if isinstance(element, Shape):
            element = element.klayout_object
        cell.shapes(layer).insert(element)


class Layout(pya.Layout):
    """
    Custom Layout class inherited from pya.Layout with additional methods and attributes.
    Sets the default database unit of the layout.dbu to the DBU_UM value found in qdl_klayout_extension.constants.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.dbu = DBU_UM

    def update_dbu(self):
        """ Updates the database unit of the layout to the current DBU_UM from qdl_klayout_extension.constants. """
        self.dbu = DBU_UM

    def create_cell(self, *args, **kwargs):
        """ Create a new cell and add the easy_insert() method to it. """
        cell: pya.Cell = super().create_cell(*args, **kwargs)
        cell.easy_insert = MethodType(easy_insert, cell)
        return cell

