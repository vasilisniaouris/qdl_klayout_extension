"""
The errors module contains custom exceptions used in this package.
"""

from typing import Sequence

from pint import Quantity as Qty


class XYLengthError(ValueError):
    """ Exception raised when the length of x and y sequences do not match. """
    def __init__(self, x: Sequence, y: Sequence):
        super().__init__(f'The two arrays {x} and {y} are not of the same length, '
                         f'with lengths {len(x)} and {len(y)}, respectively')


class OperationTypeError(TypeError):
    """ Exception raised when the attempted operation is not defined between two objects. """
    def __init__(self, op1, op2, symbol):
        super().__init__(f"Unsupported operand type(s) for {symbol}: {str(type(op1)).replace('class ', '')[1:-1]} and "
                         f"{str(type(op2)).replace('class ', '')[1:-1]}")


class ShapePointCountError(ValueError):
    """ Exception raised when a shape was initialized with a different amount of points than the expected one. """
    def __init__(self, shape: str, target_length: int, given_length: int):
        super().__init__(f"{shape} needs to {target_length} points, not {given_length} points.")


class ShapeEdgeLengthError(ValueError):
    """ Exception raised when the lengths of a shape's edges do not match the definition of the shape. """
    def __init__(self, shape: str, edge_lengths: Sequence | Qty):
        super().__init__(f"The provided coordinates yield a shape with edge lengths of {edge_lengths}, "
                         f"that do not much the criteria for a {shape}.")


class ShapePointAngleError(ValueError):
    """ Exception raised when the angles between adjacent edges do not match the definition of the shape. """
    def __init__(self, shape: str, point_angles: Sequence | Qty):
        super().__init__(f"The provided coordinates yield a shape with relative point angles of {point_angles}, "
                         f"that do not much the criteria for a {shape}.")


