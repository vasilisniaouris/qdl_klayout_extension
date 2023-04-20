"""
This module defines constants and types used throughout this package.

Attributes
----------
user_units: Dict[str, str]
    This dictionary defines the default user units for various physical types. Change this before creating any
    objects from this package, to be consistent between all objects. Changing it after declaring an object,
    it will not update the object's definition and may result in errors.
DBU_UM: float
    Database unit of the produced layout in micrometers. In this case, 1 um would be 1000 database units.
    Change this before creating any objects from this package, to be consistent between all objects.
    Changing it after declaring an object, it will not update the object's definition and may result in errors.
num: type
    A union of a `float` or `int`. This type is used to represent values
    that are a simple numerical value.
num_ext: type
    A union of a `float` or `int` with a `pint.Quantity`. This type is used to represent values
    that can be either a simple numerical value or a more complex quantity with a unit, depending on context.
    Usually, when a function takes a num_ext as input, it will automatically convert the `int` or `float` to
    pint.Quantity, using the `user_units` dictionary.
multi_num: type
    A union of a single numerical value (`float` or `int`) and a `list` of numerical values.
multi_qty: type
    A union of a single `pint.Quantity` and a `list` of `pint.Quantity` objects.
multi_num_ext: type
    A union of a single `num_ext` value and a `list` of `num_ext` values.
"""


from pint import Quantity as Qty
from typing import Sequence

user_units = {'length': 'um', 'angle': 'deg'}
DBU_UM = 0.001

# declaring types
num = float | int
num_ext = num | Qty

multi_num = num | Sequence[num]
multi_qty = Qty | Sequence[Qty]
multi_num_ext = num_ext | Sequence[num_ext]
