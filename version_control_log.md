# Version Control Log

- 04/20/2022; v0.1.0.1: Fixed build issues.
  - Changed `setup.cfg` to follow the [setuptools](https://setuptools.pypa.io/en/latest/userguide/declarative`config.html)
    guide.
  - Simplified `setup.py`.
  - Updated `__init__.py` files to import modules from the same directory.

- 04/20/2022; v0.1.0: Initial commit.
  - Defined basic `utils`, `constants` and `errors` used throughout the package.
  - Defined core class `Coordinates`, `Line`, and `CoordinatesList`.
  - Defined basing geometric shapes.
  - Defined `SimplePolygon` and `Patterns` classes that are heavily dependent on `Coordinates`-based classes.
  - Defined basing geometric shapes to assist in `SimplePolygon` creation.
  - Defined Layout, an extension of KLayout.Layout that assists with automatic database units assignment.
  - Defined coordinate transformation methods.
  - Defined visualization methods for shapes and patterns based on `CoordinatesList`.
  - Provided with working examples of the current functionality. 
