# Change Log

- 4/24/2023; v0.1.2.0: Visualization updates.
  - Restructured visualizations:
    - Removed `plot_coords` and replaced with `plot_polygon` and `plot_path` that use `SimplePolygon` and `SimplePath`
      objects instead of `CoordinatesList`s.
    - Defined `get_xy_from_coords` to use across plotting methods.
    - Defined `Line2DDataUnits` to `plot_path` easily.
    - Modified `plot_pattern` to take as input a `Pattern` and a shape or a 2D array of shapes that muches the Pattern 
      array.
  - `Pattern.get_pattern_array()` now returns a 2D array of `Coordinates`, instead of a `CoordinatesList`.
  - Updated `geometries` to return sorted `CoordinatesList`s. This fixes the warning issue you get with the definition
    of a `CoordinatesList`, but creates a new issue with the actual `CoordinatesList` elements. There needs to be a
    specified rotationally sorting angle, so that we can overcome this issue.
  - Added transformation methods to `Shape`. 
  - Updated examples to use new visualization methods.

- 4/22/2023; v0.1.1.0: Created new shapes.
  - Added a `Shape` superclass for the existing `SimplePolygon` class and the new `SimplePath` class.
  - Added new geometries (circles, rings, arcs).
  - Moved line segment length calculation for `Line` within the class for consistency.
  - Added `distance_to_point()` method in `Line` to accommodate SimplePath.contains_point().

- 4/21/2023; v0.1.0.2: Minor fixes and refactorings.
  - Changed setup to use `pyproject.toml`, updated `setup.cfg` and `setup.py`.
  - Updated conf.py to include source code links in the API reference (hopefully).
  - Updated `__init__.py` files to import modules from the same directory (successfully this time).

- 04/20/2023; v0.1.0.1: Fixed build issues.
  - Changed `setup.cfg` to follow the [setuptools](https://setuptools.pypa.io/en/latest/userguide/declarative`config.html)
    guide.
  - Simplified `setup.py`.
  - Updated `__init__.py` files to import modules from the same directory.

- 04/20/2023; v0.1.0: Initial commit.
  - Defined basic `utils`, `constants` and `errors` used throughout the package.
  - Defined core class `Coordinates`, `Line`, and `CoordinatesList`.
  - Defined basing geometric shapes.
  - Defined `SimplePolygon` and `Patterns` classes that are heavily dependent on `Coordinates`-based classes.
  - Defined basing geometric shapes to assist in `SimplePolygon` creation.
  - Defined Layout, an extension of KLayout.Layout that assists with automatic database units assignment.
  - Defined coordinate transformation methods.
  - Defined visualization methods for shapes and patterns based on `CoordinatesList`.
  - Provided with working examples of the current functionality. 
