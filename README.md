# dacin21_exact_geometry

Basic routines for computational geometry operating on integers only.

Code has not been tested much yet.

## Features

- Compile time fixed size multiprecision
- 2D convex hull, Minkowski sum of convex polygons
- 2D randomized incremental Delaunay triangulation
- Plotting stuff to .svg

## Stuff planned

- 3D convex hull with gift wrapping
- 3D convex hull with divide and conquer
- Segment intersection sweep-line
- Faster multiprecision

## versions

- ``geom_64`` is faster, but it requries ``__int128``. 
- ``geom_32`` does not need ``__int128``. This is useful for sites such as codeforces
  which compile on a 32-bit system.
