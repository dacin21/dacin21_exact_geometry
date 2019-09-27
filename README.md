# dacin21_exact_geometry

Basic routines for computational geometry operating on integers only.

Code has not been tested much yet.

## Features

- Compile time fixed size multiprecision
- 2D geometry stuff with point and angles.
- 2D convex hull, Minkowski sum of convex polygons
- 2D randomized incremental Delaunay triangulation
- Plotting stuff to .svg

## Stuff planned

- 3D convex hull with gift wrapping
- 3D convex hull with divide and conquer
- Segment intersection sweep-line
- Faster multiprecision

## versions

- The branch `with_int128` requires `__int128`, so it won't work on sites such as codeforces.
