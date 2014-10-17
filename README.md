# VoronoiDelaunay

[![Build Status](https://travis-ci.org/skariel/VoronoiDelaunay.jl.svg?branch=master)](https://travis-ci.org/skariel/VoronoiDelaunay.jl)
[![Coverage Status](https://img.shields.io/coveralls/skariel/VoronoiDelaunay.jl.svg)](https://coveralls.io/r/skariel/VoronoiDelaunay.jl)

#WIP WIP WIP
Fast, robust construction of 2D Delaunay and Voronoi tessellations on generic point types.
Implementation follows algorithms described in the [Arepo paper](http://arxiv.org/abs/0901.4107)
and used (for e.g.) in the [Illustris Simulation](http://www.illustris-project.org/). License: MIT. Bug reports welcome!

How does it work?
--------------------
Incrementally insert points to a valid Delaunay tessallation, while restoring Delaunayhood by flipping triangles.
Point location (i.e. which triangle should it devide into three) is accelerated by spatial sorting.
Spatial sorting allows to add points which are close in space thus walking the tesselation is fast.
Initial tessalletion includes two triangles built by 4 points which are outside of the allowed region for users.
These "external" triangles are skipped when iterating over Delaunay/Voronoy edges. Fast and robust predicates are
provided by the [GeometricalPredicates](https://github.com/skariel/GeometricalPredicates.jl) package. Benchmarks suggest this package is a bit faster than CGAL, see [here](https://gist.github.com/skariel/3d2018f9341a058e00fc) benchmark of an older version which is also a bit slower than current.

Current limitations
--------------------
* Due to numerical restrictions the point coordinates must be within `min_coord <= x <= max_coord` where `min_coord=1.0+eps(Float64)` and `max_coord=2.0-2eps(Float64)`. Note this is a bit different than what is required by the  `GeometricalPredicates` package.
* Only 2D, but 3D is in the todo list!

How to use?
--------------
###Installation
```Julia
Pkg.add("https://github.com/skariel/VoronoiDelaunay.jl")
```

#WIP WIP WIP
