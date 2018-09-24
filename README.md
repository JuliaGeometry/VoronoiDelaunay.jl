# VoronoiDelaunay.jl

[![Build Status](https://travis-ci.org/JuliaGeometry/VoronoiDelaunay.jl.svg?branch=master)](https://travis-ci.org/JuliaGeometry/VoronoiDelaunay.jl)
[![VoronoiDelaunay](http://pkg.julialang.org/badges/VoronoiDelaunay_0.6.svg)](http://pkg.julialang.org/detail/VoronoiDelaunay)
[![Coverage Status](https://coveralls.io/repos/JuliaGeometry/VoronoiDelaunay.jl/badge.svg?branch=master)](https://coveralls.io/r/JuliaGeometry/VoronoiDelaunay.jl?branch=master)![Alt VoronoiDelaunay.jl](http://i.imgur.com/lh8VLZ5.png5 "VoronoiDelaunay.jl")

Fast, robust construction of 2D Delaunay and Voronoi tessellations on generic point types.
Implementation follows algorithms described in the [Arepo paper](http://arxiv.org/abs/0901.4107)
and used (for e.g.) in the [Illustris Simulation](http://www.illustris-project.org/). License: MIT. Bug reports welcome!

How does it work?
--------------------
Incrementally insert points to a valid Delaunay tessallation, while restoring Delaunayhood by flipping triangles.
Point location (i.e. which triangle should it divide into three) is accelerated by spatial sorting.
Spatial sorting allows to add points which are close in space thus walking the tesselation is fast.
Initial tessalletion includes two triangles built by 4 points which are outside of the allowed region for users.
These "external" triangles are skipped when iterating over Delaunay/Voronoy edges. Fast and robust predicates are
provided by the [GeometricalPredicates](https://github.com/skariel/GeometricalPredicates.jl) package. Benchmarks suggest this package is a bit faster than CGAL, see [here](https://gist.github.com/skariel/3d2018f9341a058e00fc) benchmark of an older version which is also a bit slower than current.

Current limitations
--------------------
* Due to numerical restrictions the point coordinates must be within `min_coord <= x <= max_coord` where `min_coord=1.0+eps(Float64)` and `max_coord=2.0-2eps(Float64)`. Note this is a bit different than what is required by the  `GeometricalPredicates` package.
* The following features are not implemented, but are in the TODO list; In order of priority: centroid tessellations (Lloy's method), Weighted generators (both power and sum), bounding, maybe restricting. Hierarchal tessellations for fast random locatings; Distributed tessellation construction. 3D. Order of priority may change of course :)

How to use?
--------------
### Installation
```Julia
]add VoronoiDelaunay
```
For Julia 0.6 and below, type
```Julia
Pkg.add("VoronoiDelaunay")
```

### Building a tessellation
Define and push individual points like this:
```Julia
using  VoronoiDelaunay
tess = DelaunayTessellation()
push!(tess, Point(1.5, 1.5))
```
creation of points is explained in the [GeometricalPredicates](https://github.com/skariel/GeometricalPredicates.jl) package documentation.

Pushing arrays of points is more efficient:
```Julia
width = max_coord - min_coord
a= Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:100]
push!(tess, a)
```
notice care taken for correct range of coordinates. `min_coord` and `max_coord` are defined in the package. We can further optimize by giving a `sizehint` at time of construction:
```Julia
tess = DelaunayTessellation(100)
```
or at any later point:
```Julia
sizehint(tess, 100)
```
### Iterating
Delaunay tesselations need at least 3 points to be well defined. Voronoi need 4. Remember this when iterating or plotting.
Iterating over Delaunay edges is done like this:
```Julia
i = 0
for edge in delaunayedges(tess)
    i += 1
    # or, do something more useful :)
end
```
a `DelaunayEdge` contains two points a and b, they can be accesse with `geta(edge)` and `getb(edge)`.
Iterating over Voronoi edges is similar:
```Julia
i = 0
for edge in voronoiedges(tess)
    i += 1
    # or, do something more useful :)
end
```
a `VoronoiEdge` is a bit different than a `DelaunayEdge`: here `a` and `b` are `Point2D` and not the generators, as they have different coordinates. To get the generators use `getgena(edge)` and `getgenb(edge)` these give the relevant `AbstractPoint2D` which were used to create the edge.

If the generators are not needed when iterating over the Voronoi edges (e.g. when plotting) then a more efficient way to iterate is:
```Julia
i = 0
e=Nothing
for edge in voronoiedgeswithoutgenerators(tess)
    i += 1
    # do something more useful here :)
end
```
here `edge` is a `VoronoiEdgeWithoutGenerators`, the points `a` and `b` can be accessed as usual.

Iterating over Delaunay triangles:
```Julia
i = 0
for delaunaytriangle in tess
    i += 1
    # or, do something more useful :)
end
```
`delaunaytriangle` here is of type `DelaunayTriagle` which is s subtype of `AbstractNegativelyOrientedTriangle`. To get the generators of this triangle use the `geta`, `getb`, and `getc` methods. You can do all other operations and predicate tests on this triangle as explained in [GeometricalPredicates](https://github.com/skariel/GeometricalPredicates.jl)

### Navigating
Locating a point, i.e. finding the triangle it is inside:
```Julia
t = locate(tess, Point(1.2, 1.3))
```
if the point is outside of the tessellation then `isexternal(t) == true` holds. This is good for type stability, at least better than returning a `Nothing`. It is assumed that the point we want to locate is actually in the allowed points region. Performance is best when locating points close to each other (this is also why spatial sorting is used). Future versions may implement a hierarchal approach for fast random locations.

Navigating from a triangle to its neighbours is done like this:
```Julia
t = movea(tess, t)  # move to the direction infront of generator a
t = moveb(tess, t)  # move to the direction infront of generator b
t = movec(tess, t)  # move to the direction infront of generator c
```

### Plotting
The following retrieves a couple of vectors ready to plot Voronoi edges:
```Julia
x, y = getplotxy(voronoiedges(tess))
```
and for Delaunay edges:
```Julia
x, y = getplotxy(delaunayedges(tess))
```
Now plotting can be done with your favorite plotting package, for e.g.:
```Julia
using Gadfly
plot(x=x, y=y, Geom.path)
```
To make a nice looking plot remember to limit the axes and aspect ratio. For e.g.:
```Julia
set_default_plot_size(15cm, 15cm)
plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.0, maxvalue=2.0))
```

### From an image
You can create a tesselation from an image, just like the tesselation of the
julia logo at the top of this README. This was created from a png with `from_file`
(see `examples/img_to_vorono.jl`):
```Julia
import Images: imread
img = imread("julia.png")
tess = from_image(img, 25000)
```
