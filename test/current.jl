using VoronoiDelaunay
using Random
Random.seed!(1337)

function create_tess(n)
    tess = DelaunayTessellation(n)

    width = max_coord - min_coord
    a = Point2D[Point(min_coord + rand() * width, min_coord + rand() * width) for i in 1:n]
    push!(tess, a)
    return tess
end

function iterate_over_tess_channel(tess)
    i = 0
    for edge in voronoiedges(tess)
        i += 1
    end
end


function iterate_over_tess_iterator(tess)
    i = 0
    for edge in voronoiedges_iterator(tess)
        i += 1
    end
end