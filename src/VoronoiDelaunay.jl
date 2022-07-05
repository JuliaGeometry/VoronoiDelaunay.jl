module VoronoiDelaunay

# Fast, robust 2D Voronoi/Delaunay tesselation
# Implementation follows algorithms described in http://arxiv.org/abs/0901.4107
# and used (for e.g.) in the Illustris Simulation
# http://www.illustris-project.org/
#
# Author: Ariel Keselman (skariel@gmail.com)
# License: MIT
# Bug reports welcome!

export
DelaunayTessellation, DelaunayTessellation2D, sizehint!, isexternal,
min_coord, max_coord, locate, movea, moveb, movec,
delaunayedges, voronoiedges, voronoiedgeswithoutgenerators,
iterate, findindex, push!,
Point, Point2D, AbstractPoint2D, getx, gety, geta, getb, getc,
getgena, getgenb, getplotxy

using GeometricalPredicates
import GeometricalPredicates: geta, getb, getc

import Base: push!, iterate, copy, sizehint!
import Colors: RGB, RGBA
using Random: shuffle!

const min_coord = GeometricalPredicates.min_coord + eps(Float64)
const max_coord = GeometricalPredicates.max_coord - eps(Float64)

mutable struct DelaunayTriangle{T<:AbstractPoint2D} <: AbstractNegativelyOrientedTriangle
    _a::T; _b::T; _c::T
    _bx::Float64; _by::Float64
    _cx::Float64; _cy::Float64
    _px::Float64; _py::Float64
    _pr2::Float64
    _neighbour_a::Int64
    _neighbour_b::Int64
    _neighbour_c::Int64

    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 na::Int64, nb::Int64, nc::Int64) where T
        t = new(pa, pb, pc, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, na, nb, nc)
        clean!(t)
        t
    end

    # this constructor is good for making copies
    function DelaunayTriangle{T}(pa::T, pb::T, pc::T,
                                 bx::Float64, by::Float64,
                                 cx::Float64, cy::Float64,
                                 px::Float64, py::Float64,
                                 pr2::Float64,
                                 na::Int64, nb::Int64, nc::Int64) where T
        new(pa, pb, pc, bx, by, cx, cy, px, py, pr2, na, nb, nc)
    end
end
function DelaunayTriangle(pa::T, pb::T, pc::T,
                          na::Integer, nb::Integer, nc::Integer) where {T<:AbstractPoint2D}
    DelaunayTriangle{T}(pa, pb, pc, Int64(na), Int64(nb), Int64(nc))
end
function DelaunayTriangle(pa::T, pb::T, pc::T,
                          bx::Float64, by::Float64,
                          cx::Float64, cy::Float64,
                          px::Float64, py::Float64,
                          pr2::Float64,
                          na::Integer, nb::Integer, nc::Integer) where {T<:AbstractPoint2D}
    DelaunayTriangle{T}(pa, pb, pc, bx, by, cx, cy, px, py, pr2, Int64(na), Int64(nb), Int64(nc))
end

function copy(t::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    DelaunayTriangle(
                     t._a, t._b, t._c,
                     t._bx, t._by,
                     t._cx, t._cy,
                     t._px, t._py,
                     t._pr2,
                     t._neighbour_a, t._neighbour_b, t._neighbour_c
                     )
end

function isexternal(t::DelaunayTriangle)
    isexternal(geta(t)) || isexternal(getb(t)) || isexternal(getc(t))
end
function isexternal(p::AbstractPoint2D)
    # TODO: verify if checking only getx(p) is sufficient?
    isexternal(getx(p)) || isexternal(gety(p))
end
function isexternal(z::Real)
    z < min_coord || z > max_coord
end

mutable struct DelaunayTessellation2D{T<:AbstractPoint2D}
    _trigs::Vector{DelaunayTriangle{T}}
    _last_trig_index::Int64
    _edges_to_check::Vector{Int64}
    _total_points_added::Int64

    function DelaunayTessellation2D{T}(n::Int = 100) where {T}
        a = T(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
        b = T(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
        c = T(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
        d = T(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
        t1 = DelaunayTriangle{T}(d, c, b, Int64(2), Int64(1), Int64(1))
        t2 = DelaunayTriangle{T}(a, b, c, Int64(3), Int64(1), Int64(1))
        t3 = DelaunayTriangle{T}(d, c, b, Int64(2), Int64(1), Int64(1))
        _trigs = DelaunayTriangle{T}[t1, t2, t3]
        t = new(_trigs, Int64(3), Int64[], Int64(0))
        sizehint!(t._edges_to_check, 1000)
        sizehint!(t, n)
    end
end
DelaunayTessellation2D(n::Integer) = DelaunayTessellation2D{Point2D}(n)
DelaunayTessellation2D(n::Integer, ::T) where {T<:AbstractPoint2D} = DelaunayTessellation2D{T}(n)
DelaunayTessellation(n::Integer=100) = DelaunayTessellation2D(n)

function sizehint!(t::DelaunayTessellation2D{T}, n::Integer) where {T<:AbstractPoint2D}
    required_total_size = Int64(2n) + Int64(10)
    required_total_size <= length(t._trigs) && return
    sizehint!(t._trigs, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
    end
    t
end

# growing strategy
function sizefit_at_least(t::DelaunayTessellation2D{T}, n::Integer) where T<:AbstractPoint2D
    minimal_acceptable_actual_size = Int64(2n) + Int64(10)
    minimal_acceptable_actual_size <= length(t._trigs) && return
    required_total_size::Int64 = length(t._trigs)
    while required_total_size < minimal_acceptable_actual_size
        required_total_size += required_total_size >>> 1
    end
    sizehint!(t._trigs, required_total_size)
    while length(t._trigs) < required_total_size
        push!(t._trigs, copy(t._trigs[end]))
    end
    t
end

struct DelaunayEdge{T<:AbstractPoint2D}
    _a::T
    _b::T
end
geta(e::DelaunayEdge{T}) where {T<:AbstractPoint2D} = e._a
getb(e::DelaunayEdge{T}) where {T<:AbstractPoint2D} = e._b

struct VoronoiEdge{T<:AbstractPoint2D}
    _a::Point2D
    _b::Point2D
    _generator_a::T
    _generator_b::T
end
geta(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._a
getb(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._b
getgena(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._generator_a
getgenb(e::VoronoiEdge{T}) where {T<:AbstractPoint2D} = e._generator_b

struct VoronoiEdgeWithoutGenerators
    _a::Point2D
    _b::Point2D
end
geta(e::VoronoiEdgeWithoutGenerators) = e._a
getb(e::VoronoiEdgeWithoutGenerators) = e._b

struct DelaunayEdgeIterator{T <: DelaunayTessellation2D}
    t::T
end

Base.IteratorSize(::DelaunayEdgeIterator) = Base.SizeUnknown()
Base.eltype(::DelaunayEdgeIterator{DelaunayTessellation2D{T}}) where T = DelaunayEdge{T}

function iterate(v::DelaunayEdgeIterator, state=(Int64(2), 1))
    ix::Int64, tx = state
    t = v.t
    while ix <= t._last_trig_index
        tr = t._trigs[ix]
        if isexternal(tr)
            ix += 1
            tx = 1
            continue
        end
        if tx == 1
            ix_new = tr._neighbour_a
            if ix_new > ix || isexternal(t._trigs[ix_new])
                return DelaunayEdge(getb(tr), getc(tr)), (ix, 2)
            end
        elseif tx == 2
            ix_new = tr._neighbour_b
            if ix_new > ix || isexternal(t._trigs[ix_new])
                return DelaunayEdge(geta(tr), getc(tr)), (ix, 3)
            end
        else # tx == 3
            ix_new = tr._neighbour_c
            if ix_new > ix || isexternal(t._trigs[ix_new])
                return DelaunayEdge(geta(tr), getb(tr)), (ix+1, 1)
            end
        end
        if tx < 3
            tx += 1
        else
            ix += 1
            tx = 1
        end
    end
    return nothing
end

delaunayedges(t::DelaunayTessellation2D) = DelaunayEdgeIterator(t)

struct VoronoiEdgeIterator{T <: DelaunayTessellation2D}
    t::T
end
Base.IteratorSize(::VoronoiEdgeIterator) = Base.SizeUnknown()
Base.eltype(::VoronoiEdgeIterator{<:DelaunayTessellation2D{T}}) where T = VoronoiEdge{T}
function iterate(v::VoronoiEdgeIterator, state=(Int64(2), 1))
    ix::Int64, tx = state
    t = v.t
    while ix <= t._last_trig_index
        tr = t._trigs[ix]
        cc = circumcenter(tr)

        ix_na = tr._neighbour_a
        if tx == 1 && ix_na > ix
            nb = t._trigs[ix_na]
            return (VoronoiEdge(cc, circumcenter(nb), getb(tr), getc(tr)), (ix, 2))
        end
        ix_nb = tr._neighbour_b
        if tx <= 2 && ix_nb > ix
            nb = t._trigs[ix_nb]
            return (VoronoiEdge(cc, circumcenter(nb), geta(tr), getc(tr)), (ix, 3))
        end
        ix_nc = tr._neighbour_c
        if tx <= 3 && ix_nc > ix
            nb = t._trigs[ix_nc]
            return (VoronoiEdge(cc, circumcenter(nb), geta(tr), getb(tr)), (ix+1, 1))
        end
        tx = 1
        ix += 1
    end
    nothing
end
voronoiedges(t::DelaunayTessellation2D) = VoronoiEdgeIterator(t)

struct VoronoiEdgeIteratorWithoutGenerator{T <: DelaunayTessellation2D}
    t::T
end
Base.IteratorSize(::VoronoiEdgeIteratorWithoutGenerator) = Base.SizeUnknown()
Base.eltype(::VoronoiEdgeIteratorWithoutGenerator{<:DelaunayTessellation2D{T}}) where T = VoronoiEdge{T}
function iterate(v::VoronoiEdgeIteratorWithoutGenerator, state=(Int64(2), 1))
    ix::Int64, tx = state
    t = v.t
    while ix <= t._last_trig_index
        tr = t._trigs[ix]
        cc = circumcenter(tr)

        ix_na = tr._neighbour_a
        if tx <= 1 && ix_na > ix
            nb = t._trigs[ix_na]
            return (VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)), (ix, 2))
        end
        ix_nb = tr._neighbour_b
        if tx <= 2 && ix_nb > ix
            nb = t._trigs[ix_nb]
            return (VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)), (ix, 3))
        end
        ix_nc = tr._neighbour_c
        if tx <= 3 && ix_nc > ix
            nb = t._trigs[ix_nc]
            return (VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)), (ix+1, 1))
        end
        tx = 1
        ix += 1
    end
    nothing
end
function voronoiedgeswithoutgenerators(t::DelaunayTessellation2D)
    VoronoiEdgeIteratorWithoutGenerator(t)
end

# TODO: for v0.5, remove TrigIter
mutable struct TrigIter
    ix::Int64
end

# TODO: for v0.5, replace it by ix::Int
function iterate(t::DelaunayTessellation2D, it::TrigIter=TrigIter(Int64(2)))
    while it.ix <= t._last_trig_index && isexternal(@inbounds t._trigs[it.ix])
        it.ix += 1
    end
    if it.ix > t._last_trig_index
        return nothing
    end
    trig = t._trigs[it.ix]
    it.ix += 1
    return (trig, it)
end

function findindex(tess::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D}
    i = tess._last_trig_index
    while true
        @inbounds w = intriangle(tess._trigs[i], p)
        w > 0 && return i
        @inbounds tr = tess._trigs[i]
        if w == -1
            i = tr._neighbour_a
        elseif w == -2
            i = tr._neighbour_b
        else
            i = tr._neighbour_c
        end
    end
end

locate(t::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D} = t._trigs[findindex(t, p)]

movea(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_a]
moveb(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_b]
movec(tess::DelaunayTessellation2D{T},
      trig::DelaunayTriangle{T}) where {T<:AbstractPoint2D} = tess._trigs[trig._neighbour_c]

function _pushunfixed!(tess::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D}
    i = findindex(tess, p)
    ltrigs1 = tess._last_trig_index + 1
    ltrigs2 = tess._last_trig_index + 2

    @inbounds t1 = tess._trigs[i]
    old_t1_a = geta(t1)
    seta(t1, p)
    old_t1_b = t1._neighbour_b
    old_t1_c = t1._neighbour_c
    t1._neighbour_b = ltrigs1
    t1._neighbour_c = ltrigs2

    @inbounds t2 = tess._trigs[ltrigs1]
    setabc(t2, old_t1_a, p, getc(t1))
    t2._neighbour_a = i
    t2._neighbour_b = old_t1_b
    t2._neighbour_c = ltrigs2

    @inbounds t3 = tess._trigs[ltrigs2]
    setabc(t3, old_t1_a, getb(t1), p)
    t3._neighbour_a = i
    t3._neighbour_b = ltrigs1
    t3._neighbour_c = old_t1_c

    @inbounds nt2 = tess._trigs[t2._neighbour_b]
    if nt2._neighbour_a == i
        nt2._neighbour_a = ltrigs1
    elseif nt2._neighbour_b == i
        nt2._neighbour_b = ltrigs1
    else
        nt2._neighbour_c = ltrigs1
    end

    @inbounds nt3 = tess._trigs[t3._neighbour_c]
    if nt3._neighbour_a == i
        nt3._neighbour_a = ltrigs2
    elseif nt3._neighbour_b == i
        nt3._neighbour_b = ltrigs2
    else
        nt3._neighbour_c = ltrigs2
    end

    tess._last_trig_index += 2

    i
end

function _flipa!(tess::DelaunayTessellation2D, ix1::Int64, ix2::Int64)
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipaa!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipab!(tess, ix1, ix2, ot1, ot2)
    else
        _flipac!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipa!(tess::DelaunayTessellation2D,
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle, ot2::DelaunayTriangle)
    @inbounds n1 = tess._trigs[ot1._neighbour_a]
    if n1._neighbour_a == ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b == ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_c]
    if n2._neighbour_a == ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b == ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end

function _flipaa!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_b = getb(ot1)
    setb(ot1, geta(ot2))
    ot1._neighbour_a = ot2._neighbour_c
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    newc = geta(ot2)
    setabc(ot2, geta(ot1), old_ot1_geom_b, newc)
    ot2._neighbour_a = ot2._neighbour_b
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess, ix1, ix2, ot1, ot2)
end

function _flipab!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_b = getb(ot1)
    setb(ot1, getb(ot2))
    ot1._neighbour_a = ot2._neighbour_a
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    newc = getb(ot2)
    setabc(ot2, geta(ot1), old_ot1_geom_b, newc)
    ot2._neighbour_a = ot2._neighbour_c
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess, ix1, ix2, ot1, ot2)
end

function _flipac!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_b = getb(ot1)
    setb(ot1, getc(ot2))
    ot1._neighbour_a = ot2._neighbour_b
    old_ot1_neighbour_c = ot1._neighbour_c
    ot1._neighbour_c = ix2

    setab(ot2, geta(ot1), old_ot1_geom_b)
    ot2._neighbour_b = ix1
    ot2._neighbour_c = old_ot1_neighbour_c

    _endflipa!(tess, ix1, ix2, ot1, ot2)
end

########################

function _flipb!(tess::DelaunayTessellation2D{T},
                 ix1::Int64, ix2::Int64) where {T<:AbstractPoint2D}
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipba!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipbb!(tess, ix1, ix2, ot1, ot2)
    else
        _flipbc!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipb!(tess::DelaunayTessellation2D{T},
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    @inbounds n1 = tess._trigs[ot1._neighbour_b]
    if n1._neighbour_a == ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b == ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_a]
    if n2._neighbour_a == ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b == ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end

function _flipba!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_c = getc(ot1)
    setc(ot1, geta(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_c

    setbc(ot2, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_c = ix1

    _endflipb!(tess, ix1, ix2, ot1, ot2)
end

function _flipbb!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_c = getc(ot1)
    setc(ot1, getb(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_a

    newa = getb(ot2)
    setabc(ot2, newa, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_b = ot2._neighbour_c
    ot2._neighbour_c = ix1

    _endflipb!(tess, ix1, ix2, ot1, ot2)
end

function _flipbc!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_c = getc(ot1)
    setc(ot1, getc(ot2))
    old_ot1_neighbour_a = ot1._neighbour_a
    ot1._neighbour_a = ix2
    ot1._neighbour_b = ot2._neighbour_b

    newa = getc(ot2)
    setabc(ot2, newa, getb(ot1), old_ot1_geom_c)
    ot2._neighbour_b = ot2._neighbour_a
    ot2._neighbour_a = old_ot1_neighbour_a
    ot2._neighbour_c = ix1

    _endflipb!(tess, ix1, ix2, ot1, ot2)
end

########################

function _flipc!(tess::DelaunayTessellation2D{T},
                 ix1::Int64, ix2::Int64) where {T<:AbstractPoint2D}
    @inbounds ot1 = tess._trigs[ix1]
    @inbounds ot2 = tess._trigs[ix2]
    if ot2._neighbour_a == ix1
        _flipca!(tess, ix1, ix2, ot1, ot2)
    elseif ot2._neighbour_b == ix1
        _flipcb!(tess, ix1, ix2, ot1, ot2)
    else
        _flipcc!(tess, ix1, ix2, ot1, ot2)
    end
end

function _endflipc!(tess::DelaunayTessellation2D{T},
                    ix1::Int64, ix2::Int64,
                    ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    @inbounds n1 = tess._trigs[ot1._neighbour_c]
    if n1._neighbour_a == ix2
        n1._neighbour_a = ix1
    elseif n1._neighbour_b == ix2
        n1._neighbour_b = ix1
    else
        n1._neighbour_c = ix1
    end
    @inbounds n2 = tess._trigs[ot2._neighbour_b]
    if n2._neighbour_a == ix1
        n2._neighbour_a = ix2
    elseif n2._neighbour_b == ix1
        n2._neighbour_b = ix2
    else
        n2._neighbour_c = ix2
    end
end

function _flipca!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_a = geta(ot1)
    seta(ot1, geta(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_c

    newb = geta(ot2)
    setabc(ot2, old_ot1_geom_a, newb, getc(ot1))
    ot2._neighbour_a = ix1
    ot2._neighbour_c = ot2._neighbour_b
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess, ix1, ix2, ot1, ot2)
end

function _flipcb!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_a = geta(ot1)
    seta(ot1, getb(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_a

    setac(ot2, old_ot1_geom_a, getc(ot1))
    ot2._neighbour_a = ix1
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess, ix1, ix2, ot1, ot2)
end

function _flipcc!(tess::DelaunayTessellation2D{T},
                  ix1::Int64, ix2::Int64,
                  ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T}) where {T<:AbstractPoint2D}
    old_ot1_geom_a = geta(ot1)
    seta(ot1, getc(ot2))
    old_ot1_neighbour_b = ot1._neighbour_b
    ot1._neighbour_b = ix2
    ot1._neighbour_c = ot2._neighbour_b

    newb = getc(ot2)
    setabc(ot2, old_ot1_geom_a, newb, getc(ot1))
    ot2._neighbour_c = ot2._neighbour_a
    ot2._neighbour_a = ix1
    ot2._neighbour_b = old_ot1_neighbour_b

    _endflipc!(tess, ix1, ix2, ot1, ot2)
end

function _restoredelaunayhood!(tess::DelaunayTessellation2D{T},
                               ix_trig::Int64) where {T<:AbstractPoint2D}
    @inbounds center_pt = geta(tess._trigs[ix_trig])

    # `A` - edge
    push!(tess._edges_to_check, ix_trig)
    @inbounds while length(tess._edges_to_check) > 0
        trix = tess._edges_to_check[end]
        tr_i = tess._trigs[trix]
        nb_a = tr_i._neighbour_a
        if nb_a > 1
            tr_f = tess._trigs[nb_a]
            if incircle(tr_f, center_pt) > 0
                _flipa!(tess, trix, nb_a)
                push!(tess._edges_to_check, nb_a)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end

    # `B` - edge
    push!(tess._edges_to_check, tess._last_trig_index-1)
    @inbounds while length(tess._edges_to_check) > 0
        trix = tess._edges_to_check[end]
        tr_i = tess._trigs[trix]
        nb_b = tr_i._neighbour_b
        if nb_b > 1
            tr_f = tess._trigs[nb_b]
            if incircle(tr_f, center_pt) > 0
                _flipb!(tess, trix, nb_b)
                push!(tess._edges_to_check, nb_b)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end

    # `C` - edge
    push!(tess._edges_to_check, tess._last_trig_index)
    @inbounds while length(tess._edges_to_check) > 0
        trix = tess._edges_to_check[end]
        tr_i = tess._trigs[trix]
        nb_c = tr_i._neighbour_c
        if nb_c > 1
            tr_f = tess._trigs[nb_c]
            if incircle(tr_f, center_pt) > 0
                _flipc!(tess, trix, nb_c)
                push!(tess._edges_to_check, nb_c)
                continue
            end
        end
        pop!(tess._edges_to_check)
    end
end

# push a single point. Grows tessellation as needed
function push!(tess::DelaunayTessellation2D{T}, p::T) where {T<:AbstractPoint2D}
    tess._total_points_added += 1
    sizefit_at_least(tess, tess._total_points_added)
    i = _pushunfixed!(tess, p)
    _restoredelaunayhood!(tess, i)
end

# push an array in given order
function _pushunsorted!(tess::DelaunayTessellation2D{T}, a::Vector{T}) where {T<:AbstractPoint2D}
    sizehint!(tess, length(a))
    for p in a
        push!(tess, p)
    end
end

# push an array but sort it first for better performance
function push!(tess::DelaunayTessellation2D{T}, a::Vector{T}) where {T<:AbstractPoint2D}
    shuffle!(a)
    mssort!(a)
    _pushunsorted!(tess, a)
end

intensity(c::RGB)  = c.b
intensity(c::RGBA) = c.b
intensity(c) 	   = getfield(c, 1) # Workaround. Gray needs to be imported from images, which would take to long.

# Create DelaunayTessellation with npts points from an image
function from_image(img, npts)

    # placing points in places that represent the image
    pts = Point2D[]
    for i in 1:npts
        x = rand()
        y = rand()
        if intensity(img[Int64(floor(x * size(img)[1])) + 1, Int64(floor(y * size(img)[2])) + 1]) > 0.5
            if rand() < 0.100
                push!(pts, Point2D(1.0 + rand(), 1.0 + rand()))
            end
            continue
        end
        x /= 2
        y /= 2
        push!(pts, Point2D(1.0 + x + 0.3, 2.0 - y*2/3 - 0.3))
    end

    tess = DelaunayTessellation(npts)
    push!(tess, pts)

    tess
end


function getplotxy(edges)
    edges = collect(edges)
    x = Float64[]
    y = Float64[]
    for e in edges
        push!(x, getx(geta(e)))
        push!(x, getx(getb(e)))
        push!(x, NaN)
        push!(y, gety(geta(e)))
        push!(y, gety(getb(e)))
        push!(y, NaN)
    end
    (x, y)
end

end # module VoronoiDelaunay
