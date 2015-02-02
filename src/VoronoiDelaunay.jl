module VoronoiDelaunay

# Fast, robust 2D Voronoi/Delaunay tesselation
# Implementation follows algorithms described in http://arxiv.org/abs/0901.4107
# and used (for e.g.) in the Illustris Simulation
# http://www.illustris-project.org/
#
# Author: Ariel Keselman (skariel@gmail.com)
# License: MIT
# Bug reportss welcome!

export
	DelaunayTessellation, DelaunayTessellation2D, sizehint, isexternal,
	min_coord, max_coord, locate, movea, moveb, movec,
	delaunayedges, voronoiedges, voronoiedgeswithoutgenerators,
	start, next, done,
	findindex, push!,
	Point, Point2D, AbstractPoint2D, getx, gety, geta, getb, getc,
	getplotxy


using GeometricalPredicates
import GeometricalPredicates
import GeometricalPredicates: geta, getb, getc

import Base:push!, start, done, next, sizehint, copy

const min_coord = GeometricalPredicates.min_coord + eps(Float64)
const max_coord = GeometricalPredicates.max_coord - eps(Float64)

type DelaunayTriangle{T<:AbstractPoint2D} <: AbstractNegativelyOrientedTriangle
    _a::T; _b::T; _c::T
    _bx::Float64; _by::Float64
    _cx::Float64; _cy::Float64
    _px::Float64; _py::Float64
    _pr2::Float64
	_neighbour_a::Int64
	_neighbour_b::Int64
	_neighbour_c::Int64
    function DelaunayTriangle(pa::T, pb::T, pc::T, na::Int64, nb::Int64, nc::Int64)
        t = new(pa, pb, pc, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, na, nb, nc)
        clean!(t)
        t
    end
    # this constructor is good for making copies
    DelaunayTriangle(pa::T, pb::T, pc::T, bx::Float64, by::Float64, cx::Float64, cy::Float64, px::Float64, py::Float64, pr2::Float64, na::Int64, nb::Int64, nc::Int64) =
        new(pa, pb, pc, bx, by, cx, cy, px, py, pr2, na, nb, nc)
end
DelaunayTriangle{T<:AbstractPoint2D}(pa::T, pb::T, pc::T, na::Int64, nb::Int64, nc::Int64) = DelaunayTriangle{T}(pa, pb, pc, na, nb, nc)
DelaunayTriangle{T<:AbstractPoint2D}(pa::T, pb::T, pc::T, bx::Float64, by::Float64, cx::Float64, cy::Float64, px::Float64, py::Float64, pr2::Float64, na::Int64, nb::Int64, nc::Int64) =
    DelaunayTriangle{T}(pa, pb, pc, bx, by, cx, cy, px, py, pr2, na, nb, nc)

copy{T<:AbstractPoint2D}(t::DelaunayTriangle{T}) =
	DelaunayTriangle(
		t._a, t._b, t._c,
		t._bx, t._by,
		t._cx, t._cy,
		t._px, t._py,
		t._pr2,
		t._neighbour_a, t._neighbour_b, t._neighbour_c
	)

isexternal{T<:AbstractPoint2D}(t::DelaunayTriangle{T}) = 
	getx(geta(t)) < min_coord || getx(geta(t)) > max_coord ||
	getx(getb(t)) < min_coord || getx(getb(t)) > max_coord ||
	getx(getc(t)) < min_coord || getx(getc(t)) > max_coord


type DelaunayTessellation2D{T<:AbstractPoint2D}
	_trigs::Array{DelaunayTriangle{T}, 1}
	_last_trig_index::Int64
	_edges_to_check::Array{Int64, 1}
	_total_points_added::Int64
	function DelaunayTessellation2D(n::Int64 = 100)		
		const a = T(GeometricalPredicates.min_coord, GeometricalPredicates.min_coord)
		const b = T(GeometricalPredicates.min_coord, GeometricalPredicates.max_coord)
		const c = T(GeometricalPredicates.max_coord, GeometricalPredicates.min_coord)
		const d = T(GeometricalPredicates.max_coord, GeometricalPredicates.max_coord)
		const _trigs = DelaunayTriangle{T}[DelaunayTriangle{T}(d,c,b, 2,1,1), DelaunayTriangle{T}(a,b,c, 3,1,1), DelaunayTriangle{T}(d,c,b, 2,1,1)]
		t = new(_trigs, 3, Int64[], 0)
		sizehint(t._edges_to_check, 1000)
		sizehint(t, n)
	end
end
DelaunayTessellation2D(n::Int64) = DelaunayTessellation2D{Point2D}(n)
DelaunayTessellation2D{T<:AbstractPoint2D}(n::Int64, ::T) = DelaunayTessellation2D{T}(n)
DelaunayTessellation(n::Int64=100) = DelaunayTessellation2D(n)

function sizehint{T<:AbstractPoint2D}(t::DelaunayTessellation2D{T}, n::Int64)
	const required_total_size = 2*n+10
	required_total_size <= length(t._trigs) && return
	sizehint(t._trigs, required_total_size)
	while length(t._trigs) < required_total_size
		push!(t._trigs, copy(t._trigs[end]))
	end	 
	t
end

# growing strategy
function sizefit_at_least{T<:AbstractPoint2D}(t::DelaunayTessellation2D{T}, n::Int64)
	const minimal_acceptable_actual_size = 2*n+10
	minimal_acceptable_actual_size <= length(t._trigs) && return
	required_total_size = length(t._trigs)
	while required_total_size < minimal_acceptable_actual_size
		required_total_size += required_total_size >>> 1
	end
	sizehint(t._trigs, required_total_size)
	while length(t._trigs) < required_total_size
		push!(t._trigs, copy(t._trigs[end]))
	end	 
	t
end

immutable DelaunayEdge{T<:AbstractPoint2D}
	_a::T
	_b::T
end
geta{T<:AbstractPoint2D}(e::DelaunayEdge{T}) = e._a
getb{T<:AbstractPoint2D}(e::DelaunayEdge{T}) = e._b

immutable VoronoiEdge{T<:AbstractPoint2D}
	_a::Point2D
	_b::Point2D
	_generator_a::T
	_generator_b::T
end
geta{T<:AbstractPoint2D}(e::VoronoiEdge{T}) = e._a
getb{T<:AbstractPoint2D}(e::VoronoiEdge{T}) = e._b
getgena{T<:AbstractPoint2D}(e::VoronoiEdge{T}) = e._generator_a
getgenb{T<:AbstractPoint2D}(e::VoronoiEdge{T}) = e._generator_b

immutable VoronoiEdgeWithoutGenerators
	_a::Point2D
	_b::Point2D
end
geta(e::VoronoiEdgeWithoutGenerators) = e._a
getb(e::VoronoiEdgeWithoutGenerators) = e._b

# TODO: is an iterator faster?
function delaunayedges(t::DelaunayTessellation2D)
	visited = zeros(Bool, t._last_trig_index)
	function delaunayiterator()
		@inbounds for ix in 2:t._last_trig_index
			const tr = t._trigs[ix]
			isexternal(tr) && continue
			visited[ix] && continue
			visited[ix] = true
			const ix_na = tr._neighbour_a
			if !visited[ix_na]
				produce(DelaunayEdge(getb(tr), getc(tr)))
			end
			const ix_nb = tr._neighbour_b
			if !visited[ix_nb] 
				produce(DelaunayEdge(geta(tr), getc(tr)))
			end
			const ix_nc = tr._neighbour_c
			if !visited[ix_nc] 
				produce(DelaunayEdge(geta(tr), getb(tr)))
			end
		end
	end
	Task(delaunayiterator)
end

# TODO: is an iterator faster?
function voronoiedges(t::DelaunayTessellation2D)
	visited = zeros(Bool, t._last_trig_index)
	visited[1] = true
	function voronoiiterator()
		@inbounds for ix in 2:t._last_trig_index
			visited[ix] && continue
			const tr = t._trigs[ix]
			visited[ix] = true
			isexternal(tr) && continue
			const cc = circumcenter(tr)

			const ix_na = tr._neighbour_a
			if !visited[ix_na] && !isexternal(t._trigs[ix_na])
				const nb = t._trigs[ix_na]
				produce(VoronoiEdge(cc, circumcenter(nb), geta(tr), nb._neighbour_a==ix? geta(nb) : nb._neighbour_b==ix? getb(nb): getc(nb)))
			end
			const ix_nb = tr._neighbour_b
			if !visited[ix_nb] && !isexternal(t._trigs[ix_nb])
				const nb = t._trigs[ix_nb]
				produce(VoronoiEdge(cc, circumcenter(nb), getb(tr), nb._neighbour_a==ix? geta(nb) : nb._neighbour_b==ix? getb(nb): getc(nb)))
			end
			const ix_nc = tr._neighbour_c
			if !visited[ix_nc] && !isexternal(t._trigs[ix_nb])
				const nb = t._trigs[ix_nc]
				produce(VoronoiEdge(cc, circumcenter(nb), getc(tr), nb._neighbour_a==ix? geta(nb) : nb._neighbour_b==ix? getb(nb): getc(nb)))
			end
		end
	end
	Task(voronoiiterator)
end

# TODO: is an iterator faster?
function voronoiedgeswithoutgenerators(t::DelaunayTessellation2D)
	visited = zeros(Bool, t._last_trig_index)
	visited[1] = true
	function voronoiiterator()
		@inbounds for ix in 2:t._last_trig_index
			visited[ix] && continue
			const tr = t._trigs[ix]
			visited[ix] = true
			isexternal(tr) && continue
			const cc = circumcenter(tr)

			const ix_na = tr._neighbour_a
			if !visited[ix_na] && !isexternal(t._trigs[ix_na])
				const nb = t._trigs[ix_na]
				produce(VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
			end
			const ix_nb = tr._neighbour_b
			if !visited[ix_nb] && !isexternal(t._trigs[ix_nb])
				const nb = t._trigs[ix_nb]
				produce(VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
			end
			const ix_nc = tr._neighbour_c
			if !visited[ix_nc] && !isexternal(t._trigs[ix_nc])
				const nb = t._trigs[ix_nc]
				produce(VoronoiEdgeWithoutGenerators(cc, circumcenter(nb)))
			end
		end
	end
	Task(voronoiiterator)
end

type TrigIter
	ix::Int64
end
start(t::DelaunayTessellation2D) = TrigIter(2)
function done(t::DelaunayTessellation2D, it::TrigIter)
	while isexternal(t._trigs[it.ix]) && it.ix <= t._last_trig_index
		it.ix += 1
	end	
	it.ix > t._last_trig_index
end
function next(t::DelaunayTessellation2D, it::TrigIter)
	@inbounds const trig = t._trigs[it.ix]
	it.ix += 1
	(trig, it)
end

function findindex{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, p::T)
	i::Int64 = tess._last_trig_index
	while true		
		@inbounds const w = intriangle(tess._trigs[i], p)
		w > 0 && return i
		@inbounds const tr = tess._trigs[i]
		i = w==-1? tr._neighbour_a : w==-2? tr._neighbour_b : tr._neighbour_c
	end
end

locate{T<:AbstractPoint2D}(t::DelaunayTessellation2D{T}, p::T) = t._trigs[findindex(t, p)]

movea{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, trig::DelaunayTriangle{T}) = tess._trigs[trig._neighbour_a]
moveb{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, trig::DelaunayTriangle{T}) = tess._trigs[trig._neighbour_b]
movec{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, trig::DelaunayTriangle{T}) = tess._trigs[trig._neighbour_c]

function _pushunfixed!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, p::T)
	i = findindex(tess, p)
	const ltrigs1::Int64 = tess._last_trig_index+1
	const ltrigs2::Int64 = tess._last_trig_index+2

	@inbounds const t1 = tess._trigs[i]
	old_t1_a = geta(t1)
	seta(t1, p)
	old_t1_b = t1._neighbour_b
	old_t1_c = t1._neighbour_c
	t1._neighbour_b = ltrigs1
	t1._neighbour_c = ltrigs2

	@inbounds const t2 = tess._trigs[ltrigs1]
	setabc(t2, old_t1_a, p, getc(t1))
	t2._neighbour_a = i
	t2._neighbour_b = old_t1_b
	t2._neighbour_c = ltrigs2

	@inbounds const t3 = tess._trigs[ltrigs2]
	setabc(t3, old_t1_a, getb(t1), p)
	t3._neighbour_a = i
	t3._neighbour_b = ltrigs1
	t3._neighbour_c = old_t1_c

	@inbounds const nt2 = tess._trigs[t2._neighbour_b]
	if nt2._neighbour_a==i
		nt2._neighbour_a = ltrigs1
	elseif nt2._neighbour_b==i
		nt2._neighbour_b = ltrigs1
	else
		nt2._neighbour_c = ltrigs1
	end

	@inbounds const nt3 = tess._trigs[t3._neighbour_c]
	if nt3._neighbour_a==i
		nt3._neighbour_a = ltrigs2
	elseif nt3._neighbour_b==i
		nt3._neighbour_b = ltrigs2
	else
		nt3._neighbour_c = ltrigs2
	end

	tess._last_trig_index += 2

	i
end

function _flipa!(tess::DelaunayTessellation2D, ix1::Int64, ix2::Int64)
	@inbounds const ot1 = tess._trigs[ix1]
	@inbounds const ot2 = tess._trigs[ix2]
	if ot2._neighbour_a == ix1
		_flipaa!(tess, ix1, ix2, ot1, ot2)
	elseif ot2._neighbour_b == ix1
		_flipab!(tess, ix1, ix2, ot1, ot2)
	else
		_flipac!(tess, ix1, ix2, ot1, ot2)
	end
end

function _endflipa!(tess::DelaunayTessellation2D, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle, ot2::DelaunayTriangle)
	@inbounds const n1 = tess._trigs[ot1._neighbour_a]
	if n1._neighbour_a==ix2
		n1._neighbour_a = ix1
	elseif n1._neighbour_b==ix2
		n1._neighbour_b = ix1
	else
		n1._neighbour_c = ix1
	end
	@inbounds const n2 = tess._trigs[ot2._neighbour_c]
	if n2._neighbour_a==ix1
		n2._neighbour_a = ix2
	elseif n2._neighbour_b==ix1
		n2._neighbour_b = ix2
	else
		n2._neighbour_c = ix2
	end
end

function _flipaa!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_b = getb(ot1)
	setb(ot1, geta(ot2))
	ot1._neighbour_a = ot2._neighbour_c
	old_ot1_neighbour_c = ot1._neighbour_c
	ot1._neighbour_c = ix2

	setabc(ot2, geta(ot1), old_ot1_geom_b, geta(ot2))
	ot2._neighbour_a = ot2._neighbour_b
	ot2._neighbour_b = ix1
	ot2._neighbour_c = old_ot1_neighbour_c	

	_endflipa!(tess,ix1,ix2,ot1,ot2)
end

function _flipab!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_b = getb(ot1)

	setb(ot1, getb(ot2))
	ot1._neighbour_a = ot2._neighbour_a
	old_ot1_neighbour_c = ot1._neighbour_c
	ot1._neighbour_c = ix2


	setabc(ot2, geta(ot1), old_ot1_geom_b, getb(ot2))
	ot2._neighbour_a = ot2._neighbour_c
	ot2._neighbour_b = ix1
	ot2._neighbour_c = old_ot1_neighbour_c

	_endflipa!(tess,ix1,ix2,ot1,ot2)
end

function _flipac!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_b = getb(ot1)
	setb(ot1, getc(ot2))
	ot1._neighbour_a = ot2._neighbour_b
	old_ot1_neighbour_c = ot1._neighbour_c
	ot1._neighbour_c = ix2

	setab(ot2, geta(ot1), old_ot1_geom_b)
	ot2._neighbour_b = ix1
	ot2._neighbour_c = old_ot1_neighbour_c

	_endflipa!(tess,ix1,ix2,ot1,ot2)
end

########################

function _flipb!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64)
	@inbounds const ot1 = tess._trigs[ix1]
	@inbounds const ot2 = tess._trigs[ix2]
	if ot2._neighbour_a == ix1
		_flipba!(tess, ix1, ix2, ot1, ot2)
	elseif ot2._neighbour_b == ix1
		_flipbb!(tess, ix1, ix2, ot1, ot2)
	else
		_flipbc!(tess, ix1, ix2, ot1, ot2)
	end
end

function _endflipb!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	@inbounds const n1 = tess._trigs[ot1._neighbour_b]
	if n1._neighbour_a==ix2
		n1._neighbour_a = ix1
	elseif n1._neighbour_b==ix2
		n1._neighbour_b = ix1
	else
		n1._neighbour_c = ix1
	end
	@inbounds const n2 = tess._trigs[ot2._neighbour_a]
	if n2._neighbour_a==ix1
		n2._neighbour_a = ix2
	elseif n2._neighbour_b==ix1
		n2._neighbour_b = ix2
	else
		n2._neighbour_c = ix2
	end
end

function _flipba!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_c = getc(ot1)
	setc(ot1, geta(ot2))
	old_ot1_neighbour_a = ot1._neighbour_a
	ot1._neighbour_a = ix2
	ot1._neighbour_b = ot2._neighbour_c

	setbc(ot2, getb(ot1), old_ot1_geom_c)
	ot2._neighbour_a = old_ot1_neighbour_a
	ot2._neighbour_c = ix1

	_endflipb!(tess,ix1,ix2,ot1,ot2)
end

function _flipbb!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_c = getc(ot1)
	setc(ot1, getb(ot2))
	old_ot1_neighbour_a = ot1._neighbour_a
	ot1._neighbour_a = ix2
	ot1._neighbour_b = ot2._neighbour_a

	setabc(ot2, getb(ot2), getb(ot1), old_ot1_geom_c)
	ot2._neighbour_a = old_ot1_neighbour_a
	ot2._neighbour_b = ot2._neighbour_c
	ot2._neighbour_c = ix1

	_endflipb!(tess,ix1,ix2,ot1,ot2)
end

function _flipbc!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_c = getc(ot1)
	setc(ot1, getc(ot2))
	old_ot1_neighbour_a = ot1._neighbour_a
	ot1._neighbour_a = ix2
	ot1._neighbour_b = ot2._neighbour_b

	setabc(ot2, getc(ot2), getb(ot1), old_ot1_geom_c)
	ot2._neighbour_b = ot2._neighbour_a
	ot2._neighbour_a = old_ot1_neighbour_a
	ot2._neighbour_c = ix1

	_endflipb!(tess,ix1,ix2,ot1,ot2)
end

########################

function _flipc!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64)
	@inbounds const ot1 = tess._trigs[ix1]
	@inbounds const ot2 = tess._trigs[ix2]
	if ot2._neighbour_a == ix1
		_flipca!(tess, ix1, ix2, ot1, ot2)
	elseif ot2._neighbour_b == ix1
		_flipcb!(tess, ix1, ix2, ot1, ot2)
	else
		_flipcc!(tess, ix1, ix2, ot1, ot2)
	end
end

function _endflipc!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	@inbounds const n1 = tess._trigs[ot1._neighbour_c]
	if n1._neighbour_a==ix2
		n1._neighbour_a = ix1
	elseif n1._neighbour_b==ix2
		n1._neighbour_b = ix1
	else
		n1._neighbour_c = ix1
	end
	@inbounds const n2 = tess._trigs[ot2._neighbour_b]
	if n2._neighbour_a==ix1
		n2._neighbour_a = ix2
	elseif n2._neighbour_b==ix1
		n2._neighbour_b = ix2
	else
		n2._neighbour_c = ix2
	end
end

function _flipca!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_a = geta(ot1)
	seta(ot1, geta(ot2))
	old_ot1_neighbour_b = ot1._neighbour_b
	ot1._neighbour_b = ix2
	ot1._neighbour_c = ot2._neighbour_c

	setabc(ot2, old_ot1_geom_a, geta(ot2), getc(ot1))
	ot2._neighbour_a = ix1
	ot2._neighbour_c = ot2._neighbour_b
	ot2._neighbour_b = old_ot1_neighbour_b

	_endflipc!(tess,ix1,ix2,ot1,ot2)
end

function _flipcb!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_a = geta(ot1)
	seta(ot1, getb(ot2))
	old_ot1_neighbour_b = ot1._neighbour_b
	ot1._neighbour_b = ix2
	ot1._neighbour_c = ot2._neighbour_a

	setac(ot2, old_ot1_geom_a, getc(ot1))
	ot2._neighbour_a = ix1
	ot2._neighbour_b = old_ot1_neighbour_b

	_endflipc!(tess,ix1,ix2,ot1,ot2)
end

function _flipcc!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix1::Int64, ix2::Int64, ot1::DelaunayTriangle{T}, ot2::DelaunayTriangle{T})
	old_ot1_geom_a = geta(ot1)
	seta(ot1, getc(ot2))
	old_ot1_neighbour_b = ot1._neighbour_b
	ot1._neighbour_b = ix2
	ot1._neighbour_c = ot2._neighbour_b

	setabc(ot2, old_ot1_geom_a, getc(ot2), getc(ot1))
	ot2._neighbour_c = ot2._neighbour_a
	ot2._neighbour_a = ix1
	ot2._neighbour_b = old_ot1_neighbour_b

	_endflipc!(tess,ix1,ix2,ot1,ot2)
end

function _restoredelaunayhood!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, ix_trig::Int64)
	@inbounds const center_pt = geta(tess._trigs[ix_trig])

	# `A` - edge
	push!(tess._edges_to_check, ix_trig)
    while length(tess._edges_to_check) > 0
		@inbounds const trix = tess._edges_to_check[end]
		@inbounds const tr_i = tess._trigs[trix]
		@inbounds const nb_a = tr_i._neighbour_a
		if nb_a > 1
			@inbounds const tr_f = tess._trigs[nb_a]
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
    while length(tess._edges_to_check) > 0
		@inbounds const trix = tess._edges_to_check[end]
		@inbounds const tr_i = tess._trigs[trix]
		@inbounds const nb_b = tr_i._neighbour_b
		if nb_b > 1
			@inbounds const tr_f = tess._trigs[nb_b]
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
    while length(tess._edges_to_check) > 0
		@inbounds const trix = tess._edges_to_check[end]
		@inbounds const tr_i = tess._trigs[trix]
		@inbounds const nb_c = tr_i._neighbour_c
		if nb_c > 1
			@inbounds const tr_f = tess._trigs[nb_c]
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
function push!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, p::T)
	tess._total_points_added += 1
	sizefit_at_least(tess, tess._total_points_added)
	const i = _pushunfixed!(tess, p)
	_restoredelaunayhood!(tess, i)
end

# push an array in given order
function _pushunsorted!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, a::Array{T, 1})
	sizehint(tess, length(a))
	for p in a
		push!(tess, p)
	end
end

# push an array but sort it first for better performance
function push!{T<:AbstractPoint2D}(tess::DelaunayTessellation2D{T}, a::Array{T, 1})
	shuffle!(a)
	mssort!(a)
	_pushunsorted!(tess, a)
end

# Create DelaunayTessellation with npts points from an image
function from_image(img, npts)

	# placing points in places that represent the image
	pts = Point2D[]
	for i in 1:npts
		x = rand()
		y = rand()
		if img[int64(floor(x * size(img)[1])) + 1, int64(floor(y * size(img)[2])) + 1].c.b > 0.5
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
