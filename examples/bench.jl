using VoronoiDelaunay
println(workers())

function test(tess, points)
	push!(tess, points)
end
function test()
	numPoints = 0
	width 	  = max_coord - min_coord
	time_map  = Dict{Int, Float64}() 
	numPoints = 10^6
	for i=1:10
		vInputPoints = Point2D[Point(min_coord+rand()*width, min_coord+rand()*width) for i in 1:numPoints]
		tess = DelaunayTessellation(numPoints)
		@time test(tess, vInputPoints)
	end	

end

@fastmath test()

#=
elapsed time: 2.140358087 seconds (18 MB allocated)
elapsed time: 1.252718545 seconds (513 kB allocated)
elapsed time: 1.186685373 seconds (456 kB allocated)
elapsed time: 1.268425153 seconds (412 kB allocated)
elapsed time: 1.210739505 seconds (447 kB allocated)
elapsed time: 1.222889958 seconds (447 kB allocated)
elapsed time: 1.200753356 seconds (504 kB allocated)
elapsed time: 1.21738593 seconds (552 kB allocated)
elapsed time: 1.217354692 seconds (465 kB allocated)
elapsed time: 1.224980633 seconds (443 kB allocated)
=#
#=
elapsed time: 1.815125503 seconds (18 MB allocated)
elapsed time: 1.310488535 seconds (434 kB allocated)
elapsed time: 1.265986406 seconds (556 kB allocated)
elapsed time: 1.258479167 seconds (512 kB allocated)
elapsed time: 1.263162099 seconds (373 kB allocated)
elapsed time: 1.275858315 seconds (469 kB allocated)
elapsed time: 1.218306987 seconds (439 kB allocated)
elapsed time: 1.192661087 seconds (499 kB allocated)
elapsed time: 1.24761614 seconds (456 kB allocated)
elapsed time: 1.263454391 seconds (469 kB allocated)
=#