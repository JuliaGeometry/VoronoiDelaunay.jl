using VoronoiDelaunay, Gadfly

points = [
 Point2D(-1.0563841812533212, -1.4606363138997696) 
 Point2D(0.06312982975128989, -0.48031801152366027)
 Point2D(0.1624918689993189, -0.19919450833195906) 
 Point2D(-1.5293344878962758, -0.7657808444340142) 
 Point2D(0.5319064220493406, 0.6107808132476504)   
 Point2D(-0.3670342825169435, 0.8207427582546951)  
 Point2D(-1.9797019290444364, -0.5066353099040788) 
 Point2D(-1.5811584920674324, 1.0346409888830976)  
 Point2D(1.2185165319349451, 1.4177209167374605)   
 Point2D(-1.5991536318191626, -1.3063986775765466) ];

tess, ranges = DelaunayTessellation( points )
xy = getplotxy( delaunayedges( tess, ranges ) )
l1 = layer( x=getx.(points), y=gety.(points), Theme(default_color=colorant"orange"), Geom.point )
l2 = layer( x=xy[1], y=xy[2], Geom.path )
plot( l1, l2, Coord.cartesian(fixed=true) )