# importing all the stuff we need.
# You'll see many warnings because of using Images and Gadfly together
import VoronoiDelaunay: Point2D, DelaunayTessellation, push!, voronoiedges
import Images: imread
import Gadfly: set_default_plot_size, plot, Geom, Scale, cm

# placing points in places that represent the image...
img = imread("julia.png");
pts = Point2D[]
for i in 1:25000
    x = rand()
    y = rand()
    if img[int64(floor(x*size(img)[1]))+1,int64(floor(y*size(img)[2]))+1].c.b > 0.5
        if rand() < 0.100
            push!(pts, Point2D(1.0+rand(),1.0+rand()))
        end
        continue
    end
    x /= 2
    y /= 2
    push!(pts, Point2D(1.0+x+0.3,2.0-y*2/3-0.3))
end

# building the tesselation...
tess = DelaunayTessellation(length(pts))
push!(tess, pts)

# making the plot
set_default_plot_size(30cm,15cm)
x, y = getplotxy(voronoiedges(tess))
p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.25, maxvalue=1.75))

# save as SVG
draw(SVG("voroimage.svg", 8inch, 4inch), p)