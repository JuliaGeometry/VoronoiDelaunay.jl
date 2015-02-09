# importing all the stuff we need.
# You'll see many warnings because of using Images and Gadfly together
import VoronoiDelaunay: from_image, voronoiedges, getplotxy
import Images: imread
import Gadfly: set_default_plot_size, plot, Geom, Scale, cm, draw, SVG, inch

img = imread("julia.png")
# placing points in places that represent the image.
tess = from_image(img, 25000)

# making the plot
set_default_plot_size(30cm,15cm)
x, y = getplotxy(voronoiedges(tess))
p = plot(x=x, y=y, Geom.path, Scale.x_continuous(minvalue=1.0, maxvalue=2.0), Scale.y_continuous(minvalue=1.25, maxvalue=1.75))

# save as SVG
draw(SVG("voroimage.svg", 8inch, 4inch), p)
