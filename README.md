# ![icon](icon.png) Polygon Divider QGIS Plugin

### Produced by [Roy Ferguson Consultancy Ltd](https://royferguson.co.uk/) for [Zero Waste Scotland Ltd](http://www.zerowastescotland.org.uk/).

**Polygon Divider** is a plugin for [QGIS](http://www.qgis.org/en/site/) that takes a polygon and efficiently divides it into a number of 'squareish' polygons of a defined size, which is useful for a multitude of applications such as land parceling, environmental sampling, and so on.

As a simple worked example, you can take a polygon like this:

![hull](images/hull.png)

...and divide it into a number of smaller 'squareish' polygons of about 1000m<sup>2</sup> *(it would be exactly 1000m<sup>2</sup> if the polygon happened to have an area that precisely divides by 1000)*:

![divided hull](images/dividedhull.png)

Each output polygon inherits all of the attributes from its parent, as well as the following additional attributes:

* `ps_id`: a unique integer ID for each output polygon
* `ps_uuid`: a version 4 [uuid](https://en.wikipedia.org/wiki/Universally_unique_identifier)
* `ps_area`: the area of the polygon
* `ps_repPointX`: the X coordinate of a point gurranteed to be within the resulting polygon *(not necessarily the geometric centroid, as this is not gurranteed to be within the resulting polygon)*
* `ps_repPointY`: the X coordinate of a point gurranteed to be within the resulting polygon *(not necessarily the geometric centroid, as this is not gurranteed to be within the resulting polygon)*

The software should work well on some quite complex polygons:

![complex division example 1](images/complex1.png)

...even if they are very large:

![complex division example 1](images/complex2.png)

In cases where the algorithm finds the geometry difficult to divide, it will make the polygons slightly less square and more rectangular (as is illustrated in both of the abpve examples). If you find that it doesn't manage to divide a certain geometry at all, then you can normally remedy this by simplifying it a little.

Please do get in touch with a copy of any polygons that do not work, they will help us continue to improve this plugin!

A standalone version of this plugin *(i.e. not dependent upon QGIS)* will be available shortly.

#### Acknowledgements:
The production of this plugin was funded by [Zero Waste Scotland Ltd.](http://www.zerowastescotland.org.uk/). Development was greatly assisted by the accepted answer to [this](http://gis.stackexchange.com/questions/5300/dividing-polygon-into-specific-sizes-using-arcgis) forum post and by the [pyroots](https://pypi.python.org/pypi/pyroots/0.1.0) Python implementations of Brent's method.