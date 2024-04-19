# ![icon](icon.png) Polygon Divider QGIS Plugin

### Now available for QGIS 3 (QGIS2 version available in [this branch](https://github.com/jonnyhuck/RFCL-PolygonDivider/tree/QGIS2-version))

**Polygon Divider** is a plugin for [QGIS](http://www.qgis.org/en/site/) that takes a polygon and efficiently divides it into a number of 'squareish' polygons of a defined size, which is useful for a multitude of applications such as land parceling, environmental sampling, and so on.

As a simple worked example, you can take a polygon like this:

![hull](images/hull.png)

...and divide it into a number of smaller 'squareish' polygons of about 1000m<sup>2</sup> *(it would be exactly 1000m<sup>2</sup> if the polygon happened to have an area that precisely divides by 1000)*:

![divided hull](images/dividedhull.png)

There are two options available to the user when dividing a polygon. The above uses the **absorb** method, whereby all of the polygons are slightly larger than the requested size, in order to 'absorb' any odd-sized 'offcuts' that would otherwise be left behind. The other alternative would be the **offcut** method. An example of this is given below, in which all of the polygons would be the precise size as requested, except for the light green one at the very top, which represente the *'offcut'*:

![offcut hull](images/hulloffcut.png)

The choice of which is best from the above will vary depending upon the specific requirements of the user, and the simplicity of the polygon in question. 

Each of the above cutting methods can be undertaken in 4 directions: **left-right**, **right-left**, **bottom-top** and **top-bottom**. Again, depending upon the shape of the specific polygon to be divided, better results might be achieved in some directions than others.

Each output polygon inherits all of the attributes from its parent, as well as the following additional attributes:

* `ps_id`: a unique integer ID for each output polygon
* `ps_uuid`: a version 4 [uuid](https://en.wikipedia.org/wiki/Universally_unique_identifier)
* `ps_area`: the area of the polygon
* `ps_repPointX`: the X coordinate of a point guaranteed to be within the resulting polygon *(not necessarily the geometric centroid, as this is not gurranteed to be within the resulting polygon)*
* `ps_repPointY`: the Y coordinate of a point guaranteed to be within the resulting polygon *(not necessarily the geometric centroid, as this is not gurranteed to be within the resulting polygon)*

The software should work well on some quite complex polygons:

![complex division example 1](images/complex1.png)

...even if they are very large:

![complex division example 1](images/complex2.png)

In cases where the algorithm finds the geometry difficult to divide, it will make the polygons slightly less square and more rectangular (as is illustrated in both of the above examples). If you find that it doesn't manage to divide a certain geometry at all, then you can normally remedy this by simplifying it a little.

### Rotation

![rotation example original](images/rotation_example_original.png)

consider the following polygon:
if we divide it without rotating it prior, we get the following result:

![rotation example divided without rotation](images/rotation_example_divided_not_rotated.png)

if we rotate it 45 degrees, we get the following result:

![rotation example divided with rotation](images/rotation_example_divided_rotated.png)

Please do get in touch with a copy of any polygons that do not work, they will help us continue to improve this plugin!

#### Data Considerations:
The Polygon Divider expects planar geometry (x,y coordinates) and will not divide geographical coordinates (degrees).  If you are using geographical coordinates you must save your dataset using a Projected Coordinate System (It is recommended that you use either a local or equal area projection).  Note that while QGIS will allow you to change the CRS for the project/layer, doing so is not enough as it does not save the data and only reprojects the view for the user.

#### Acknowledgements:
The QGIS2 (original) version of this plugin was funded by [Zero Waste Scotland Ltd.](http://www.zerowastescotland.org.uk/). The conversion to QGIS3 was funded by [Deutsche Forestservice GMBH](https://www.dfs-online.de/).

Development was greatly assisted by the accepted answer to [this](http://gis.stackexchange.com/questions/5300/dividing-polygon-into-specific-sizes-using-arcgis) forum post and the [pyroots](https://pypi.python.org/pypi/pyroots/0.1.0)  implementation of Brent's method.
