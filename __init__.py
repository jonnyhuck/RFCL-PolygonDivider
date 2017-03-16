# -*- coding: utf-8 -*-
"""
/***************************************************************************
 PolygonDivider
                                 A QGIS plugin
 Divides Polygons
                             -------------------
        begin                : 2017-01-22
        copyright            : (C) 2017 by Roy Ferguson Consulting
        email                : jonnyhuck@gmail.com
        git sha              : $Format:%H$
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/
 This script initializes the plugin, making it known to QGIS.
"""


# noinspection PyPep8Naming
def classFactory(iface):  # pylint: disable=invalid-name
    """Load PolygonDivider class from file PolygonDivider.

    :param iface: A QGIS interface instance.
    :type iface: QgsInterface
    """
    #
    from .polygondivider import PolygonDivider
    return PolygonDivider(iface)
