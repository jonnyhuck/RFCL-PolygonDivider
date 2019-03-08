# -*- coding: utf-8 -*-

"""
	/***************************************************************************
	PolygonDivider
								 A QGIS plugin
	Divides Polygons
							-------------------
		begin				 : 2017-01-22
		git sha				 : $Format:%H$
		copyright			 : (C) 2017 by Roy Ferguson Consultancy
		email				 : jonnyhuck@gmail.com
	***************************************************************************/

	/***************************************************************************
	*																		 *
	*	 This program is free software; you can redistribute it and/or modify  *
	*	 it under the terms of the GNU General Public License as published by  *
	*	 the Free Software Foundation; either version 2 of the License, or	 *
	*	 (at your option) any later version.								 *
	*																		 *
	***************************************************************************/

	* This script divides a polygon into squareish sections of a specified size
	* This version is threaded as per: https://snorfalorpagus.net/blog/2013/12/07/multithreading-in-qgis-python-plugins/
	*
	* ERROR WE CAN FIX BY ADJUSTING TOLERANCE / N_SUBDIVISIONS:
	*  - Bracket is smaller than tolerance: the shape got smaller? OR got dramatically bigger. Can we check this? This is where we just want to cut at the last location that worked and re-calculate the division stuff.
	*
	* TODO'S:
	*  - Where / how often should we calculate the desired area?
	*  - How should we be dealing with reversing direction for subdivision? Undoing changes seems to make it worse...
	*  - Need to re-think how to deal with problems in subdivision - perhaps need to calculate all subdivisions then save all at once so we can roll back?
	*  - Look at saving last good bounds to narrow search interval after adjusting tolerance? Maybe use bisection to minimise the adjustment in tolerance?
	*  - More minor TODO's throughout the code...
	*
	* @author jonnyhuck
	*
"""

from PyQt4 import QtCore
from PyQt4.QtCore import QVariant, QObject, pyqtSignal, QSettings, QTranslator, qVersion, QCoreApplication, Qt, QThread
from PyQt4.QtGui import QAction, QIcon, QFileDialog, QProgressBar, QPushButton
import resources, os.path, traceback
import qgis.utils, sys
from uuid import uuid4
from math import sqrt
from polygondivider_dialog import PolygonDividerDialog
import psycopg2
from qgis.core import *		#TODO: Should be more specific here
# from qgis.core import QgsVectorLayer, QgsMapLayerRegistry, QgsMessageLog
from qgis.gui import QgsMessageBar


'''

***************************** STUFF FOR THE WORKER THREAD ********************************

'''


class AbstractWorker(QtCore.QObject):
	"""Abstract worker, inherit from this and implement the work method"""

	# available signals to be used in the concrete worker
	finished = QtCore.pyqtSignal(object)
	error = QtCore.pyqtSignal(Exception, basestring)
	progress = QtCore.pyqtSignal(float)
	toggle_show_progress = QtCore.pyqtSignal(bool)
	set_message = QtCore.pyqtSignal(str)
	toggle_show_cancel = QtCore.pyqtSignal(bool)
	
	# private signal, don't use in concrete workers this is automatically emitted if the result is not None
	successfully_finished = QtCore.pyqtSignal(object)

	def __init__(self):
		QtCore.QObject.__init__(self)
		self.killed = False

	def run(self):
		try:
			result = self.work()
			self.finished.emit(result)
		except UserAbortedNotification:
			self.finished.emit(None)
		except Exception, e:
			# forward the exception upstream
			self.error.emit(e, traceback.format_exc())
			self.finished.emit(None)

	# this is overridden below
	def work(self):
		""" Reimplement this putting your calculation here
			available are:
				self.progress.emit(0-100)
				self.killed
			:returns a python object - use None if killed is true
		"""
		raise NotImplementedError

	def kill(self):
		self.killed = True
		self.set_message.emit('Aborting...')
		self.toggle_show_progress.emit(False)


class BrentError(Exception):
	"""
	* Simple class for exceptions from Brent's Method.
	"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


class CoreWorker(AbstractWorker):
	"""
	* Core worker thread - created by Exegesis SDM to process in batches
	"""

	# add additional signal to say batch finished and trigger new batch
	batch_finished = QtCore.pyqtSignal(bool)
	
	def __init__(self, iface, inLayer, outputType, outFilePath, pgDetails, chunkSize, targetArea, absorbFlag, direction):
		"""
		* Initialise Thread
		"""

		# superclass
		AbstractWorker.__init__(self)

		# import args to thread
		self.iface = iface
		self.layer = inLayer
		self.outputType = outputType
		self.outFilePath = outFilePath
		self.pgDetails = pgDetails
		self.chunk_size = chunkSize
		self.target_area = targetArea
		self.absorb_flag = absorbFlag
		self.direction = direction

	#------------------------- Additional PostGIS methods -------------------------------#
	def createDBConnection(self):
		# create DB connection - fail if connection details invalid
		try:
			self.dbConn = psycopg2.connect( database = self.pgDetails['database'],
											user = self.pgDetails['user'],
											password = self.pgDetails['password'],
											host = self.pgDetails['host'],
											port = self.pgDetails['port'])
			self.dbConn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_READ_COMMITTED)
		except Exception as e:
			QgsMessageLog.logMessage("PostGreSQL database connection details invalid. {0}".format(e), level=QgsMessageLog.CRITICAL)
			raise Exception("PostGreSQL database connection details invalid. {0}".format(e))


	def createTable(self, layer, fieldList):
		# create new PostGIS table to write the results to
		dp = layer.dataProvider()
		tmp = self.pgDetails['table'].split('.')
		
		# add double quotes if schema/table in upper case
		if len(tmp) == 1:
			schema = 'public'
			lowerCase = (self.pgDetails['table'] == self.pgDetails['table'].lower())			
			if lowerCase:
				table = self.pgDetails['table']
			else:
				table = '"{0}"'.format(self.pgDetails['table'])
		else:		
			lowerCase = (tmp[0] == tmp[0].lower())			
			if lowerCase:
				schema = tmp[0]
			else:
				schema = '"{0}"'.format(tmp[0])
			
			lowerCase = (tmp[1] == tmp[1].lower())			
			if lowerCase:
				table = tmp[1]
			else:
				table = '"{0}"'.format(tmp[1])
		
		# build SQL command and field name string
		sqlCmd = []
		sqlCmd.append('id serial primary key')	
		sqlCmd.append('geom geometry(Polygon,{0})'.format(self.crs))
		fieldNames = []	
		fieldNames.append('geom')
		idList = ['id']
		idCount = 0           
		for field in fieldList:
			if dp.name() == 'postgres' and field.typeName() != '':
				if field.length() == -1:
					fType = field.typeName()
				else:
					fType = '{0}({1})'.format(field.typeName(), field.length())
			else:
				if field.type() == 1:
					fType = 'boolean'
				elif field.type() == 2:
					if field.typeName() == 'bit':
						fType = 'boolean'
					else: 
						fType = 'int'
				elif field.type() == 4:
					if field.length() > 10:
						fType = 'bigint'
					else:
						fType = 'int'
				elif field.type() == 6:
					fType = 'float8'
				elif field.type() == 10:
					if field.typeName() == 'bool' or field.typeName() == 'boolean':
						fType = 'boolean'
					else:
						fType = 'varchar'
						if field.length() > 0:
							fType += '({0})'.format(field.length())
				elif field.type() == 14:
					fType = 'date'
				elif field.type() == 16:
					fType = 'timestamp'
   			# ensure we do not get multiple id columns with same name    
			if field.name().lower().startswith('id'):
				if field.name().lower() in idList:
					newName = '{0}_{1}'.format(field.name().lower(), idCount)        
					while newName in idList:
						idCount += 1
						newName = '{0}_{1}'.format(field.name().lower(), idCount)
					idList.append(newName)
					fieldNames.append(newName)
					sqlCmd.append('{0} {1}'.format(newName, fType))
				else:
					idList.append(field.name().lower())
					fieldNames.append(field.name().lower())
					sqlCmd.append('{0} {1}'.format(field.name().lower(), fType))
			else:
				fieldNames.append(field.name().lower())
				sqlCmd.append('{0} {1}'.format(field.name().lower(), fType))
		createCmd = 'CREATE TABLE {0}.{1} ({2})'.format(schema,table,','.join(sqlCmd))
		self.fieldStr = ','.join(fieldNames)
		
		try:
			curs = self.dbConn.cursor()
			curs.execute(createCmd)
			curs.close() 
			self.dbConn.commit()			
		except Exception as e:
			self.dbConn.close()	
			QgsMessageLog.logMessage("Could not create PostGIS table. {0} Command text: {1}".format(e, createCmd), level=QgsMessageLog.CRITICAL)
			raise Exception("Could not create PostGIS table. {0}".format(e))

# ---------------------------------- MAIN FUNCTION ------------------------------------- #

	def work(self):
		"""
		* Actually do the processing
		"""
		
		# setup for progress bar and message
		self.toggle_show_cancel.emit(True)
		self.toggle_show_progress.emit(True)
		self.set_message.emit('Dividing Polygons')
		
		# get fields from the input shapefile
		layer = self.layer
		fieldList = layer.fields()

		# add new fields for this tool
		if fieldList.fieldNameIndex('POLY_ID') == -1:
			fieldList.append(QgsField('POLY_ID',QVariant.Int))
		if fieldList.fieldNameIndex('UNIQUE_ID') == -1:
			fieldList.append(QgsField('UNIQUE_ID', QVariant.String))
		if fieldList.fieldNameIndex('AREA') == -1:
			fieldList.append(QgsField('AREA', QVariant.Double))
		if fieldList.fieldNameIndex('POINTX') == -1:
			fieldList.append(QgsField('POINTX',QVariant.Int))
		if fieldList.fieldNameIndex('POINTY') == -1:
			fieldList.append(QgsField('POINTY',QVariant.Int))

		if self.outputType == 'PostGIS':
			self.createDBConnection()
			self.crs = layer.crs().postgisSrid()
			self.createTable(layer, fieldList)
			self.dbConn.close()
			writer = None
		else:
			# create a new shapefile to write the results to
			writer = QgsVectorFileWriter(self.outFilePath, "CP1250", fieldList, QGis.WKBPolygon, layer.crs(), "ESRI Shapefile")
		
		# how many features / sections will we have (for progress bar)
		iter = self.layer.getFeatures()
		featCount = 0
		totalArea = 0
		for feat in iter:
			featCount += 1
			totalArea += feat.geometry().area()
		del iter
		
		if self.chunk_size < 1:
			QgsMessageLog.logMessage("Chunk size invalid, defaulting to 50", level=QgsMessageLog.CRITICAL)
			self.chunk_size = 50
		
		# calculate no. of batches / progress bar divisions
		noBatches = featCount // self.chunk_size
		totalDivisions = totalArea // self.target_area
		
		# initialise progress counter
		noCompleted = 0

		# check if you've been killed		
		if self.killed:
			self.cleanup()
			raise UserAbortedNotification('USER Killed')

		# filter rows in layer
		if self.layer.dataProvider().name() == 'ogr':
			key_field = 'fid'
		else:
			uri = QgsDataSourceURI(self.layer.dataProvider().dataSourceUri())
			key_field = uri.keyColumn()
			if key_field == '':
				QgsMessageLog.logMessage("Layer must have a integer key column", level=QgsMessageLog.CRITICAL)
				raise Exception("Layer must have a integer key column. {0}".format(e))
			
			if self.layer.fields()[self.layer.fieldNameIndex(key_field)].type() != 2:
				QgsMessageLog.logMessage("Layer must have a integer key column", level=QgsMessageLog.CRITICAL)
				raise Exception("Layer must have an integer key column")
		
		for i in range(noBatches):
			batchMin = i * self.chunk_size
			batchMax = (i + 1) * self.chunk_size
			filtered = self.layer.setSubsetString('{0} > {1} AND {0} <= {2}'.format(key_field, batchMin, batchMax))
			if filtered:
			# run example worker
				self.example_worker = ExampleWorker(self, layer, fieldList, self.outputType, writer, self.pgDetails, self.target_area, self.absorb_flag, self.direction, totalDivisions, noCompleted)
				result = self.example_worker.work()
				if result != None:
					noCompleted = result
					QgsMessageLog.logMessage("No completed: " + str(noCompleted))
			else:
				QgsMessageLog.logMessage("Layer could not be filtered.", level=QgsMessageLog.CRITICAL)
				raise Exception("Layer could not be filtered.")
		
		if self.killed:
			self.cleanup()
			raise UserAbortedNotification('USER Killed')
		
		# finally, open the resulting file and return it
		if self.outputType == 'PostGIS':
			uri = QgsDataSourceURI()
			uri.setConnection(self.pgDetails['host'], self.pgDetails['port'], self.pgDetails['database'], self.pgDetails['user'], self.pgDetails['password'])
			tmp = self.pgDetails['table'].split('.')
			if len(tmp) == 1:		
				uri.setDataSource('public', self.pgDetails['table'],'geom','')
			else:
				uri.setDataSource(tmp[0], tmp[1],'geom','')	
				uri.setKeyColumn('id')
				uri.setWkbType(QgsWKBTypes.Polygon)		
				uri.setSrid(str(self.crs))
				layer = QgsVectorLayer(uri.uri(), 'Divided Polygon', 'postgres')
		else:
			layer = QgsVectorLayer(outFilePath, 'Divided Polygon', 'ogr')
		
		if layer.isValid():
			return layer
		else:
			return None


	def cleanup(self):
#		 print "cleanup here"
		try:
			self.dbConn.close() 
		except:
			pass
	
	def kill(self):
		self.killed = True
		try:
			self.example_worker.dbConn.close() 
		except:
			pass
		self.example_worker.killed = True
		self.set_message.emit('Aborting...')
		self.toggle_show_progress.emit(False)


class ExampleWorker():
	def __init__(self, parent, inLayer, fieldList, outputType, writer, pgDetails, targetArea, absorbFlag, direction, totalDivisions, noCompleted):

		# import args
		self.parent = parent
		self.layer = inLayer
		self.fieldList = fieldList
		self.outputType = outputType
		self.writer = writer
		self.pgDetails = pgDetails
		self.target_area = targetArea
		self.absorb_flag = absorbFlag
		self.direction = direction
		self.totalDivisions = totalDivisions
		self.noCompleted = noCompleted
		self.killed = False


	def brent(self, xa, xb, xtol, ftol, max_iter, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward):
		"""
		* Brent's Method, following Wikipedia's article algorithm.
		* 
		* - xa is the lower bracket of the interval of the solution we search.
		* - xb is the upper bracket of the interval of the solution we search.
		* - xtol is the minimum	 permitted width (in map units) of the interval before we give up.
		* - ftol is the required precision of the solution.
		* - max_iter is the maximum allowed number of iterations.
		* - geom is the geometry we are dividing.
		* - fixedCoord1 is the the first coordinate in the fixed dimension.
		* - fixedCoord2 is the the second coordinate in the fixed dimension.
		* - targetArea is the desired area of the section to cut off geom.
		* - horizontal or vertical cut - True / False respectively
		* - forward (left-right or bottom top) cut or backward (right-left or top-bottom) - True / False respectively
		"""

		## SET SOME VALUES

		# standard for iterative algorithms
		EPS = sys.float_info.epsilon

		## BASIC ERROR CHECKING (INTERVAL VALIDITY)

		# check that the bracket's interval is sufficiently big for this computer to work with.
		if abs(xb - xa) < EPS:
			raise BrentError("Initial bracket smaller than system epsilon.")

		# check lower bound
		fa = self.f(xa, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward)		 # first function call
		if abs(fa) < ftol:
			raise BrentError("Root is equal to the lower bracket")

		# check upper bound
		fb = self.f(xb, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward)		 # second function call
		if abs(fb) < ftol:
			raise BrentError("Root is equal to the upper bracket")
	
		# check if the root is bracketed.
		if fa * fb > 0.0:	# this is checking for different signs (to be sure we are either side of 0)
			raise BrentError("Root is not bracketed.")

		## START CALCULATION

		# if the area from a is smaller than b, switch the values
		if abs(fa) < abs(fb):
			xa, xb = xb, xa
			fa, fb = fb, fa

		# initialise c at a (therefore at the one with the biggest area to the right of it)
		xc, fc = xa, fa

		# init mflag
		mflag = True

		# do until max iterations is reached
		for i in range(max_iter):
 
			# try to calculate `xs` by using inverse quadratic interpolation...
			if fa != fc and fb != fc:
				xs = (xa * fb * fc / ((fa - fb) * (fa - fc)) + xb * fa * fc / ((fb - fa) * (fb - fc)) + xc * fa * fb / ((fc - fa) * (fc - fb)))
			else:
				# ...if you can't, use the secant rule.
				xs = xb - fb * (xb - xa) / (fb - fa)

			# check if the value of `xs` is acceptable, if it isn't use bisection.
			if ((xs < ((3 * xa + xb) / 4) or xs > xb) or
				(mflag == True	and (abs(xs - xb)) >= (abs(xb - xc) / 2)) or 
				(mflag == False and (abs(xs - xb)) >= (abs(xc - d) / 2)) or
				(mflag == True	and (abs(xb - xc)) < EPS) or 
				(mflag == False and (abs(xc - d)) < EPS)):

				# overwrite unacceptable xs value with result from bisection
				xs = (xa + xb) / 2
				mflag = True
			else:
				mflag = False
		
			## THE ABOVE BLOCK USED BRENT'S METHOD TO GET A SUGGESTED VALUE FOR S, THE BELOW BLOCK CHECKS IF IT IS GOOD, AND IF NOT SEEKS A NEW VALUE

			# get the value from f using the new xs value
			fs = self.f(xs, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward) # repeated function call

			# if the value (ideally 0) is less than the specified tolerance, return
			if abs(fs) < ftol:
				return xs

			# if the bracket has become smaller than the tolerance (but the value wasn't reached, something is wrong)
			# this can indicate the 'W' condition, where decreasing the interval increases the size of the resulting  area
			if abs(xb - xa) < xtol:
				raise BrentError("Bracket is smaller than tolerance.")

			# d is assigned for the first time here; it won't be used above on the first iteration because mflag is set
			d = xc	# it is just used in Brent's checks, not in the calculation per se
	
			# move c to b
			xc, fc = xb, fb
	
			# move one of the interval edges to the new point, such that zero remains within the interval
			# if the areas from a and s (current result) are same sign, move b to s, otherwise, move a to s
			if fa * fs < 0:			# different signs
				xb, fb = xs, fs
			else:					# same sign
				xa, fa = xs, fs

			# if the area from a is smaller than b, switch the values
			if abs(fa) < abs(fb):
				xa, xb = xb, xa
				fa, fb = fb, fa

		# this isn't as good (ran out of iterations), but seems generally fine
		return xs	# NB: increasing the number of iterations doesn't seem to get any closer


	def splitPoly(self,polygon, splitter, horizontal, forward):
		"""
		* Split a Polygon with a LineString
		* Returns exactly two polygons notionally referred to as being 'left' and 'right' of the cutline. 
		* The relationship between them is that the 'right' polygon (the chip) is notionally cut from the 'left' one (the potato).
		"""

		# need to take a deep copy	for the incoming polygon, as splitGeometry edits it directly...
		poly = QgsGeometry(polygon)

		# split poly (polygon) by splitter (line) http://gis.stackexchange.com/questions/114414/cannot-split-a-line-using-qgsgeometry-splitgeometry-in-qgis
		res, polys, topolist = poly.splitGeometry(splitter, False)
	
		# add poly (which might be a multipolygon) to the polys array
		if poly.isMultipart():
			## TODO: I think that this is where we might be getting that odd error on no absorb
			multiGeom = QgsGeometry()
			multiGeom = poly.asMultiPolygon()
			for i in multiGeom:
				polys.append(QgsGeometry().fromPolygon(i))
		else:
			# ...OR load the feature into a list of one (it may be extended in the course of splitting if we create noncontiguous offcuts) and loop through it
			polys.append(poly)
	
		## sort right, left, residual
	
		# verify that it worked and that more than one polygon was returned
		if res == 0 and len(polys) >1:
			if forward: ### from bottom left
				if horizontal:	## cut from the bottom

					# left is the top one
					maxy = float('-inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().yMaximum()
						if p > maxy:
							maxy = p
							maxyi = i
					left = polys.pop(maxyi)
	
					# right is the bottom one
					miny = float('inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().yMinimum()
						if p < miny:
							miny = p
							minyi = i
						elif p == miny:		# if there is a tie for which is the rightest, get the rightest in the other dimension
							if polys[i].boundingBox().xMinimum() < polys[minyi].boundingBox().xMinimum():	# left
								minyi = i
					right = polys.pop(minyi)
	
				else:	## cutting from the left
	
					# left is the rightest one
					maxx = float('-inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().xMaximum()
						if p > maxx:
							maxx = p
							maxxi = i
					left = polys.pop(maxxi)
	
					# right is the leftest one
					minx = float('inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().xMinimum()
						if p < minx:
							minx = p
							minxi = i
						elif p == minx:		# if there is a tie for which is the rightest, get the rightest in the other dimension
							if polys[i].boundingBox().yMinimum() < polys[minxi].boundingBox().yMinimum():	# bottom
								minxi = i
					right = polys.pop(minxi)
		
			else:	### cut from top / right (forward_flag is false)

				if horizontal:	## cut from the top

					# left is the bottom one
					miny = float('inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().yMinimum()
						if p < miny:
							miny = p
							minyi = i
					left = polys.pop(minyi)
	
					# right is the top one
					maxy = float('-inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().yMaximum()
						if p > maxy:
							maxy = p
							maxyi = i
						elif p == maxy:		# if there is a tie for which is the rightest, get the rightest in the other dimension
							if polys[i].boundingBox().xMaximum() > polys[maxyi].boundingBox().xMaximum():
								maxyi = i
					right = polys.pop(maxyi)
	
				else:	## cutting from the right
	
					# left is the leftest one
					minx = float('inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().xMinimum()
						if p < minx:
							minx = p
							minxi = i
					left = polys.pop(minxi)
		
					# right is the rightest one
					maxx = float('-inf')
					for i in range(len(polys)):
						p = polys[i].boundingBox().xMaximum()
						if p > maxx:
							maxx = p
							maxxi = i
						elif p == maxx:		# if there is a tie for which is the rightest, get the rightest in the other dimension
							if polys[i].boundingBox().yMaximum() > polys[maxxi].boundingBox().yMaximum():
								maxxi = i
					right = polys.pop(maxxi)

			# work out if any remaining polygons are contiguous with left or not
			contiguous = []
			noncontiguous = []
			if len(polys) > 0:
				for j in polys:
					if left.touches(j):
						contiguous.append(j)
					else:
						noncontiguous.append(j)

				# join all contiguous parts back to left
				if len(contiguous) > 0:
					contiguous += [left]
					left = QgsGeometry.unaryUnion(contiguous)
		
			# return the two sections (left is the potato, right is the chip...), plus any noncontiguous polygons
			return left, right, noncontiguous
		else:
			# log error
			QgsMessageLog.logMessage("FAIL: Polygon division failed.", level=QgsMessageLog.CRITICAL)
			

	def getSliceArea(self,sliceCoord, poly, fixedCoord1, fixedCoord2, horizontal, forward):
		"""
		* Splits a polygon to a defined distance from the minimum value for the selected dimension and
		*  returns the area of the resultng polygon
		"""

		# construct a list of points representing a line by which to split the polygon
		if horizontal:
			splitter = [QgsPoint(fixedCoord1, sliceCoord), QgsPoint(fixedCoord2, sliceCoord)]	# horizontal split
		else:
			splitter = [QgsPoint(sliceCoord, fixedCoord1), QgsPoint(sliceCoord, fixedCoord2)] # vertical split

		# split the polygon
		left, right, residual = self.splitPoly(poly, splitter, horizontal, forward)
	
		# return the area of the bit you cut off
		return right.area()


	def f(self,sliceCoord, poly, fixedCoord1, fixedCoord2, targetArea, horizontal, forward):
		"""
		* Split a Polygon with a LineString, returning the area of the polygon to the right of the split
		* returns the area of the polygon on the right of the splitter line
		"""
	
		# return the difference between the resulting polygon area (right of the line) and the desired area
		return self.getSliceArea(sliceCoord, poly, fixedCoord1, fixedCoord2, horizontal, forward) - targetArea

	#------------------------- Additional PostGIS methods -------------------------------#
	def createDBConnection(self):
		# create DB connection - fail if connection details invalid
		try:
			self.dbConn = psycopg2.connect( database = self.pgDetails['database'],
											user = self.pgDetails['user'],
											password = self.pgDetails['password'],
											host = self.pgDetails['host'],
											port = self.pgDetails['port'])
			self.dbConn.set_isolation_level(psycopg2.extensions.ISOLATION_LEVEL_READ_COMMITTED)
		except Exception as e:
			QgsMessageLog.logMessage("PostGreSQL database connection details invalid. {0}".format(e), level=QgsMessageLog.CRITICAL)
			raise Exception("PostGreSQL database connection details invalid. {0}".format(e))


	def writeFeature(self, feat):
		sqlValues = []
		sqlValues.append('ST_GeomFromText(\'{0}\', {1})'.format(feat.geometry().exportToWkt(),self.crs))
		fields = feat.fields()
		attributes = feat.attributes()
		for n in range(len(fields)):
			if attributes[n] == NULL:
				sqlValues.append('{0}'.format(attributes[n]))
			else:
				if fields[n].type() == 2:
					if fields[n].typeName() == 'bit':
						if attributes[n] == 1 or attributes[n] == -1:
							sqlValues.append('true')
						elif attributes[n] == 0:
							sqlValues.append('false')
						else:				
							sqlValues.append('NULL')
					else:
						sqlValues.append('{0}'.format(attributes[n]))
				elif fields[n].type() == 4 or fields[n].type() == 6:
					sqlValues.append('{0}'.format(attributes[n]))
				elif fields[n].type() == 10:
					if fields[n].typeName() == 'bool':
						if attributes[n] == 't':
							sqlValues.append('true')
						elif attributes[n] == 'f':
							sqlValues.append('false')
						else:
							sqlValues.append('NULL')
					else:
						# insert value handling apostrophes
						sqlValues.append('\'{0}\''.format(attributes[n].replace(u'\u2019','\'').replace("'","\'\'")))
				elif fields[n].type() == 14: # date
					sqlValues.append('\'{0}\''.format(attributes[n].toString('yyyy-MM-dd')))
				elif fields[n].type() == 16: # date/time
					sqlValues.append('\'{0}\''.format(attributes[n].toString('yyyy-MM-dd hh:mm:ss')))
			
		insertCmd = 'INSERT INTO {0} ({1}) VALUES({2})'.format(self.pgDetails['table'], self.parent.fieldStr, ','.join(sqlValues))
		
		try:
			curs = self.dbConn.cursor()
			curs.execute(insertCmd)
			curs.close()
		except Exception as e:
			try:
				curs.close()	
			except:
				pass
			QgsMessageLog.logMessage("Feature could not be written to the output table. {0} Command text: {1}".format(e, insertCmd), level=QgsMessageLog.CRITICAL)

#---------------------------------------------------------------------- CL #

# ---------------------------------- MAIN FUNCTION ------------------------------------- #

	def work(self):
		"""
		* Actually do the processing
		"""
		
		# TODO: reference these properly
		layer = self.layer
		fieldList = self.fieldList
		target_area = self.target_area
		absorb_flag = self.absorb_flag
		direction = self.direction
  
		# initial settings
		t = 0.1				# tolerance for function rooting - this is flexible now it has been divorced from the buffer
		buffer = 1e-6		# this is the buffer to ensure that an intersection occurs

		# set the direction (horizontal or vertical)
		if direction < 2:
			horizontal_flag = True 
		else:
			horizontal_flag = False

		# True = cut from bottom / left, False = cut from top / right
		if (direction == 0 or direction == 2): 
			forward_flag = True
		else:
			forward_flag = False	

		# this is used to make sure we don't hit an insurmountable error and just repeatedly change direction
		ERROR_FLAG_0 = False	# tracks if increasing number of subdivisions failed 
		ERROR_FLAG_1 = False	# tracks decreasing number of subdivisions to try and work around an error
		ERROR_FLAG_2 = False	# tracks change direction from forward to backward (or vice versa) by switching forward_flag
		ERROR_FLAG_3 = False	# tracks change cutline from horizontal to backward (or vice versa) by switching horizontal_flag

		if self.outputType == 'PostGIS':
			self.createDBConnection()
			self.crs = layer.crs().postgisSrid()              

		# define this to ensure that it's global
		subfeatures = []

		# init feature counter (for ID's)
		j = self.noCompleted
		# init feature counter to save batches  
		x = 0
		totalDivisions = self.totalDivisions
		# used to control progress bar (only send signal for an increase)
		currProgress = int((j*1.0) / totalDivisions * 100)
  		QgsMessageLog.logMessage("Initialised current progress to: " + str(currProgress))

		# check if you've been killed		
		if self.killed:
			self.cleanup()
			raise UserAbortedNotification('USER Killed')

		# loop through all of the features in the input data
		 
		iter = layer.getFeatures()
		for feat in iter:
			# check if you've been killed
			if self.killed:
				break
			
			# verify that it is a polygon
			if feat.geometry().wkbType() == QGis.WKBPolygon or feat.geometry().wkbType() == QGis.WKBMultiPolygon:	## if geom.type() == QGis.Polygon:

				# get the attributes to write out
				currAttributes = feat.attributes()

				# extract the geometry and sort out self intersections etc. with a buffer of 0m
				bufferedPolygon = feat.geometry().buffer(0, 15)

				# if the buffer came back as None, skip
				if bufferedPolygon is None:
					QgsMessageLog.logMessage("A polygon could not be buffered by QGIS, ignoring", level=QgsMessageLog.WARNING)
					continue

				# make multipolygon into list of polygons...
				subfeatures = []
				if bufferedPolygon.isMultipart():
					multiGeom = QgsGeometry()
					multiGeom = bufferedPolygon.asMultiPolygon()
					for i in multiGeom:
						subfeatures.append(QgsGeometry().fromPolygon(i))
				else:
					# ...OR load the feature into a list of one (it may be extended in the course of splitting if we create noncontiguous offcuts) and loop through it
					subfeatures.append(bufferedPolygon) 
	
				#loop through the geometries
				for poly in subfeatures:
					# check if you've been killed
					if self.killed:
						break
					
					# how many polygons are we going to have to chop off?
					nPolygons = int(poly.area() // target_area)

					# (needs to be at least 1...)
					if nPolygons == 0:
						nPolygons = 1

					# adjust the targetArea to reflect absorption if required
					if absorb_flag:
						targetArea = target_area + ((poly.area() % target_area) / nPolygons)
					else:
						targetArea = target_area

					# work out the size of a square with area = targetArea if required
					sq = sqrt(targetArea)

					# until there is no more dividing to do...						
					while poly.area() > targetArea + t:
						# check if you've been killed
						if self.killed:
							break
						
						# the bounds are used for the interval
						boundsR = poly.geometry().boundingBox()
						bounds = [boundsR.xMinimum(), boundsR.yMinimum(), boundsR.xMaximum(), boundsR.yMaximum()]

						# get interval and fixed coordinates (buffer otherwise there won't be an intersection between polygon and cutline!)
						if horizontal_flag:
							interval = bounds[1] + buffer, bounds[3] - buffer
							fixedCoords = bounds[0], bounds[2]
						else:
							interval = bounds[0] + buffer, bounds[2] - buffer
							fixedCoords = bounds[1], bounds[3]

						# is the interval larger than the required square? (We know the required area is > target+t because of the while statement)
						if (interval[1]-interval[0]) > sq: 

							# this is the resulting area of a slice the width of sq from the polygon
							if forward_flag:
								sqArea = self.getSliceArea(interval[0] + sq - buffer, poly, fixedCoords[0], fixedCoords[1], horizontal_flag, forward_flag)	# cutting from bottom/left
							else:
								sqArea = self.getSliceArea(interval[1] - sq + buffer, poly, fixedCoords[0], fixedCoords[1], horizontal_flag, forward_flag)	# cutting from top/right

							# what is the nearest number of subdivisions of targetArea that could be extracted from that slice?
							nSubdivisions = int(round(sqArea / targetArea))

							# if the answer is 0, make it 1...
							if nSubdivisions == 0:
								nSubdivisions = 1
		
							# make a backup copy to reset if we move from ERROR_FLAG_0 to ERROR_FLAG_1
							nSubdivisions2 = nSubdivisions
		
							'''now use brent's method to find the optimal coordinate in the variable dimension (e.g. the y coord for a horizontal cut)'''
		
							# if it fails, try increasing nSubdivisions (k) until it works or you get a different error
							while True:
								# check if you've been killed
								if self.killed:
									break
	
								# how big must the target area be to support this many subdivisions?
								initialTargetArea = nSubdivisions * targetArea
	
								# try to split using this new value
								try:
									# try to zero the equation
									result = self.brent(interval[0], interval[1], 1e-6, t, 500, poly, fixedCoords[0], fixedCoords[1], initialTargetArea, horizontal_flag, forward_flag)

									# if it worked (no exception raised) then exit this while loop and carry on
									break

								# if it didn't work...
								except BrentError as e:
				
									# is it a W condition error?
									if e.value == "Bracket is smaller than tolerance.":
										# ...increase number of subdivisions and go around again
										nSubdivisions += 1											
										continue
				
									# if not a W condition error, just move on
									else:
				
										# set flag and stop trying to adjust nSubdivisions
										ERROR_FLAG_0 = True
										break
					
							# if that didn't work, try decreasing instead of increasing								
							if ERROR_FLAG_0:
		
								# log message
								QgsMessageLog.logMessage("Increasing number of subdivisions didn't work, try decreasing... (Division)", level=QgsMessageLog.WARNING)
		
								nSubdivisions = nSubdivisions2	# reset
								limit = 1
								while nSubdivisions >= limit:
									# check if you've been killed
									if self.killed:
										break
		
									# set the flag if it's the last time around
									if nSubdivisions == limit:
										ERROR_FLAG_1 = True
	
									# how big must the target area be to support this many subdivisions?
									initialTargetArea = nSubdivisions * targetArea
	
									# try to split using this new value
									try:
										# try to zero the equation
										result = self.brent(interval[0], interval[1], 1e-6, t, 500, poly, fixedCoords[0], fixedCoords[1], initialTargetArea, horizontal_flag, forward_flag)

										# if it worked (no exception raised) then exit this while loop and carry on
										break

									# if it didn't work...
									except BrentError as e:
		
										# ...increase number of subdivisions and go around again
										nSubdivisions -= 1											
										continue
		
							# if increasing the subdivision size didn't help, then start trying shifting directions
							if ERROR_FLAG_1:
		
								# these need resetting here otherwise it won't try to cut again, just skip to the next error!
								ERROR_FLAG_0 = False
								ERROR_FLAG_1 = False
		
								# log message
								QgsMessageLog.logMessage("Decreasing number of subdivisions didn't work, try playing with direction... (Division)", level=QgsMessageLog.WARNING)
			
								# switch the movement direction 
								if ERROR_FLAG_2 == False:
				
									# log that this has been tried
									ERROR_FLAG_2 = True
									QgsMessageLog.logMessage("Reversing movement direction (Division)", level=QgsMessageLog.WARNING)

									# reverse the direction of movement and try again
									forward_flag = not forward_flag
									continue

								# if the above didn't work, switch the direction of the cutline
								elif ERROR_FLAG_3 == False:
			
									# un-log 2, meaning that it will run again and so try the 4th direction
									ERROR_FLAG_2 = False
				
									# log that this has been tried
									ERROR_FLAG_3 = True
									QgsMessageLog.logMessage("Reversing cutline direction (Division)", level=QgsMessageLog.WARNING)

									# reverse the cutline direction and try again
									horizontal_flag = not horizontal_flag
									continue
							
								# if none of the above worked, just skip it and move to a new feature
								else:
							
									## WRITE THE UNSPLITTABLE POLYGON TO THE SHAPEFILE ANYWAY

									# make a feature with the right schema
									fet = QgsFeature()
									fet.setFields(fieldList)
				
									# populate inherited attributes
									for a in range(len(currAttributes)):
										fet[a] = currAttributes[a]
		
									# calculate representative point
									pt = poly.pointOnSurface().asPoint()
				
									# populate new attributes
									fet.setAttribute('POLY_ID', j)
									fet.setAttribute('UNIQUE_ID', str(uuid4()))
									fet.setAttribute('AREA', poly.area())
									fet.setAttribute('POINTX', pt[0])
									fet.setAttribute('POINTY', pt[1])
				
									# add the geometry to the feature
									fet.setGeometry(poly)
				
									# write the feature to the Shapefile/PostGIS table
									if self.outputType == 'PostGIS':
										self.writeFeature(fet)
									else:
										# write the feature to the out file
										writer.addFeature(fet)
				
									# increment feature counters 
									j+=1
									x+=1
															
									# commit to PostgreSQL if required
									if self.outputType == 'PostGIS' and k == self.pgDetails['batch_size']:
										self.dbConn.commit()
										x = 0 
									 
									# update progress bar if required
									if j // totalDivisions * 100 > currProgress:
										self.parent.progress.emit(j // totalDivisions * 100)
							
									# log that there was a problem
									QgsMessageLog.logMessage("There was an un-dividable polygon in this dataset.", level=QgsMessageLog.WARNING)
									
									# on to the next one, hopefully with more luck!
									continue

							# if it worked, reset the flags
							ERROR_FLAG_0 = False
							ERROR_FLAG_1 = False
							ERROR_FLAG_2 = False
							ERROR_FLAG_3 = False

							# create the desired cutline as lists of QgsPoints
							if horizontal_flag:
								line = [QgsPoint(fixedCoords[0], result), QgsPoint(fixedCoords[1], result)] # horizontal split
							else:
								line = [QgsPoint(result, fixedCoords[0]), QgsPoint(result, fixedCoords[1])] # vertical split

							# calculate the resulting polygons - poly will be sliced again, initialSlice will be subdivided
							poly, initialSlice, residuals = self.splitPoly(poly, line, horizontal_flag, forward_flag)

							# put the residuals in the list to be processed
							subfeatures += residuals

						# bounds not bigger than sq, so no division necessary, just subdivide this last one directly (nothing will happen if it can't be subdivided)
						else:
	
							# set the remainder of the polygon as the final slice, and poly to an empty polygon
							initialSlice = poly
							poly = QgsGeometry.fromPolygon([[]])

							# what is the nearest number of subdivisions of targetArea that could be extracted from that slice? (must be at least 1)
							# TODO: verify this doesn't need rounding
							nSubdivisions = int(initialSlice.area() // targetArea) # shouldn't need rounding...
							if nSubdivisions == 0:
								nSubdivisions = 1							

						#...then divide that into sections of targetArea
						for k in range(nSubdivisions-1):	# nCuts = nPieces - 1
							# check if you've been killed
							if self.killed:
								break
							
							# the bounds are used for the interval
							sliceBoundsR = initialSlice.boundingBox()
							sliceBounds = [sliceBoundsR.xMinimum(), sliceBoundsR.yMinimum(), sliceBoundsR.xMaximum(), sliceBoundsR.yMaximum()]

							# get the slice direction (opposite to main direction)
							sliceHorizontal = not horizontal_flag

							if sliceHorizontal:
								# get interval and fixed coordinates
								sliceInterval = sliceBounds[1] + buffer, sliceBounds[3] - buffer # buffer otherwise there won't be an intersection between polygon and cutline!
								sliceFixedCoords = sliceBounds[0], sliceBounds[2]
							else:
								# get interval and fixed coordinates
								sliceInterval = sliceBounds[0] + buffer, sliceBounds[2] - buffer # buffer otherwise there won't be an intersection between polygon and cutline!
								sliceFixedCoords = sliceBounds[1], sliceBounds[3]

							# restore the tolerance (may be adjusted in the below loop)
							tol = t

							# infinite loop
							while True:
								# check if you've been killed
								if self.killed:
									break
								
								# brent's method to find the optimal coordinate in the variable dimension (e.g. the y coord for a horizontal cut)
								try:
								
									# search for result
									sliceResult = self.brent(sliceInterval[0], sliceInterval[1], 1e-6, tol, 500, initialSlice, sliceFixedCoords[0], sliceFixedCoords[1], targetArea, sliceHorizontal, forward_flag)
			
									# stop searching if result is found
									break
								except BrentError as e:
				
									# if it is a W condition error, double the tolerance
									if e.value == "Bracket is smaller than tolerance.":
										#QgsMessageLog.logMessage(e.value + ": increasing tolerance (Subdivision)", level=QgsMessageLog.WARNING)
				
										# double the tolerance and try again
										tol *= 2
										continue
					
									# otherwise, give up and try something else
									else:
				
										# set the flag that this has been tried and move on
										ERROR_FLAG_1 = True
										break
					
							## if the above didn't work then we need to try some more drastic measures

							# try reversing the movement direction
							if ERROR_FLAG_1 and not ERROR_FLAG_2: # (NB: Subdivision does not use Errorflag 0)

								# log that this has been tried
								ERROR_FLAG_2 = True
								QgsMessageLog.logMessage("Reversing movement direction (Subdivision)", level=QgsMessageLog.WARNING)

								# reverse the direction of movement and try again
								forward_flag = not forward_flag
								continue
			
							# if that didn't work, switch the direction of the cutline
							elif ERROR_FLAG_1 and not ERROR_FLAG_3:
			
								# re-reverse movement direction (as it didn't work)
		#												forward_flag = not forward_flag

								# log that this has been tried
								ERROR_FLAG_3 = True
								QgsMessageLog.logMessage("Reversing cutline direction (Subdivision):", level=QgsMessageLog.WARNING)

								# reverse the cutline direction and pass back to the outer division to try again in the opposite direction (otherwise we would get long thin strips, not squares)
								horizontal_flag = not horizontal_flag
								break	# this should mean that the 'else' for this statement will never be reached
			
						# re-reverse movement direction (as it didn't work)
		#										if ERROR_FLAG_2:
		#											forward_flag = not forward_flag

							# if it worked, reset the flags
							ERROR_FLAG_1 = False
							ERROR_FLAG_2 = False
							ERROR_FLAG_3 = False
					
							# create the desired cutline as lists of QgsPoints
							if horizontal_flag:
								sliceLine = [QgsPoint(sliceResult, sliceFixedCoords[0]), QgsPoint(sliceResult, sliceFixedCoords[1])]	# horizontal split
							else:
								sliceLine = [QgsPoint(sliceFixedCoords[0], sliceResult), QgsPoint(sliceFixedCoords[1], sliceResult)]	# vertical split

							# calculate the resulting polygons - initialSlice becomes left (to be chopped again)
							initialSlice, right, residuals = self.splitPoly(initialSlice, sliceLine, sliceHorizontal, forward_flag)

							# put the residuals in the list to be processed
							subfeatures += residuals
					
							## WRITE TO SHAPEFILE

							# make a feature with the right schema
							fet = QgsFeature()
							fet.setFields(fieldList)
					
							# populate inherited attributes
							for a in range(len(currAttributes)):
								fet[a] = currAttributes[a]
					
							# calculate representative point
							pt = right.pointOnSurface().asPoint()
						
							# populate new attributes
							fet.setAttribute('POLY_ID', j)
							fet.setAttribute('UNIQUE_ID', str(uuid4()))
							fet.setAttribute('AREA', right.area())
							fet.setAttribute('POINTX', pt[0])
							fet.setAttribute('POINTY', pt[1])
					
							# add the geometry to the feature
							fet.setGeometry(right)
					
							# write the feature to the Shapefile/PostGIS table
							if self.outputType == 'PostGIS':
								self.writeFeature(fet)
							else:
								# write the feature to the out file
								writer.addFeature(fet)
					
							# increment feature counters 
							j+=1
							x+=1
													
							# commit to PostgreSQL if required
							if self.outputType == 'PostGIS' and x == self.pgDetails['batch_size']:
           							self.dbConn.commit()
								x = 0

							# update progress bar if required
							if int((j*1.0) / totalDivisions * 100) > currProgress:
								currProgress = int((j*1.0) / totalDivisions * 100)
								self.parent.progress.emit(currProgress)

						## WRITE ANY OFFCUT FROM SUBDIVISION TO SHAPEFILE

						# make a feature with the right schema
						fet = QgsFeature()
						fet.setFields(fieldList)
				
						# populate inherited attributes
						for a in range(len(currAttributes)):
							fet[a] = currAttributes[a]
				
						# calculate representative point
						pt = initialSlice.pointOnSurface().asPoint()
						
						# populate new attributes
						fet.setAttribute('POLY_ID', j)
						fet.setAttribute('UNIQUE_ID', str(uuid4()))
						fet.setAttribute('AREA', initialSlice.area())
						fet.setAttribute('POINTX', pt[0])
						fet.setAttribute('POINTY', pt[1])
				
						# add the geometry to the feature
						fet.setGeometry(initialSlice)
				
						# write the feature to the Shapefile/PostGIS table
						if self.outputType == 'PostGIS':
							self.writeFeature(fet)
						else:
							# write the feature to the out file
							writer.addFeature(fet)
				
						# increment feature counters 
						j+=1
						x+=1
												
						# commit to PostgreSQL if required
						if self.outputType == 'PostGIS' and x == self.pgDetails['batch_size']:
           						self.dbConn.commit()
							x = 0

						# update progress bar if required
						if int((j*1.0) / totalDivisions * 100) > currProgress:
							currProgress = int((j*1.0) / totalDivisions * 100)
							self.parent.progress.emit(currProgress)

					try:
			
						## WRITEANY OFFCUT FROM DIVISION TO SHAPEFILE

						# make a feature with the right schema
						fet = QgsFeature()
						fet.setFields(fieldList)
				
						# populate inherited attributes
						for a in range(len(currAttributes)):
							fet[a] = currAttributes[a]
				
						# calculate representative point
						pt = poly.pointOnSurface().asPoint()
						
						# populate new attributes
						fet.setAttribute('POLY_ID', j)
						fet.setAttribute('UNIQUE_ID', str(uuid4()))
						fet.setAttribute('AREA', poly.area())
						fet.setAttribute('POINTX', pt[0])
						fet.setAttribute('POINTY', pt[1])
				
						# add the geometry to the feature
						fet.setGeometry(poly)
				
						# write the feature to the Shapefile/PostGIS table
						if self.outputType == 'PostGIS':
          						self.writeFeature(fet)
						else:
							# write the feature to the out file
							writer.addFeature(fet)
					
						# increment feature counters 
						j+=1
						x+=1
												
						# commit to PostgreSQL if required
						if self.outputType == 'PostGIS' and x == self.pgDetails['batch_size']:
          						self.dbConn.commit()
          						x = 0

						# update progress bar if required
						if int((j*1.0) / totalDivisions * 100) > currProgress:
							currProgress = int((j*1.0) / totalDivisions * 100)
							self.parent.progress.emit(currProgress)

					except:
						# this just means that there is no offcut, which is no problem!
						pass
					
			else:
				QgsMessageLog.logMessage("Whoops! That dataset isn't polygons!", level=QgsMessageLog.CRITICAL)
				raise Exception("Whoops! That dataset isn't polygons!")
		
		if self.killed:
			self.cleanup()
			raise UserAbortedNotification('USER Killed')
		
		if self.outputType == 'PostGIS':
			self.dbConn.commit()

		return j


	def cleanup(self):
#		 print "cleanup here"
		try:
			self.dbConn.close() 
		except:
			pass	


class UserAbortedNotification(Exception):
	pass

'''

***************************** CORE WORKER METHODS ******************************

'''

def start_worker(worker, iface, message, with_progress=True):
	"""
	* Launch the core worker thread
	"""

	# configure the QgsMessageBar
	message_bar = iface.messageBar().createMessage(message)
	progress_bar = QProgressBar()
	progress_bar.setAlignment(Qt.AlignLeft | Qt.AlignVCenter)
	if not with_progress:
		progress_bar.setMinimum(0)
		progress_bar.setMaximum(0)
	cancel_button = QPushButton()
	cancel_button.setText('Cancel')
	cancel_button.clicked.connect(worker.kill)
	message_bar.layout().addWidget(progress_bar)
	message_bar.layout().addWidget(cancel_button)
	iface.messageBar().pushWidget(message_bar, iface.messageBar().INFO)

	# start the worker in a new thread
	# let Qt take ownership of the QThread
	thread = QThread(iface.mainWindow())
	worker.moveToThread(thread)
	
	worker.set_message.connect(lambda message: set_worker_message(
		message, message_bar))

	worker.toggle_show_progress.connect(lambda show: toggle_worker_progress(
		show, progress_bar))
		
	worker.toggle_show_cancel.connect(lambda show: toggle_worker_cancel(
		show, cancel_button))
		
	worker.finished.connect(lambda result: worker_finished(
		result, thread, worker, iface, message_bar))
		
	worker.error.connect(lambda e, exception_str: worker_error(
		e, exception_str, iface))
		
	worker.progress.connect(progress_bar.setValue)
	
	worker.message_bar = message_bar
	worker.progress_bar = progress_bar
	thread.started.connect(worker.run)
	
	thread.start()
	return thread, message_bar


def worker_finished(result, thread, worker, iface, message_bar):
	"""
	* Open the resulting file and clean up after the worker thread
	"""

	# remove widget from message bar
	iface.messageBar().popWidget(message_bar)
	if result is not None:
		# report the result
		iface.messageBar().pushMessage('Success!')
		worker.successfully_finished.emit(result)
		
		# add the result to the workspace
		QgsMapLayerRegistry.instance().addMapLayer(result)	
	
	else:
		# notify the user that something went wrong
		iface.messageBar().pushMessage('Sorry, something went wrong! See the message log for more information.', level=QgsMessageBar.CRITICAL)
		
	# clean up the worker and thread
	worker.deleteLater()
	thread.quit()
	thread.wait()
	thread.deleteLater()


def worker_error(e, exception_string, iface):
	# notify the user that something went wrong
	iface.messageBar().pushMessage('Something went wrong! See the message log for more information.', level=QgsMessageBar.CRITICAL, duration=3)
	QgsMessageLog.logMessage('Worker thread raised an exception: %s' % exception_string, 'SVIR worker', level=QgsMessageLog.CRITICAL)


def set_worker_message(message, message_bar_item):
	message_bar_item.setText(message)


def toggle_worker_progress(show_progress, progress_bar):
	progress_bar.setMinimum(0)
	if show_progress:
		progress_bar.setMaximum(100)
	else:
		# show an undefined progress
		progress_bar.setMaximum(0)

		
def toggle_worker_cancel(show_cancel, cancel_button):
	cancel_button.setVisible(show_cancel)


'''

********************************** STUFF FOR THE GUI *************************************

'''



class PolygonDivider:
	"""QGIS Plugin Implementation."""

	def __init__(self, iface):
		"""Constructor.

		:param iface: An interface instance that will be passed to this class
			which provides the hook by which you can manipulate the QGIS
			application at run time.
		:type iface: QgsInterface
		"""
		# Save reference to the QGIS interface
		self.iface = iface
		# initialize plugin directory
		self.plugin_dir = os.path.dirname(__file__)
		# initialize locale
		locale = QSettings().value('locale/userLocale')[0:2]
		locale_path = os.path.join(
			self.plugin_dir,
			'i18n',
			'PolygonDivider_{}.qm'.format(locale))

		if os.path.exists(locale_path):
			self.translator = QTranslator()
			self.translator.load(locale_path)

			if qVersion() > '4.3.3':
				QCoreApplication.installTranslator(self.translator)


		# Declare instance attributes
		self.actions = []
		self.menu = self.tr(u'&Polygon Divider')
		self.toolbar = self.iface.addToolBar(u'PolygonDivider')
		self.toolbar.setObjectName(u'PolygonDivider')
		

	# noinspection PyMethodMayBeStatic
	def tr(self, message):
		"""Get the translation for a string using Qt translation API.

		We implement this ourselves since we do not inherit QObject.

		:param message: String for translation.
		:type message: str, QString

		:returns: Translated version of message.
		:rtype: QString
		"""
		# noinspection PyTypeChecker,PyArgumentList,PyCallByClass
		return QCoreApplication.translate('PolygonDivider', message)


	def add_action(
		self,
		icon_path,
		text,
		callback,
		enabled_flag=True,
		add_to_menu=True,
		add_to_toolbar=True,
		status_tip=None,
		whats_this=None,
		parent=None):
		"""Add a toolbar icon to the toolbar.

		:param icon_path: Path to the icon for this action. Can be a resource
			path (e.g. ':/plugins/foo/bar.png') or a normal file system path.
		:type icon_path: str

		:param text: Text that should be shown in menu items for this action.
		:type text: str

		:param callback: Function to be called when the action is triggered.
		:type callback: function

		:param enabled_flag: A flag indicating if the action should be enabled
			by default. Defaults to True.
		:type enabled_flag: bool

		:param add_to_menu: Flag indicating whether the action should also
			be added to the menu. Defaults to True.
		:type add_to_menu: bool

		:param add_to_toolbar: Flag indicating whether the action should also
			be added to the toolbar. Defaults to True.
		:type add_to_toolbar: bool

		:param status_tip: Optional text to show in a popup when mouse pointer
			hovers over the action.
		:type status_tip: str

		:param parent: Parent widget for the new action. Defaults None.
		:type parent: QWidget

		:param whats_this: Optional text to show in the status bar when the
			mouse pointer hovers over the action.

		:returns: The action that was created. Note that the action is also
			added to self.actions list.
		:rtype: QAction
		"""

		# Create the dialog (after translation) and keep reference
		self.dlg = PolygonDividerDialog()

		icon = QIcon(icon_path)
		action = QAction(icon, text, parent)
		action.triggered.connect(callback)
		action.setEnabled(enabled_flag)

		if status_tip is not None:
			action.setStatusTip(status_tip)

		if whats_this is not None:
			action.setWhatsThis(whats_this)

		if add_to_toolbar:
			self.toolbar.addAction(action)

		if add_to_menu:
			self.iface.addPluginToMenu(
				self.menu,
				action)

		self.actions.append(action)

		return action


	def initGui(self):
		"""Create the menu entries and toolbar icons inside the QGIS GUI."""

		icon_path = ':/plugins/PolygonDivider/icon.png'
		self.add_action(
			icon_path,
			text=self.tr(u'Divide Polygons'),
			callback=self.run,
			parent=self.iface.mainWindow())
 
		# launch file browser for output file button - link to function
		self.dlg.btnBrowse.clicked.connect(self.select_output_file)


	def unload(self):
		"""Removes the plugin menu item and icon from QGIS GUI."""
		for action in self.actions:
			self.iface.removePluginMenu(
				self.tr(u'&Polygon Divider'),
				action)
			self.iface.removeToolBarIcon(action)
		# remove the toolbar
		del self.toolbar

	def select_output_file(self):
		"""
		* JJH: Open file browser
		"""

		# get filename from dialog	
		filename = QFileDialog.getSaveFileName(self.dlg, "Select output file ","", '*.shp')
	
		# verify that a name was selected
		if filename != "":

			# clear previous value
			self.dlg.outputFile.clear()

			# make sure that an extension was included
			if filename[-4:] != '.shp':
				filename += '.shp'

			# put the result in the text box on the dialog
			self.dlg.outputFile.setText(filename)

#--- CL : Extract PostGIS settings from QGIS PostgreSQL Connection
	def extract_PostGIS_connection_details(self):
	 
		 pgConName = self.dlg.cboPGConnection.currentText()

		 pgDetails = dict()
		 
		 s = QtCore.QSettings()
		 pgDetails['database'] = s.value("PostgreSQL/connections/{0}/database".format(pgConName), '')
		 if len(pgDetails['database']) == 0:
			raise Exception('The selected PostGIS connection details could not be found, please check your settings')
		 pgDetails['host'] = s.value("PostgreSQL/connections/{0}/host".format(pgConName), '')
		 pgDetails['port'] = s.value("PostgreSQL/connections/{0}/port".format(pgConName), '')
		 pgDetails['user'] = s.value("PostgreSQL/connections/{0}/user".format(pgConName), 'postgres')
		 pgDetails['password'] = s.value("PostgreSQL/connections/{0}/password".format(pgConName), '')
		 
		 if len(pgDetails['user']) == 0 or len(pgDetails['password']) == 0:
			uri = QgsDataSourceURI()
			uri.setConnection(pgDetails['host'], pgDetails['port'], pgDetails['database'], pgDetails['user'], pgDetails['password'])
			(success, pgDetails['user'], pgDetails['password']) = QgsCredentials.instance().get(uri.uri(), pgDetails['user'], pgDetails['password'])
			if not success:
				raise Exception('Either user name or password for PostgreSQL were not entered, please try again')
		
		 return pgDetails
	

	def startWorker(self, inLayer, outputType, outFilePath, pgDetails, chunkSize, targetArea, absorbFlag, direction):
		"""
		* JJH: Run the polygon division in a thread, feed back to progress bar
		"""
		
		worker = CoreWorker(self.iface, inLayer, outputType, outFilePath, pgDetails, chunkSize, targetArea, absorbFlag, direction)
		start_worker(worker, self.iface, 'Running the worker')
		

	def run(self):
		"""
		* Run method that performs all the real work
		"""

		# JJH: set up the dialog here--------------------------------------------
	
		# populate cboLayer with the active layers
		self.dlg.cboLayer.clear()	# need to clear here or it will add them all again every time the dialog is opened
		layers = QgsMapLayerRegistry.instance().mapLayers().values() # List all open map layers - visible or not
		layer_list = []
		for layer in layers:
			 layer_list.append(layer.name())
		self.dlg.cboLayer.addItems(layer_list)

		# populate cboPGConnection with configured PostgreSQL connections
		self.dlg.cboPGConnection.clear()	# need to clear here or it will add them all again every time the dialog is opened
		s = QtCore.QSettings()
		s.beginGroup('PostgreSQL/connections')
		for connectionName in s.childGroups():
			self.dlg.cboPGConnection.addItem(connectionName)
		s.endGroup()

		# populate cboCutDir with the possible directions
		self.dlg.cboCutDir.clear() # need to clear here or it will add them all again every time the dialog is opened
		self.dlg.cboCutDir.addItems(['left to right', 'right to left', 'bottom to top', 'top to bottom'])


		#----------------------------------------------------------------------JJH

		# show the dialog
		self.dlg.show()
		
		# Run the dialog event loop
		result = self.dlg.exec_()
		
		# See if OK was pressed
		if result:

			# JH: RUN THE TOOL------------------------------------------------
			
			# get user settings
			inLayer = layers[self.dlg.cboLayer.currentIndex()]
			outFilePath = None
			pgDetails = None
			if self.dlg.rbShapefile.isChecked():
				outputType = 'Shapefile'
				outFilePath = self.dlg.outputFile.text()
				if outFilePath == '':
					QgsMessageLog.logMessage("Output shapefile not specified.", level=QgsMessageLog.CRITICAL)
					raise Exception("Output shapefile not specified.")
			#--- CL : Get PostGIS settings
			if self.dlg.rbPostgreSQL.isChecked():
				outputType = 'PostGIS'
				pgDetails = self.extract_PostGIS_connection_details()
				pgDetails['table'] = self.dlg.pgTable.text()
				if pgDetails['table'] == '':
					QgsMessageLog.logMessage("PostgreSQL connection or table not specified.", level=QgsMessageLog.CRITICAL)
					raise Exception("PostgreSQL connection or table not specified")
				try:
					pgDetails['batch_size'] = int(self.dlg.batchSize.text())
				except:
					pgDetails['batch_size'] = 10
			chunkSize = int(self.dlg.chunkSize.text())
			targetArea = float(self.dlg.targetArea.text())
			absorbFlag = self.dlg.chkOffcuts.isChecked()
			direction = self.dlg.cboCutDir.currentIndex()
		
			# run the tool
			self.startWorker(inLayer, outputType, outFilePath, pgDetails, chunkSize, targetArea, absorbFlag, direction)

			#--------------------------------------------------------------JJH
