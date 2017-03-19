# ------------------------------ LOAD QGIS STANDALONE ---------------------------------- #

from PyQt4.QtCore import QVariant
from qgis.gui import QgsMessageBar
from qgis.core import *
import qgis.utils, sys
from uuid import uuid4
from math import sqrt

"""
* This script divides a polygon into squareish sections of a specified size
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

class BrentError(Exception):
	"""
	* Simple class for exceptions from Brent's Method.
	"""
	def __init__(self, value):
		self.value = value
	def __str__(self):
		return repr(self.value)


def brent(xa, xb, xtol, ftol, max_iter, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward):
	"""
	* Brent's Method, following Wikipedia's article algorithm.
	* 
	* - xa is the lower bracket of the interval of the solution we search.
	* - xb is the upper bracket of the interval of the solution we search.
	* - xtol is the minimum  permitted width (in map units) of the interval before we give up.
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
	fa = f(xa, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward)        # first function call
	if abs(fa) < ftol:
		raise BrentError("Root is equal to the lower bracket")

	# check upper bound
	fb = f(xb, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward)        # second function call
	if abs(fb) < ftol:
		raise BrentError("Root is equal to the upper bracket")
	
	# check if the root is bracketed.
	if fa * fb > 0.0: 	# this is checking for different signs (to be sure we are either side of 0)
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
			(mflag == True  and (abs(xs - xb)) >= (abs(xb - xc) / 2)) or 
			(mflag == False and (abs(xs - xb)) >= (abs(xc - d) / 2)) or
			(mflag == True  and (abs(xb - xc)) < EPS) or 
			(mflag == False and (abs(xc - d)) < EPS)):

			# overwrite unacceptable xs value with result from bisection
			xs = (xa + xb) / 2
			mflag = True
		else:
			mflag = False
		
		## THE ABOVE BLOCK USED BRENT'S METHOD TO GET A SUGGESTED VALUE FOR S, THE BELOW BLOCK CHECKS IF IT IS GOOD, AND IF NOT SEEKS A NEW VALUE

		# get the value from f using the new xs value
		fs = f(xs, geom, fixedCoord1, fixedCoord2, targetArea, horizontal, forward)	# repeated function call

		# if the value (ideally 0) is less than the specified tolerance, return
		if abs(fs) < ftol:
			return xs

		# if the bracket has become smaller than the tolerance (but the value wasn't reached, something is wrong)
		# this can indicate the 'W' condition, where decreasing the interval increases the size of the resulting  area
		if abs(xb - xa) < xtol:
			raise BrentError("Bracket is smaller than tolerance.")

		# d is assigned for the first time here; it won't be used above on the first iteration because mflag is set
		d = xc 	# it is just used in Brent's checks, not in the calculation per se
	
		# move c to b
		xc, fc = xb, fb
	
		# move one of the interval edges to the new point, such that zero remains within the interval
		# if the areas from a and s (current result) are same sign, move b to s, otherwise, move a to s
		if fa * fs < 0: 		# different signs
			xb, fb = xs, fs
		else:					# same sign
			xa, fa = xs, fs

		# if the area from a is smaller than b, switch the values
		if abs(fa) < abs(fb):
			xa, xb = xb, xa
			fa, fb = fb, fa

	# this isn't as good (ran out of iterations), but seems generally fine
	return xs	# NB: increasing the number of iterations doesn't seem to get any closer


def splitPoly(polygon, splitter, horizontal, forward):
	"""
	* Split a Polygon with a LineString
	* Returns exactly two polygons notionally referred to as being 'left' and 'right' of the cutline. 
	* The relationship between them is that the 'right' polygon (the chip) is notionally cut from the 'left' one (the potato).
	"""

	# need to take a deep copy  for the incoming polygon, as splitGeometry edits it directly...
	poly = QgsGeometry(polygon)

	# split poly (polygon) by splitter (line) http://gis.stackexchange.com/questions/114414/cannot-split-a-line-using-qgsgeometry-splitgeometry-in-qgis
	res, polys, topolist = poly.splitGeometry(splitter, False)
	
	# add poly (which might be a multipolygon) to the polys array
	if poly.isMultipart():
		print "multi!"
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
				maxy = 0
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
				maxx = 0
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
				maxy = 0
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
				maxx = 0
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
		

def getSliceArea(sliceCoord, poly, fixedCoord1, fixedCoord2, horizontal, forward):
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
	left, right, residual = splitPoly(poly, splitter, horizontal, forward)
	
	# return the area of the bit you cut off
	return right.area()


def f(sliceCoord, poly, fixedCoord1, fixedCoord2, targetArea, horizontal, forward):
	"""
	* Split a Polygon with a LineString, returning the area of the polygon to the right of the split
	* returns the area of the polygon on the right of the splitter line
	"""
	
	# return the difference between the resulting polygon area (right of the line) and the desired area
	return getSliceArea(sliceCoord, poly, fixedCoord1, fixedCoord2, horizontal, forward) - targetArea

# --------------------------------- USER SETTINGS -------------------------------------- #

def runSplit(self, layer, outFilePath, target_area, absorb_flag, direction, progress):

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
	ERROR_FLAG_0 = False 	# tracks if increasing number of subdivisions failed 
	ERROR_FLAG_1 = False 	# tracks decreasing number of subdivisions to try and work around an error
	ERROR_FLAG_2 = False 	# tracks change direction from forward to backward (or vice versa) by switching forward_flag
	ERROR_FLAG_3 = False	# tracks change cutline from horizontal to backward (or vice versa) by switching horizontal_flag

	# get fields from the input shapefile
	fieldList = layer.fields()
	
	# TODO: NEED TO CHECK IF THEY ALREADY EXIST
	# add new fields for this tool (if they don't already exist)
# 	if fieldList.lookupField('ps_id') == -1:
	fieldList.append(QgsField('ps_id',QVariant.Int))
# 	if fieldList.lookupField('ps_uuid') == -1:
	fieldList.append(QgsField('ps_uuid', QVariant.String))
# 	if fieldList.lookupField('ps_area') == -1:
	fieldList.append(QgsField('ps_area', QVariant.Double))
	
	# create a new shapefile to write the results to
	writer = QgsVectorFileWriter(outFilePath, "CP1250", fieldList, QGis.WKBPolygon, layer.crs(), "ESRI Shapefile")
	
	# define this to ensure that it's global
	subfeatures = []

	# init feature counter (for ID's)
	j = 0
	
	# how many sections will we have (for progress bar)
	iter = layer.getFeatures()
	totalArea = 0
	for feat in iter:
		totalArea += feat.geometry().area()
	totalDivisions = totalArea // target_area

	# loop through all of the features in the input data
	iter = layer.getFeatures()
	for feat in iter:

		# verify that it is a polygon
		if feat.geometry().wkbType() == QGis.WKBPolygon or feat.geometry().wkbType() == QGis.WKBMultiPolygon:	## if geom.type() == QGis.Polygon:
	
			# get the attributes to write out
			currAttributes = feat.attributes()
	
			# extract the geometry and sort out self intersections etc. with a buffer of 0m
			bufferedPolygon = feat.geometry().buffer(0, 15)
	
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
				
				print poly.area(), targetArea + t, (poly.area() > targetArea + t)

				# until there is no more dividing to do...						
				while poly.area() > targetArea + t:

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
							sqArea = getSliceArea(interval[0] + sq - buffer, poly, fixedCoords[0], fixedCoords[1], horizontal_flag, forward_flag)	# cutting from bottom/left
						else:
							sqArea = getSliceArea(interval[1] - sq + buffer, poly, fixedCoords[0], fixedCoords[1], horizontal_flag, forward_flag)	# cutting from top/right

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
		
							# how big must the target area be to support this many subdivisions?
							initialTargetArea = nSubdivisions * targetArea
		
							# try to split using this new value
							try:
								# try to zero the equation
								result = brent(interval[0], interval[1], 1e-6, t, 500, poly, fixedCoords[0], fixedCoords[1], initialTargetArea, horizontal_flag, forward_flag)

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
			
								# set the flag if it's the last time around
								if nSubdivisions == limit:
									ERROR_FLAG_1 = True
		
								# how big must the target area be to support this many subdivisions?
								initialTargetArea = nSubdivisions * targetArea
		
								# try to split using this new value
								try:
									# try to zero the equation
									result = brent(interval[0], interval[1], 1e-6, t, 500, poly, fixedCoords[0], fixedCoords[1], initialTargetArea, horizontal_flag, forward_flag)

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
								
							# if none of the above worked, dump to shapefile
							else:
								QgsMessageLog.logMessage("Tried everything but can't divide the polygons in this layer.", level=QgsMessageLog.CRITICAL)
								self.iface.messageBar().pushMessage("Sorry, I tried everything but I can't divide the polygons in this layer. Can you simplify it and try again?", level=QgsMessageBar.CRITICAL)
								raise BrentError("Failed.") # go and print the map out

						# if it worked, reset the flags
						ERROR_FLAG_0 = False
						ERROR_FLAG_1 = False
						ERROR_FLAG_2 = False
						ERROR_FLAG_3 = False

						# create the desired cutline as lists of QgsPoints
						if horizontal_flag:
							line = [QgsPoint(fixedCoords[0], result), QgsPoint(fixedCoords[1], result)]	# horizontal split
						else:
							line = [QgsPoint(result, fixedCoords[0]), QgsPoint(result, fixedCoords[1])] # vertical split

						# calculate the resulting polygons - poly will be sliced again, initialSlice will be subdivided
						poly, initialSlice, residuals = splitPoly(poly, line, horizontal_flag, forward_flag)
	
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
					for k in range(nSubdivisions-1): 	# nCuts = nPieces - 1

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
		
							# brent's method to find the optimal coordinate in the variable dimension (e.g. the y coord for a horizontal cut)
							try:									
			
								# search for result
								sliceResult = brent(sliceInterval[0], sliceInterval[1], 1e-6, tol, 500, initialSlice, sliceFixedCoords[0], sliceFixedCoords[1], targetArea, sliceHorizontal, forward_flag)
				
								# stop searching if result is found
								break
							except BrentError as e:
					
								# if it is a W condition error, double the tolerance
								if e.value == "Bracket is smaller than tolerance.":
									QgsMessageLog.logMessage(e.value + ": increasing tolerance (Subdivision)", level=QgsMessageLog.WARNING)
					
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
							QgsMessageLog.logMessage(e.value + ": Reversing movement direction (Subdivision)", level=QgsMessageLog.WARNING)

							# reverse the direction of movement and try again
							forward_flag = not forward_flag
							continue
				
						# if that didn't work, switch the direction of the cutline
						elif ERROR_FLAG_1 and not ERROR_FLAG_3:
				
							# re-reverse movement direction (as it didn't work)
	# 												forward_flag = not forward_flag

							# log that this has been tried
							ERROR_FLAG_3 = True
							QgsMessageLog.logMessage(e.value + ": Reversing cutline direction (Subdivision):", level=QgsMessageLog.WARNING)

							# reverse the cutline direction and pass back to the outer division to try again in the opposite direction (otherwise we would get long thin strips, not squares)
							horizontal_flag = not horizontal_flag
							break	# this should mean that the 'else' for this statement will never be reached
				
					# re-reverse movement direction (as it didn't work)
	# 										if ERROR_FLAG_2:
	# 											forward_flag = not forward_flag

						# if it worked, reset the flags
						ERROR_FLAG_1 = False
						ERROR_FLAG_2 = False
						ERROR_FLAG_3 = False
						
						# create the desired cutline as lists of QgsPoints
						if horizontal_flag:
							sliceLine = [QgsPoint(sliceResult, sliceFixedCoords[0]), QgsPoint(sliceResult, sliceFixedCoords[1])]	# horizontal split
						else:
							sliceLine = [QgsPoint(sliceFixedCoords[0], sliceResult), QgsPoint(sliceFixedCoords[1], sliceResult)] 	# vertical split

						# calculate the resulting polygons - initialSlice becomes left (to be chopped again)
						initialSlice, right, residuals = splitPoly(initialSlice, sliceLine, sliceHorizontal, forward_flag)
	
						# put the residuals in the list to be processed
						subfeatures += residuals
						
						## WRITE TO SHAPEFILE
	
						# make a feature with the right schema
						fet = QgsFeature()
						fet.setFields(fieldList)
						
						# populate inherited attributes
						for a in range(len(currAttributes)):
							fet[a] = currAttributes[a]
						
						# populate new attributes
						fet.setAttribute('ps_id', j)
						fet.setAttribute('ps_uuid', str(uuid4()))
						fet.setAttribute('ps_area', right.area())
						
						# add the geometry to the feature
						fet.setGeometry(right)
						
						# write the feature to the out file
						writer.addFeature(fet)
						
						# increment feature counter and 
						j+=1
						self.dlg.progressBar.setValue(j / totalDivisions * 100)
						progress.setValue(j / totalDivisions * 100)

	 				## WRITE ANY OFFCUT FROM SUBDIVISION TO SHAPEFILE
	
					# make a feature with the right schema
					fet = QgsFeature()
					fet.setFields(fieldList)
					
					# populate inherited attributes
					for a in range(len(currAttributes)):
						fet[a] = currAttributes[a]
					
					# populate new attributes
					fet.setAttribute('ps_id', j)
					fet.setAttribute('ps_uuid', str(uuid4()))
					fet.setAttribute('ps_area', initialSlice.area())
					
					# add the geometry to the feature
					fet.setGeometry(initialSlice)
					
					# write the feature to the out file
					writer.addFeature(fet)
					
					# increment feature counter and 
					j+=1
					self.dlg.progressBar.setValue(j / totalDivisions * 100)
					progress.setValue(j / totalDivisions * 100)

				try:
				
					## WRITE  ANY OFFCUT FROM DIVISION TO SHAPEFILE
	
					# make a feature with the right schema
					fet = QgsFeature()
					fet.setFields(fieldList)
					
					# populate inherited attributes
					for a in range(len(currAttributes)):
						fet[a] = currAttributes[a]
					
					# populate new attributes
					fet.setAttribute('ps_id', j)
					fet.setAttribute('ps_uuid', str(uuid4()))
					fet.setAttribute('ps_area', poly.area())
					
					# add the geometry to the feature
					fet.setGeometry(poly)
					
					# write the feature to the out file
					writer.addFeature(fet)
					
					# increment feature counter and 
					j+=1
					self.dlg.progressBar.setValue(j / totalDivisions * 100)
					progress.setValue(j / totalDivisions * 100)
					
				except:
					pass
		else:
			QgsMessageLog.logMessage("Whoops! That dataset isn't polygons!", level=QgsMessageLog.CRITICAL)
			self.iface.messageBar().pushMessage("Whoops! That dataset isn't polygons!", level=QgsMessageBar.CRITICAL)