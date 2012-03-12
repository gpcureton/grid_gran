#!/usr/bin/env python

import numpy as np
from matplotlib.colors import Colormap, normalize, LinearSegmentedColormap
from types import IntType, FloatType, ListType
from scipy import interpolate

def cmap_discretize(cmap, N):
	"""Return a discrete colormap from the continuous colormap cmap.

	cmap: colormap instance, eg. cm.jet. 
	N: Number of colors.

	Example
	x = resize(arange(100), (5,100))
	djet = cmap_discretize(cm.jet, 5)
	imshow(x, cmap=djet)
	"""

	cdict = cmap._segmentdata.copy()
	# N colors
	colors_i = np.linspace(0,1.,N)
	# N+1 indices
	indices = np.linspace(0,1.,N+1)
	for key in ('red','green','blue'):
		# Find the N colors
		D = np.array(cdict[key])
		I = interpolate.interp1d(D[:,0], D[:,1])
		#I = np.interp(indices,D[:,0], D[:,1])
		colors = I(colors_i)
		# Place these colors at the correct indices.
		A = np.zeros((N+1,3), float)
		A[:,0] = indices
		A[1:,1] = colors
		A[:-1,2] = colors
		# Create a tuple for the dictionary.
		L = []
		for l in A:
			L.append(tuple(l))
		cdict[key] = tuple(L)

	# Return colormap object.
	return LinearSegmentedColormap('colormap',cdict,1024)

def cmap_reverse(my_cmap):
	my_newdict = {}
	for items in my_cmap._segmentdata :
		item_tuple = my_cmap._segmentdata[items]
		new_item_tuple = []
		item_tuple_revIter = reversed(item_tuple)
		for segments in item_tuple :
			dummyTuple = item_tuple_revIter.next()
			new_item_tuple.append(tuple((1.-dummyTuple[0],dummyTuple[2],dummyTuple[1])))
		my_newdict[items] = tuple(new_item_tuple)

	return  LinearSegmentedColormap('my_newcolormap',my_newdict,256)

####################################################################
###   Inserts black segment into colourmap up to splitVal,       ###
###   compressing rest of colourmap towards high end. Colours    ###
###   in colourmap are preserved.                                ###
####################################################################

def cmap_preTruncate(my_cmap,splitVal):
	my_newdict = {}
	for items in my_cmap._segmentdata :
		item_tuple = my_cmap._segmentdata[items]
		new_item_tuple = []
		new_item_tuple.append(tuple((0.0,0.0,0.0)))
		new_item_tuple.append(tuple((splitVal,0.0,item_tuple[0][2])))

		item_tuple_Iter = iter(item_tuple)
		item_tuple_Iter.next()

		while 1 :
			try :
				dummyTuple = item_tuple_Iter.next()
				newVal = splitVal + dummyTuple[0]*(1.0 - splitVal)
				new_item_tuple.append(tuple((newVal,dummyTuple[2],dummyTuple[1])))
			except StopIteration :
				break

		my_newdict[items] = tuple(new_item_tuple)

	return  LinearSegmentedColormap('my_newcolormap',my_newdict,256)

####################################################################
###   Replaces colourmap upto splitVal with black, and rescales  ###
###   colourmap so that the transition occurs at startVal.       ###
###   Colours in original colourmap up to splitval are lost.     ###
####################################################################

def cmap_preTruncate2(my_cmap,startVal,splitVal):
	my_newdict = {}
	for items in my_cmap._segmentdata :
		item_tuple = my_cmap._segmentdata[items]
		colourTuple = []
		colourTupleNew = []

		# Create tuple containing tuples with colorbar coord > splitVal
		item_tuple_Iter = iter(item_tuple)
		while 1 :
			try :
				dummyTuple = item_tuple_Iter.next()
				if (dummyTuple[0] > splitVal) :
					colourTuple.append([dummyTuple[0],dummyTuple[1],dummyTuple[2]])
			except StopIteration :
				break

		# Rescale colorbar coord to range [ startVal-->1.0 ]
		oldRange = colourTuple[-1][0] - colourTuple[0][0]
		newRange = 1.0 - startVal

		for tup in np.arange(np.shape(colourTuple)[0]) :
			newCbarCoord = newRange*(colourTuple[tup][0]-colourTuple[0][0])/oldRange + startVal
			colourTupleNew.append([newCbarCoord,colourTuple[tup][1],colourTuple[tup][2]])

		colourTupleNew[0][1] = 0.0

		new_item_tuple = []
		new_item_tuple.append((0.0,0.0,0.0))
		colourTupleNew_iter = iter(colourTupleNew)

		while 1 :
			try :
				new_item_tuple.append(colourTupleNew_iter.next())
			except StopIteration :
				break

		my_newdict[items] = tuple(new_item_tuple)

	return  LinearSegmentedColormap('my_newcolormap',my_newdict,256)

