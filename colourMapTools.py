import numpy as np
from matplotlib.colors import Colormap, normalize, LinearSegmentedColormap
import matplotlib.numerix as nx
from types import IntType, FloatType, ListType

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
    colors_i = numpy.linspace(0,1.,N)
    # N+1 indices
    indices = numpy.linspace(0,1.,N+1)
    for key in ('red','green','blue'):
        # Find the N colors
        D = numpy.array(cdict[key])
        I = numpy.interp(D[:,0], D[:,1])
        colors = I(colors_i)
        # Place these colors at the correct indices.
        A = zeros((N+1,3), float)
        A[:,0] = indices
        A[1:,1] = colors
        A[:-1,2] = colors
        # Create a tuple for the dictionary.
        L = []
        for l in A:
            L.append(tuple(l))
        cdict[key] = tuple(L)
    # Return colormap object.
    return matplotlib.colors.LinearSegmentedColormap('colormap',cdict,1024)

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


####################################################################
###              Example mainline for SentinelMap                ###
####################################################################

#import pylab
#import matplotlib.colors
#import colourMapTools as cmTools

#n=100

# create a random array
#X = nx.mlab.rand(n,n)
#cmBase = pylab.cm.jet

# plot it array as an image
#pylab.figure(1)
#pylab.imshow(X, cmap=cmBase, interpolation='nearest')

# define the sentinels
#sentinel1 = -10
#sentinel2 = 10

# replace some data with sentinels
#X[int(.1*n):int(.2*n), int(.5*n):int(.7*n)]  = sentinel1
#X[int(.6*n):int(.8*n), int(.2*n):int(.3*n)]  = sentinel2

# define the colormap and norm
#rgb1 = (0.,0.,0.)
#rgb2 = (1.,0.,0.)
#cmap = cmTools.SentinelMap(cmBase, sentinels={sentinel1:rgb1,sentinel2:rgb2,})
#norm = cmTools.SentinelNorm(ignore=[sentinel1,sentinel2])

# plot with the modified colormap and norm
#pylab.figure(2)
#pylab.imshow(X, cmap = cmap, norm=norm, interpolation='nearest')

#pylab.show()

#class SentinelMap(Colormap):
	#def __init__(self, cmap, sentinels={}):
		# boilerplate stuff
		#self.N = cmap.N
		#self.name = 'SentinelMap'
		#self.cmap = cmap
		#self.sentinels = sentinels
		#for rgb in sentinels.values():
			#if len(rgb)!=3:
				#raise ValueError('sentinel color must be RGB')


	#def __call__(self, scaledImageData, alpha=1):
		# assumes the data is already normalized (ignoring sentinels)
		# clip to be on the safe side
		#rgbaValues = self.cmap(nx.clip(scaledImageData, 0.,1.))

		#replace sentinel data with sentinel colors
		#for sentinel,rgb in self.sentinels.items():
			#r,g,b = rgb
			#rgbaValues[:,:,0] =  nx.where(scaledImageData==sentinel, r, rgbaValues[:,:,0])
			#rgbaValues[:,:,1] =  nx.where(scaledImageData==sentinel, g, rgbaValues[:,:,1])
			#rgbaValues[:,:,2] =  nx.where(scaledImageData==sentinel, b, rgbaValues[:,:,2])
			#rgbaValues[:,:,3] =  nx.where(scaledImageData==sentinel, alpha, rgbaValues[:,:,3])

		#return rgbaValues

	#_lut = rgbaValues

#class SentinelNorm(normalize):
        #"""
        #Leave the sentinel unchanged
        #"""
        #def __init__(self, ignore=[], vmin=None, vmax=None):
                #self.vmin=vmin
                #self.vmax=vmax

                #if type(ignore) in [IntType, FloatType]:
                        #self.ignore = [ignore]
                #else:
                        #self.ignore = list(ignore)


        #def __call__(self, value):

                #vmin = self.vmin
                #vmax = self.vmax

                #if type(value) in [IntType, FloatType]:
                        #vtype = 'scalar'
                        #val = array([value])
                #else:
                        #vtype = 'array'
                        #val = nx.asarray(value)

                # if both vmin is None and vmax is None, we'll automatically
                # norm the data to vmin/vmax of the actual data, so the
                # clipping step won't be needed.
                #if vmin is None and vmax is None:
                        #needs_clipping = False
                #else:
                        #needs_clipping = True

                #if vmin is None or vmax is None:
                        #rval = nx.ravel(val)
                        #do this if sentinels (values to ignore in data)
                        #if self.ignore:
                                #sortValues=nx.sort(rval)
                                #if vmin is None:
                                        # find the lowest non-sentinel value
                                        #for thisVal in sortValues:
                                                #if thisVal not in self.ignore:
                                                        #vmin=thisVal #vmin is the lowest non-sentinel value
                                                        #break
                                        #else:
                                                #vmin=0.
                                #if vmax is None:
                                        #for thisVal in sortValues[::-1]:
                                                #if thisVal not in self.ignore:
                                                        #vmax=thisVal #vmax is the greatest non-sentinel value
                                                        #break
                                        #else:
                                                #vmax=0.
                        #else:
                                #if vmin is None: vmin = min(rval)
                                #if vmax is None: vmax = max(rval)
                #if vmin > vmax:
                        #raise ValueError("minvalue must be less than or equal to maxvalue")
                #elif vmin==vmax:
                        #return 0.*value
                #else:
                        #if needs_clipping:
                                #val = nx.clip(val,vmin, vmax)
                        #result = (1.0/(vmax-vmin))*(val-vmin)

                # replace sentinels with original (non-normalized) values
                #for thisIgnore in self.ignore:
                        #result = nx.where(val==thisIgnore,thisIgnore,result)

                #if vtype == 'scalar':
                        #result = result[0]
                #return result

