#!/usr/bin/env python

### System libraries
import time, string, sys
from os import path

import numpy as np
from bisect import bisect,bisect_right,bisect_left

def binary_search_SP(cpsval,size,cpsArr):
	status = 0
	if ((cpsval > cpsArr[0]) and (cpsval < cpsArr[-1])):
		#-- This is for unit offset
		#l_index = 0
		#u_index = size+1
		#-- This is for zero offset
		l_index = -1
		u_index = size
		while ((u_index - l_index) > 1):
			search_index = (u_index + l_index)/2
			if (cpsval <= cpsArr[search_index]):
				u_index = search_index
			else :
				l_index = search_index

	elif (cpsval <= cpsArr[0]):
		print "cpsval <= minimum in LUT"
		l_index = 0
		u_index = 1
		status = 1
	elif (cpsval >= cpsArr[-1]):
		print "cpsval >= maximum in LUT"
		l_index = size-2
		u_index = size-1
		status = 1
	else :
		print "Some situation not captured, passing!!!"
		pass

	return	[u_index, l_index, status]

def COTconvert(nCOT,COTlut,CPS,COTvis):

	cpsArr = COTlut[:,0]
	visIrEfficRatio = COTlut[:,1]

	# If CPS exceeds largest value in LUT, assign to COTir
	# the last COT value in LUT

	if ( CPS > cpsArr[-1]) :
		s = visIrEfficRatio[-1]
	else :
		u_index,l_index,status = binary_search_SP(CPS, nCOT, cpsArr)

		print "Finished binary search..."
		print "l_index = ",l_index
		print "u_index = ",u_index

		m = (visIrEfficRatio[u_index]-visIrEfficRatio[l_index])/(cpsArr[u_index]-cpsArr[l_index])
		b = visIrEfficRatio[u_index]-m*cpsArr[u_index]
		s = m*CPS+b

		print "cpsArr[l_index] = ",cpsArr[l_index]
		print "cpsArr[u_index] = ",cpsArr[u_index]
		print "visIrEfficRatio[l_index] = ",visIrEfficRatio[l_index]
		print "visIrEfficRatio[u_index] = ",visIrEfficRatio[u_index]
		print "s = ",s

	# Compute the IR COT
	COTir = COTvis/s
	return COTir

def readCOTfile(file):
	results = []
	file=str(file)
	f = open(file,'r')

	lines = f.readlines()[1:]
	f.close()
	for l in lines:
		fields = l.split()
		CPS = float(fields[0])
		COT = float(fields[1])
		all = [CPS,COT]
		results.append(all)
	
	results = np.array(results)
	nCOT = np.shape(results)[0]
	return results,nCOT

def main(cot,cps):
	file = '/data/geoffc/MODIS_NPOESS/trunk/LEOCAT/data/VIIRS/LUT/CTP/COTlut.db'
	COTlut,nCOT, = readCOTfile(file)
	COTir = COTconvert(nCOT,COTlut,cps,cot)
	print "COTir = ",COTir

if __name__=='__main__':
	cot   = float(sys.argv[1])
	cps   = float(sys.argv[2])
	print "Input cloud optical thickness: ",cot
	print "Input cloud particle size: ",cps,"microns"

	main(cot,cps)
