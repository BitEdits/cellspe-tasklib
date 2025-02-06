#!/usr/bin/env python
import numpy
import pylab
import struct
import array
import math
import numarray
import sys
import string

if len(sys.argv) != 4:
	print "swsgif.py <file> <FFTpoints> <samplingfreq>"
	sys.exit()

f = open(sys.argv[1], 'r')
Nfft = string.atof(sys.argv[2])
fs = string.atof(sys.argv[3])
print f
print Nfft
print fs

Npts=math.floor(Nfft/2) + 1
data = array.array('f')	# 'f' for float
data.read(f, math.trunc(2*Npts))
spec = numpy.zeros([Npts], complex);
xvals = numpy.zeros([Npts], float);
for i in range(0, math.trunc(Npts-1)):
	spec[i] = complex(data[2*i+0], data[2*i+1])
	xvals[i] = i * (fs/Npts)

pylab.figure(figsize=(24, 12))
pylab.semilogy(xvals, abs(spec), 'b-') # poor because it "looses" the peaks at low print resolution
pylab.ylabel('abs(FFT)')
pylab.xlabel('Hertz')
pylab.grid(True)
pylab.savefig('testfig.png')
pylab.show()
