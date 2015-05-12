#!/bin/env python
import manode
import pylab as pl
import numpy as np
import xraylib_np as xrl
z = np.array([64],dtype=np.int)
vtube = 100
itube = 6
estep = 0.1
nstep = int(vtube/estep)
x = np.arange(0.0001,100,estep,dtype=np.double)
y = np.zeros(x.shape[0])
manode.anode.pydata(vtube,itube,estep)
for i in range(0,x.shape[0]):
	a = manode.anode.anode_int(z, np.take(x,i))
	np.put(y, i, a)
x = np.reshape(x, x.shape[0])
y = np.reshape(y, x.shape[0])
print x
print y
pl.plot(x, y)
pl.show()
