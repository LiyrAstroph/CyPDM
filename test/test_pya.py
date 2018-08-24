import sys
import numpy as np 
from PyAstronomy.pyTiming import pyPDM

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

nbins=10
covers=3

S = pyPDM.Scanner(minVal=1.0/(50.0*365.0), maxVal=1.0e-2, dVal=1.0e-5, mode="frequency")
P = pyPDM.PyPDM(con[:, 0], con[:, 1])
f1, t1 = P.pdmEquiBinCover(nbins, covers, S)
f2, t2 = P.pdmEquiBin(nbins, S)
