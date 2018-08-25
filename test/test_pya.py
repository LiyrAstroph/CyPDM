import sys
import numpy as np 
from PyAstronomy.pyTiming import pyPDM

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

nbins=10
covers=3

nd = 1000
P = 10.0
jd = np.linspace(0.0, 100.0, nd, dtype=np.double)
fs = np.sin(jd/P * 2.0*np.pi) + np.random.randn(nd)*0.1

S = pyPDM.Scanner(minVal=1.0/(50.0*365.0), maxVal=1.0e-2, dVal=1.0e-5, mode="frequency")
P = pyPDM.PyPDM(jd, fs)
f1, t1 = P.pdmEquiBinCover(nbins, covers, S)
f2, t2 = P.pdmEquiBin(nbins, S)
