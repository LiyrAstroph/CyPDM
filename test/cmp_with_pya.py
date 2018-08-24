import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming import pyPDM
import sys;
sys.path.append("../")
import cypdm 

nbins = 10
covers = 3

con = np.loadtxt("con_all.txt")
t0 = con[0, 0]
con[:, 0] -= t0

S = pyPDM.Scanner(minVal=1.0/(50.0*365.0), maxVal=1.0e-2, dVal=1.0e-5, mode="frequency")
P = pyPDM.PyPDM(con[:, 0], con[:, 1])
f1, t1 = P.pdmEquiBinCover(nbins, covers, S)
f2, t2 = P.pdmEquiBin(nbins, S)

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], nbins, covers)
myp = 1.0/f1
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)

fig = plt.figure()
ax = fig.add_subplot(211)
plt.plot(1.0/f1, t1, color='k')
plt.plot(myp, myt1,color='g', marker='o', markersize=2)

ax = fig.add_subplot(212)
plt.plot(1.0/f2, t2, color='b')
plt.plot(myp, myt2,color='r', marker='o', markersize=2)

plt.show()
