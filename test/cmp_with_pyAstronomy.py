import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming import pyPDM

my1 = np.loadtxt("out1")
my2 = np.loadtxt("out2")

con = np.loadtxt("con_all.txt")
t0 = con[0, 0]
con[:, 0] -= t0

S = pyPDM.Scanner(minVal=1.0/(50.0*365.0), maxVal=1.0e-2+1.0/(50.0*365.0), dVal=1.0e-5, mode="frequency")

P = pyPDM.PyPDM(con[:, 0], con[:, 1])

f1, t1 = P.pdmEquiBinCover(10, 3, S)
# For comparison, carry out PDM analysis using 10 bins equidistant
# bins (no covers).
f2, t2 = P.pdmEquiBin(10, S)

plt.plot(1.0/f1, t1, color='k')
plt.plot(1.0/f2, t2, color='b')
plt.plot(my1[:, 0], my1[:, 1],color='g')
plt.plot(my2[:, 0], my2[:, 1],color='r')
plt.show()
