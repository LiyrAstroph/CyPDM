import numpy as np
import matplotlib.pyplot as plt
from PyAstronomy.pyTiming import pyPDM
import sys
sys.path.append("../")
import cypdm 

plt.rc('text', usetex=True)
plt.rc('font', family='serif', size=12)

nbins = 5
covers = 3

nd = 1000
P = 10.0
jd = np.linspace(0.0, 100.0, nd, dtype=np.double)
fs = np.sin(jd/P * 2.0*np.pi) + np.random.randn(nd)*0.1


S = pyPDM.Scanner(minVal=1.0/(50.0), maxVal=1.0e0, dVal=1.0e-4, mode="frequency")
P = pyPDM.PyPDM(jd, fs)
f1, t1 = P.pdmEquiBinCover(nbins, covers, S)
f2, t2 = P.pdmEquiBin(nbins, S)

pdm = cypdm.PyPDM(jd, fs, nbins, covers)
myp = 1.0/f1
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)

fig = plt.figure(figsize=(8, 9))
ax = fig.add_subplot(311)
ax.plot(jd, fs)
ax.set_xlabel(r"Time")
ax.set_ylabel(r"Flux")

ax = fig.add_subplot(312)
plt.plot(1.0/f1, t1, color='k', label="PyAstronomy")
plt.plot(myp, myt1,color='g', marker='o', markersize=2, label="CyPDM")

ax.legend()
ax.set_xlabel(r"Period")
ax.set_ylabel(r"$\Theta_{\rm PDM}$")

ax.text(40.0, 0.45, r"EquiBinCover")

ax = fig.add_subplot(313)
plt.plot(1.0/f2, t2, color='k', label="PyAstronomy")
plt.plot(myp, myt2,color='g', marker='o', markersize=2, label="CyPDM")

ax.legend()

ax.set_xlabel(r"Period")
ax.set_ylabel(r"$\Theta_{\rm PDM}$")

ax.text(40.0, 0.45, r"EquiBin")


fig.savefig("cypdm_cmp.jpg", bbox_inches='tight')
plt.show()
