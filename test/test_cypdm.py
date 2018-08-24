import sys
sys.path.append("../")
import cypdm 
import numpy as np 
import matplotlib.pyplot as plt 

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], -10)
myp = 1.0/np.linspace(1.0/50.0/365.0, 1.0e-2, 995, dtype=np.double)
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)


