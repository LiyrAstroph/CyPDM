# CyPDM

A fast package to apply the phase disperion minimization (PDM) algorithm, based on PyAstronomy module PyPDM (https://github.com/sczesla/PyAstronomy).

The package is written in C and wrapped with Cython, so that it can be called in both C and Python code.

The PDM alogrithm refers to the reference:
http://adsabs.harvard.edu/abs/1978ApJ...224..953S

 
Author:
Yan-Rong Li, liyanrong@mail.ihep.ac.cn

# Compiling

```bash
python setup.py build_ext -i
```
This will generate a library ``cypdm.so``, which can be directly imported in Python.

# Usage

```python
import numpy as np
import cypdm 

con = np.loadtxt("con_all.txt", usecols=(0, 1))
t0 = con[0, 0]
con[:, 0] -= t0

pdm = cypdm.PyPDM(con[:, 0], con[:, 1], 10, 3)
myp = 1.0/np.linspace(1.0/50.0/365.0, 1.0e-2, 995, dtype=np.double)
myt1 = pdm.getPDM_EquiBinCover(myp)
myt2 = pdm.getPDM_EquiBin(myp)
```

This loads data from file "con_all.txt" and calculates PDMs.


