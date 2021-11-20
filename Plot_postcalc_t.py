import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.sans-serif'] = ['simsun']
plt.rcParams['axes.unicode_minus'] = False

data = pd.read_csv("DATA_postcalc\\Fortran_dt_-8\\PROFILE_T\\PROFILE_T.csv")
t = np.array(data.iloc[:, 0])
ca = np.array(data.iloc[:, 1])
f = np.array(data.iloc[:, 2])
fpsf = np.array(data.iloc[:, 3])
cac = np.array(data.iloc[:, 4])
capsfc = np.array(data.iloc[:, 5])

plt.figure()
plt.grid()
plt.title("profile_t.csv  Fortran  dt=2*10**-8")
plt.xlim((0, max(t)))
l1, = plt.plot(t, ca, linewidth=1)
l2, = plt.plot(t, f, linewidth=1)
l3, = plt.plot(t, fpsf, linewidth=1)
l4, = plt.plot(t, cac, linewidth=1)
l5, = plt.plot(t, capsfc, linewidth=1)
plt.legend(handles=[l1, l2, l3, l4, l5], labels=["Ca", "F", "FPSF", "CaC", "CaPSFC"])
plt.savefig("Figure\\PROFILE_T_Fortran.jpg")
plt.show()
