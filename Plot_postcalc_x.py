import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

plt.rcParams['font.sans-serif'] = ['simsun']
plt.rcParams['axes.unicode_minus'] = False

data = pd.read_csv("DATA_postcalc\\Fortran_dt_-8\\PROFILE_X\\PROFILE_X_01000000.csv")
x = np.array(data.iloc[:, 0])
ca = np.array(data.iloc[:, 1])
f = np.array(data.iloc[:, 2])
fpsf = np.array(data.iloc[:, 3])
cac = np.array(data.iloc[:, 4])
capsfc = np.array(data.iloc[:, 5])

plt.figure()
plt.grid()
plt.title("profile_x.csv  Fortran  dt=2*10**-8")
print(max(x))
plt.xlim((min(x),max(x)))
l1, = plt.plot(x, ca, linewidth=1)
l2, = plt.plot(x, f, linewidth=1)
l3, = plt.plot(x, fpsf, linewidth=1)
l4, = plt.plot(x, cac, linewidth=1)
l5, = plt.plot(x, capsfc, linewidth=1)
plt.legend(handles=[l1, l2, l3, l4, l5], labels=["Ca", "F", "FPSF", "CaC", "CaPSFC"])
plt.savefig("Figure\\PROFILE_X_Fortran.jpg")
plt.show()
