import os
from math import asin

import numpy as np

path = "SPARK_Fortran_dt_-8_saveInterval_1000"


def plot():
    os.chdir("DATA\\" + path + "\\SPARK")
    for file in os.listdir():
        if '.csv' in file:
            data = np.loadtxt(file)
            c_ca_cyt = data[:, 0]
            print(c_ca_cyt)


def average(concentration):
    '''
    加权体积求平均值
    '''
    PI = 2 * asin(1.0)
    NR = 1001
    DR = 10
    RADIUS = np.zeros(NR)
    for i in range(0, NR):
        radius = i * DR
        RADIUS[i] = radius

    total_Concentration = 0.0
    total_V = (4 * PI * RADIUS[NR - 1] ** 3) / 3.0
    for i in range(0, NR):
        if i == 0:
            r3 = (RADIUS[0] + RADIUS[1]) / 2.0
            r1 = 0
        elif i == NR - 1:
            r3 = RADIUS[NR - 1]
            r1 = (RADIUS[NR - 1] + RADIUS[NR - 2]) / 2.0
        else:
            r3 = (RADIUS[i] + RADIUS[i + 1]) / 2.0
            r1 = (RADIUS[i - 1] + RADIUS[i]) / 2.0

        ctrl_V = (4 * PI / 3) * (r3 ** 3 - r1 ** 3)
        total_Concentration = total_Concentration + concentration[i] * ctrl_V
    average_Concentration = total_Concentration / total_V
    return average_Concentration

plot()