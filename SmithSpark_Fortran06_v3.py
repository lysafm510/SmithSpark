from math import fmod
import numpy as np
import os
from Parameter.t0 import PI, NR, NR0, NR1, DR, MIUPUMP, KPUMP, MPUMP, CCACYTREST, DCACYT, KRYR2, KPCAF, KMCAF, BCAF, \
    KDCAF, DCAF, KPCAM, KMCAM, BCAM, KDCAM, KPTRC, KMTRC, BTRC, KDTRC, KPSRM, KMSRM, BSRM, KDSRM, KPSLM, KMSLM, BSLM, \
    KDSLM


def parameter():
    global CCASTORE
    global DT, RADIUS, RELEASE_TIME, STOP_TIME, SAVE_INTERVAL, ISTART, savePath
    global CCACYT, MEDCACYT, NEWCACYT, COESPK1, COESPK2
    global CCAF, MEDCAF, NEWCAF, COECAF1, COECAF2
    global CCAM, MEDCAM, NEWCAM, CTRC, MEDTRC, NEWTRC, CSRM, MEDSRM, NEWSRM, CSLM, MEDSLM, NEWSLM

    CCASTORE = 1.0

    # control parameter
    DT = 2 * 10 ** -6
    RELEASE_TIME = 2.0 * 10 ** -2
    STOP_TIME = 0.1
    SAVE_INTERVAL = 1  # 保存间隔
    ISTART = 0  # 起始步数，默认从0开始
    savePath = 'SPARK_Fortran_dt_-6_saveInterval_1'

    # 初始化
    # mesh
    RADIUS = np.zeros(NR)
    for i in range(0, NR):
        radius = i * DR
        RADIUS[i] = radius
    # cytosolic spark
    CCACYT = np.zeros(NR, float)
    MEDCACYT = np.zeros(NR, float)
    NEWCACYT = np.zeros(NR, float)
    COESPK1 = np.zeros(NR, float)
    COESPK2 = np.zeros(NR, float)
    # Fluo-3 parameter
    CCAF = np.zeros(NR, float)
    MEDCAF = np.zeros(NR, float)
    NEWCAF = np.zeros(NR, float)
    COECAF1 = np.zeros(NR, float)
    COECAF2 = np.zeros(NR, float)
    # Calmodulin parameter
    CCAM = np.zeros(NR, float)
    MEDCAM = np.zeros(NR, float)
    NEWCAM = np.zeros(NR, float)
    # Troponin C parameter
    CTRC = np.zeros(NR, float)
    MEDTRC = np.zeros(NR, float)
    NEWTRC = np.zeros(NR, float)
    # SR membrane parameter
    CSRM = np.zeros(NR, float)
    MEDSRM = np.zeros(NR, float)
    NEWSRM = np.zeros(NR, float)
    # SL membrane parameter
    CSLM = np.zeros(NR, float)
    MEDSLM = np.zeros(NR, float)
    NEWSLM = np.zeros(NR, float)


# *****************************    Fortran翻译    *****************************
def cytosolic_ca_equation():
    global TMPCOEA, TMPCOEB, TMPCOEC, TMPCOEV
    global NEWCACYT, NEWCAM, NEWCAF, NEWTRC, NEWSRM, NEWSL
    global MEDCACYT, MEDCAM, MEDCAF, MEDTRC, MEDSRM, MEDSLM
    TMPCOEA = np.zeros(NR, float)
    TMPCOEB = np.zeros(NR, float)
    TMPCOEC = np.zeros(NR, float)
    TMPCOEV = np.zeros(NR, float)

    # pre-estimate
    Jleak = - MIUPUMP / (1 + (KPUMP / CCACYTREST) ** MPUMP)

    TMPCOEV[0] = (KMCAF * CCAF[0] + KMCAM * CCAM[0] + KMTRC * CTRC[0] + KMSRM * CSRM[0] + KMSLM * CSLM[0] + MIUPUMP / (
            1 + (KPUMP / CCACYT[0]) ** MPUMP) + Jleak) * DT / 2 + CCACYT[0] + KRYR2 * CCASTORE * 3 * DT / (
                         PI * DR * DR * DR)
    TMPCOEA[0] = (KPCAF * (BCAF - CCAF[0]) + KPCAM * (BCAM - CCAM[0]) + KPTRC * (BTRC - CTRC[0]) + KPSRM * (
            BSRM - CSRM[0]) + KPSLM * (BSLM - CSLM[0])) * DT / 2 + COESPK2[0] + KRYR2 * 3 * DT / (
                         PI * DR * DR * DR) + 1.0
    TMPCOEB[0] = -COESPK2[0]
    TMPCOEC[0] = -COESPK1[0]
    for I in range(1, NR):
        TMPCOEV[I] = (KMCAF * CCAF[I] + KMCAM * CCAM[I] + KMTRC * CTRC[I] + KMSRM * CSRM[I] + KMSLM * CSLM[
            I] + MIUPUMP / (1 + (KPUMP / CCACYT[I]) ** MPUMP) + Jleak) * DT / 2 + CCACYT[I]
        TMPCOEA[I] = (KPCAF * (BCAF - CCAF[I]) + KPCAM * (BCAM - CCAM[I]) + KPTRC * (BTRC - CTRC[I]) + KPSRM * (
                BSRM - CSRM[I]) + KPSLM * (BSLM - CSLM[I])) * DT / 2 + COESPK1[I] + COESPK2[I] + 1.0
        TMPCOEB[I] = -COESPK2[I]
        TMPCOEC[I] = -COESPK1[I]

    trdiag(TMPCOEA, TMPCOEB, TMPCOEC, MEDCACYT, TMPCOEV, NR)

    for I in range(0, NR):
        TMPCOEV[I] = KPCAF * (BCAF - CCAF[I]) * MEDCACYT[I] - KMCAF * CCAF[I]
        TMPCOEV[I] = TMPCOEV[I] * DT * 0.5 + CCAF[I]
        TMPCOEA[I] = COECAF1[I] + COECAF2[I] + 1.0
        TMPCOEB[I] = -COECAF2[I]
        TMPCOEC[I] = -COECAF1[I]

    trdiag(TMPCOEA, TMPCOEB, TMPCOEC, MEDCAF, TMPCOEV, NR)

    for I in range(0, NR):
        MEDCAM[I] = CCAM[I] + (KPCAM * (BCAM - CCAM[I]) * MEDCACYT[I] - KMCAM * CCAM[I]) * DT * 0.5
        MEDTRC[I] = CTRC[I] + (KPTRC * (BTRC - CTRC[I]) * MEDCACYT[I] - KMTRC * CTRC[I]) * DT * 0.5
        MEDSRM[I] = CSRM[I] + (KPSRM * (BSRM - CSRM[I]) * MEDCACYT[I] - KMSRM * CSRM[I]) * DT * 0.5
        MEDSLM[I] = CSLM[I] + (KPSLM * (BSLM - CSLM[I]) * MEDCACYT[I] - KMSLM * CSLM[I]) * DT * 0.5

    # correct
    TMPCOEA[0] = (KPCAF * (BCAF - MEDCAF[0]) + KPCAM * (BCAM - MEDCAM[0]) + KPTRC * (BTRC - MEDTRC[0]) + KPSRM * (
            BSRM - MEDSRM[0]) + KPSLM * (BSLM - MEDSLM[0])) * DT / 2 + COESPK2[0] + KRYR2 * 3 * DT / (
                         PI * DR * DR * DR)
    TMPCOEV[0] = (KMCAF * MEDCAF[0] + KMCAM * MEDCAM[0] + KMTRC * MEDTRC[0] + KMSRM * MEDSRM[0] + KMSLM * MEDSLM[
        0] + MIUPUMP / (1 + (KPUMP / MEDCACYT[0]) ** MPUMP) + Jleak) * DT + (1 - TMPCOEA[0]) * CCACYT[0] + COESPK2[0] * \
                 CCACYT[1] + KRYR2 * CCASTORE * 6 * DT / (PI * DR * DR * DR)
    TMPCOEA[0] = TMPCOEA[0] + 1.0
    TMPCOEB[0] = -COESPK2[0]
    TMPCOEC[0] = -COESPK1[0]
    for I in range(1, NR1):
        TMPCOEA[I] = (KPCAF * (BCAF - MEDCAF[I]) + KPCAM * (BCAM - MEDCAM[I]) + KPTRC * (BTRC - MEDTRC[I]) + KPSRM * (
                BSRM - MEDSRM[I]) + KPSLM * (BSLM - MEDSLM[I])) * DT * 0.5 + COESPK2[I] + COESPK1[I]
        TMPCOEV[I] = (KMCAF * MEDCAF[I] + KMCAM * MEDCAM[I] + KMTRC * MEDTRC[I] + KMSRM * MEDSRM[I] + KMSLM * MEDSLM[
            I] + MIUPUMP / (1 + (KPUMP / MEDCACYT[I]) ** MPUMP) + Jleak) * DT + (1 - TMPCOEA[I]) * CCACYT[I] + COESPK2[
                         I] * CCACYT[I + 1] + COESPK1[I] * CCACYT[I - 1]
        TMPCOEA[I] = TMPCOEA[I] + 1.0
        TMPCOEB[I] = -COESPK2[I]
        TMPCOEC[I] = -COESPK1[I]
    I = NR - 1
    TMPCOEA[I] = (KPCAF * (BCAF - MEDCAF[I]) + KPCAM * (BCAM - MEDCAM[I]) + KPTRC * (BTRC - MEDTRC[I]) + KPSRM * (
            BSRM - MEDSRM[I]) + KPSLM * (BSLM - MEDSLM[I])) * DT * 0.5 + COESPK2[I] + COESPK1[I]
    TMPCOEV[I] = (KMCAF * MEDCAF[I] + KMCAM * MEDCAM[I] + KMTRC * MEDTRC[I] + KMSRM * MEDSRM[I] + KMSLM * MEDSLM[
        I] + MIUPUMP / (1 + (KPUMP / MEDCACYT[I]) ** MPUMP) + Jleak) * DT + (1 - TMPCOEA[I]) * CCACYT[I] + COESPK1[I] * \
                 CCACYT[I - 1]
    TMPCOEA[I] = TMPCOEA[I] + 1.0
    TMPCOEB[I] = -COESPK2[I]
    TMPCOEC[I] = -COESPK1[I]

    trdiag(TMPCOEA, TMPCOEB, TMPCOEC, NEWCACYT, TMPCOEV, NR)

    I = 0
    TMPCOEV[I] = (0.5 * KPCAF * (BCAF - MEDCAF[I]) * (CCACYT[I] + NEWCACYT[I]) - KMCAF * MEDCAF[I]) * DT + (
            1.0 - COECAF1[I] - COECAF2[I]) * CCAF[I] + CCAF[I + 1] * COECAF2[I]
    TMPCOEA[I] = COECAF1[I] + COECAF2[I] + 1.0
    TMPCOEB[I] = -COECAF2[I]
    TMPCOEC[I] = -COECAF1[I]
    for I in range(1, NR1):
        TMPCOEV[I] = (0.5 * KPCAF * (BCAF - MEDCAF[I]) * (CCACYT[I] + NEWCACYT[I]) - KMCAF * MEDCAF[I]) * DT + (
                1.0 - COECAF1[I] - COECAF2[I]) * CCAF[I] + CCAF[I - 1] * COECAF1[I] + CCAF[I + 1] * COECAF2[I]
        TMPCOEA[I] = COECAF1[I] + COECAF2[I] + 1.0
        TMPCOEB[I] = -COECAF2[I]
        TMPCOEC[I] = -COECAF1[I]
    I = NR - 1
    TMPCOEV[I] = (0.5 * KPCAF * (BCAF - MEDCAF[I]) * (CCACYT[I] + NEWCACYT[I]) - KMCAF * MEDCAF[I]) * DT + (
            1.0 - COECAF1[I] - COECAF2[I]) * CCAF[I] + CCAF[I - 1] * COECAF1[I]
    TMPCOEA[I] = COECAF1[I] + COECAF2[I] + 1.0
    TMPCOEB[I] = -COECAF2[I]
    TMPCOEC[I] = -COECAF1[I]

    trdiag(TMPCOEA, TMPCOEB, TMPCOEC, NEWCAF, TMPCOEV, NR)

    for I in range(0, NR):
        NEWCAM[I] = CCAM[I] + (KPCAM * (BCAM - MEDCAM[I]) * MEDCACYT[I] - KMCAM * MEDCAM[I]) * DT * 0.5
        NEWTRC[I] = CTRC[I] + (KPTRC * (BTRC - MEDTRC[I]) * MEDCACYT[I] - KMTRC * MEDTRC[I]) * DT * 0.5
        NEWSRM[I] = CSRM[I] + (KPSRM * (BSRM - MEDSRM[I]) * MEDCACYT[I] - KMSRM * MEDSRM[I]) * DT * 0.5
        NEWSLM[I] = CSLM[I] + (KPSLM * (BSLM - MEDSLM[I]) * MEDCACYT[I] - KMSLM * MEDSLM[I]) * DT * 0.5


def trdiag(A, B, C, X, F, N):
    B[0] = B[0] / A[0]
    F[0] = F[0] / A[0]
    for I in range(1, N):
        A[I] = A[I] - C[I] * B[I - 1]
        B[I] = B[I] / A[I]
        F[I] = (F[I] - F[I - 1] * C[I]) / A[I]

    X[N - 1] = F[N - 1]
    for I in range(N - 2, -1, -1):
        X[I] = F[I] - B[I] * X[I + 1]


def spark_coefficient_matrix():
    global COESPK1, COESPK2, COECAF1, COECAF2
    COESPK1 = np.zeros(NR, float)
    COESPK2 = np.zeros(NR, float)
    COECAF1 = np.zeros(NR, float)
    COECAF2 = np.zeros(NR, float)

    const1 = 1.5 * DCACYT * DT / (DR * DR)
    const2 = 1.5 * DCAF * DT / (DR * DR)
    COESPK1[0] = 0
    COESPK2[0] = const1 * 2
    COECAF1[0] = 0
    COECAF2[0] = const2 * 2

    for I in range(NR0, NR + 1):
        COESPK1[I - 1] = const1 * (I - 1.5) * (I - 1.5) / (3.25 + 3 * I * (I - 2))
        COESPK2[I - 1] = const1 * (I - 0.5) * (I - 0.5) / (3.25 + 3 * I * (I - 2))
        COECAF1[I - 1] = const2 * (I - 1.5) * (I - 1.5) / (3.25 + 3 * I * (I - 2))
        COECAF2[I - 1] = const2 * (I - 0.5) * (I - 0.5) / (3.25 + 3 * I * (I - 2))

    COESPK1[NR - 1] = COESPK1[NR - 1] + COESPK2[NR - 1]
    COESPK2[NR - 1] = 0
    COECAF1[NR - 1] = COECAF1[NR - 1] + COECAF2[NR - 1]
    COECAF2[NR - 1] = 0


# *****************************************************************
def doloop_diffusion():
    global times, current_step
    global CCACYT, CCAF, CCAM, CTRC, CSRM, CSLM
    global RADIUS

    # #生成网格
    # file = open("DATA\\SPARKMESH.dat", "w")
    # file.write(str(NR) + "\n")
    # for I in range(0, NR):
    #     radius = R0 + I * DR
    #     file.write(str(radius) + "\n")
    # file.close()
    path = "DATA\\" + savePath + "\\SPARK"
    if ISTART == 0:
        times = 0.0  # 当前时间
        current_step = 0  # 当前步数

        INICAF = 0.0001 * BCAF / (0.0001 + KDCAF)
        INICAM = 0.0001 * BCAM / (0.0001 + KDCAM)
        INITRC = 0.0001 * BTRC / (0.0001 + KDTRC)
        INISRM = 0.0001 * BSRM / (0.0001 + KDSRM)
        INISLM = 0.0001 * BSLM / (0.0001 + KDSLM)

        for I in range(0, NR):
            CCACYT[I] = 0.0001
            CCAF[I] = INICAF
            CCAM[I] = INICAM
            CTRC[I] = INITRC
            CSRM[I] = INISRM
            CSLM[I] = INISLM

        average_fn_gn()

        print("***** 数据开始写入 SPARK00000000.csv *****")
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        file = open(path + "\\SPARK00000000.csv", "w")
        for I in range(0, NR):
            file.write(str(CCACYT[I]) + "," + str(CCAF[I]) + "," + str(CCAM[I]) + "," + str(CTRC[I]) + "," + str(
                CSRM[I]) + "," + str(CSLM[I]) + "\n")
        file.write(str(average_ccacyt) + "," + str(average_ccaf) + "\n")
        file.write(str(current_step) + "," + str(times))
        file.close()
    else:  # 不是第一步开始，先读取先前保存的Spark
        file_name = "SPARK" + str(ISTART - 1).zfill(8) + ".csv"
        file = open(path + "\\" + file_name, "r")
        I = 0
        for line in file.readlines():
            current_line = line.strip("\n").split(",")
            if I < NR:
                CCACYT[I] = float(current_line[0])
                CCAF[I] = float(current_line[1])
                CCAM[I] = float(current_line[2])
                CTRC[I] = float(current_line[3])
                CSRM[I] = float(current_line[4])
                CSLM[I] = float(current_line[5])
            else:
                if I == NR + 1:
                    current_step = int(current_line[0])
                    times = float(current_line[1])
            I = I + 1
        file.close()
        print("读取上一步SPARK文件：", file_name)
        print()

    if (times >= (RELEASE_TIME - DT / 2)):
        KRYR2 = 0
        print("关闭RyR通道")

    spark_coefficient_matrix()

    while (times <= STOP_TIME):

        current_step = current_step + 1
        times = times + DT
        print("计算第" + str(current_step) + "步")

        cytosolic_ca_equation()
        average_fn_gn()
        for I in range(0, NR):
            CCACYT[I] = NEWCACYT[I]
            CCAF[I] = NEWCAF[I]
            CCAM[I] = NEWCAM[I]
            CTRC[I] = NEWTRC[I]
            CSRM[I] = NEWSRM[I]
            CSLM[I] = NEWSLM[I]

        if (fmod(current_step, SAVE_INTERVAL) == 0) or (times == STOP_TIME):

            save_name = "SPARK" + str(str(current_step).zfill(8)) + ".csv"
            file = open(path + "\\" + save_name, "w")
            print("***** 数据开始写入 " + save_name + " *****")
            for I in range(0, NR):
                file.write(str(CCACYT[I]) + "," + str(CCAF[I]) + "," + str(CCAM[I]) + "," + str(CTRC[I]) + "," + str(
                    CSRM[I]) + "," + str(CSLM[I]) + "\n")
            file.write(str(average_ccacyt) + "," + str(average_ccaf) + "\n")
            file.write(str(current_step) + "," + str(times))
            file.close()

        if (abs(times - RELEASE_TIME) <= DT / 2):
            KRYR2 = 0
            print("关闭RyR通道")


def average_fn_gn():
    '''
    加权体积求平均值
    '''
    global RADIUS, CCACYT, CCAF, average_ccacyt, average_ccaf

    total_f = 0.0
    total_g = 0.0
    average_ccacyt = 0.0
    average_ccaf = 0.0
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
        total_f = total_f + CCACYT[i] * ctrl_V
        total_g = total_g + CCAF[i] * ctrl_V
    average_ccacyt = total_f / total_V
    average_ccaf = total_g / total_V


if __name__ == '__main__':
    parameter()
    doloop_diffusion()
