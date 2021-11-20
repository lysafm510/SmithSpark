import os
from math import asin
import numpy as np


def spark_parameter():
    """
    spark参数
    """
    global KPCAF, KMCAF, BCAF, KDCAF
    global KPCAM, KMCAM, BCAM, KDCAM
    global KPTRC, KMTRC, BTRC, KDTRC
    global KPSRM, KMSRM, BSRM, KDSRM
    global KPSLM, KMSLM, BSLM, KDSLM
    global CCACYT, CCAF, CCAM, CTRC, CSRM, CSLM, NEWCCAF, DCAF
    global KRyR2, DCACYT
    global Jdye, Jbuffers, NEWCCACYT
    global END_TIME, RELEASE_TIME, DT, DR, PI
    global CCASTORE, MEDCCACYT, MEDCCAF

    # Fluo-3
    KPCAF = 80000  # K_plus_Ca_Fluno-3  结合速率
    KMCAF = 90  # K_minus_Ca_Fluno-3    解离速率
    BCAF = 0.05  # [F]T_Ca_Fluno-3     Total concentration for fluo-3
    KDCAF = KMCAF / KPCAF  # Kn = Kn-/Kn+
    DCAF = 2 * 10 ** 7  # 扩散系数

    # Calmodulin   钙调蛋白
    KPCAM = 100000
    KMCAM = 38
    BCAM = 0.024  # [Bn]T
    KDCAM = KMCAM / KPCAM

    # Troponin C   肌钙蛋白C
    KPTRC = 39000
    KMTRC = 20
    BTRC = 0.07
    KDTRC = KMTRC / KPTRC

    # SR membrane
    KPSRM = 115000
    KMSRM = 100
    BSRM = 0.047
    KDSRM = KMSRM / KPSRM

    # SL membrane
    KPSLM = 115000
    KMSLM = 1000
    BSLM = 1.124
    KDSLM = KMSLM / KPSLM

    KRyR2 = 1.22522 * 10 ** 10
    DCACYT = 3.5 * 10 ** 8  # 胞浆的Ca扩散系数
    CCASTORE = 1.0

    END_TIME = 0.1  # 结束时间
    RELEASE_TIME = 0.02  # 关闭RyR通道时间
    DT = 2 * 10 ** -8  # dt
    PI = 2 * asin(1.0)

    CCACYT = np.zeros(NR, float)
    NEWCCACYT = np.zeros(NR, float)
    MEDCCACYT = np.zeros(NR, float)
    CCAF = np.zeros(NR, float)
    NEWCCAF = np.zeros(NR, float)
    MEDCCAF = np.zeros(NR, float)
    CCAM = np.zeros(NR, float)
    CTRC = np.zeros(NR, float)
    CSRM = np.zeros(NR, float)
    CSLM = np.zeros(NR, float)
    Jdye = np.zeros(NR, float)
    Jbuffers = np.zeros(NR, float)


def load_spark_gridinfo():
    """
    钙火花网格
    """
    global NR, DR
    global RADIUS

    # 网格参数
    NR = 1001  # 1001个点，从0到10000
    DR = 10.0  # 间距
    RADIUS = np.zeros(NR)

    file = open("Grid\\SPARKMESH.dat", "w")
    # file.write(str(NR) + "\n")
    for i in range(0, NR):
        radius = i * DR
        RADIUS[i] = radius
        file.write(str(radius) + "\n")
    print("Spark网格生成")
    file.close()


def not_empty(s):
    return s and s.strip()


# **************************************************************************************
def cal_dye():
    """
    计算ca扩散方程的dye项,CCACYT，CCAF用n时刻的值
    """
    global Jdye
    for i in range(0, NR):
        Jdye[i] = -KPCAF * CCACYT[i] * (BCAF - CCAF[i]) + KMCAF * CCAF[i]


def cal_buffers():
    """
    计算ca扩散方程的buffers项，CCACYT，CCAM,CTRC,CSRM,CSLM用n时刻的值
    """
    global Jbuffers
    for i in range(0, NR):
        J1 = -KPCAM * CCACYT[i] * (BCAM - CCAM[i]) + KMCAM * CCAM[i]
        J2 = -KPTRC * CCACYT[i] * (BTRC - CTRC[i]) + KMTRC * CTRC[i]
        J3 = -KPSRM * CCACYT[i] * (BSRM - CSRM[i]) + KMSRM * CSRM[i]
        J4 = -KPSLM * CCACYT[i] * (BSLM - CSLM[i]) + KMSLM * CSLM[i]
        Jbuffers[i] = J1 + J2 + J3 + J4


def cytosolic_buffers_equation():
    '''
    更新[CaBi]
    '''
    global CCAM, CTRC, CSRM, CSLM
    for i in range(0, NR):
        # CCACYT用的上一步的值
        CCAM[i] = CCAM[i] + (KPCAM * CCACYT[i] * (BCAM - CCAM[i]) - KMCAM * CCAM[i]) * DT
        CTRC[i] = CTRC[i] + (KPTRC * CCACYT[i] * (BTRC - CTRC[i]) - KMTRC * CTRC[i]) * DT
        CSRM[i] = CSRM[i] + (KPSRM * CCACYT[i] * (BSRM - CSRM[i]) - KMSRM * CSRM[i]) * DT
        CSLM[i] = CSLM[i] + (KPSLM * CCACYT[i] * (BSLM - CSLM[i]) - KMSLM * CSLM[i]) * DT


def cytosolic_ca_equation():
    global RADIUS, fn_coeff, fn_const, CCACYT

    fn_coeff = np.zeros((NR, NR))
    fn_const = np.zeros(NR)

    coeff = 3 * DCACYT * DT / DR

    r3 = (RADIUS[0] + RADIUS[1]) / 2.0
    r1 = 0
    fn_coeff[0][0] = 1 + coeff * r3 * r3 / (r3 ** 3 - r1 ** 3) + 3 * DT * KRyR2 / (4 * PI * r3 ** 3)
    fn_coeff[0][1] = -coeff * r3 * r3 / (r3 ** 3 - r1 ** 3)
    fn_const[0] = CCACYT[0] + DT * Jdye[0] + DT * Jbuffers[0] + 3 * DT * KRyR2 * CCASTORE / (4 * PI * r3 ** 3)

    for i in range(1, NR - 1):
        r3 = (RADIUS[i] + RADIUS[i + 1]) / 2.0
        r1 = (RADIUS[i - 1] + RADIUS[i]) / 2.0
        fn_coeff[i][i - 1] = -coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
        fn_coeff[i][i] = 1 + coeff * (r3 * r3 + r1 * r1) / (r3 ** 3 - r1 ** 3)
        fn_coeff[i][i + 1] = -coeff * r3 * r3 / (r3 ** 3 - r1 ** 3)
        fn_const[i] = CCACYT[i] + DT * Jdye[i] + DT * Jbuffers[i]

    r3 = RADIUS[NR - 1]
    r1 = (RADIUS[NR - 1] + RADIUS[NR - 2]) / 2.0
    fn_coeff[NR - 1][NR - 2] = -coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
    fn_coeff[NR - 1][NR - 1] = 1 + coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
    fn_const[NR - 1] = CCACYT[NR - 1] + DT * Jdye[NR - 1] + DT * Jbuffers[NR - 1]

    CCACYT = np.linalg.solve(fn_coeff, fn_const)


def cytosolic_gn_equation():
    global RADIUS, gn_coeff, gn_const, CCAF
    gn_coeff = np.zeros((NR, NR))
    gn_const = np.zeros(NR)

    coeff = 3 * DCAF * DT / DR

    r3 = (RADIUS[0] + RADIUS[1]) / 2.0
    r1 = 0
    gn_coeff[0][0] = 1 + coeff * r3 * r3 / (r3 ** 3 - r1 ** 3)
    gn_coeff[0][1] = -coeff * r3 * r3 / (r3 ** 3 - r1 ** 3)
    gn_const[0] = CCAF[0] - DT * Jdye[0]

    for i in range(1, NR - 1):
        r3 = (RADIUS[i] + RADIUS[i + 1]) / 2.0
        r1 = (RADIUS[i - 1] + RADIUS[i]) / 2.0

        gn_coeff[i][i - 1] = -coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
        gn_coeff[i][i] = 1 + coeff * (r3 * r3 + r1 * r1) / (r3 ** 3 - r1 ** 3)
        gn_coeff[i][i + 1] = -coeff * r3 * r3 / (r3 ** 3 - r1 ** 3)
        gn_const[i] = CCAF[i] - DT * Jdye[i]

    r3 = RADIUS[NR - 1]
    r1 = (RADIUS[NR - 1] + RADIUS[NR - 2]) / 2.0

    gn_coeff[NR - 1][NR - 2] = -coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
    gn_coeff[NR - 1][NR - 1] = 1 + coeff * r1 * r1 / (r3 ** 3 - r1 ** 3)
    gn_const[NR - 1] = CCAF[NR - 1] - DT * Jdye[NR - 1]

    CCAF = np.linalg.solve(gn_coeff, gn_const)


# **************************************************************************

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


# *********************************  main  *********************************


def doloop_diffusion():
    global DCARYR, KRyR2
    global CCAMYO
    global AVG_CA_JSR, ICARYR, ICAFSR, AVG_Gn_JSR
    global current_step, times
    global CCACYT, CCAF

    START_STEP = 1  # 从多少步开始，默认从第1步开始
    SAVE_INTERVAL = 1  # 文件保存间隔
    savePath = "SPARK_solveEquations_dt_-8_saveInterval_1000"

    path = "DATA\\" + savePath + "\\SPARK"
    if START_STEP == 1:
        times = 0.0  # 当前时间
        current_step = 0
        # F0 初始值
        INICAF = 0.0001 * BCAF / (0.0001 + KDCAF)
        INICAM = 0.0001 * BCAM / (0.0001 + KDCAM)
        INITRC = 0.0001 * BTRC / (0.0001 + KDTRC)
        INISRM = 0.0001 * BSRM / (0.0001 + KDSRM)
        INISLM = 0.0001 * BSLM / (0.0001 + KDSLM)

        for i in range(0, NR):
            CCACYT[i] = 0.0001
            CCAF[i] = INICAF
            CCAM[i] = INICAM
            CTRC[i] = INITRC
            CSRM[i] = INISRM
            CSLM[i] = INISLM

        average_fn_gn()

        print("***** 数据开始写入 SPARK00000000.csv *****")
        folder = os.path.exists(path)
        if not folder:
            os.makedirs(path)
        file_s0 = open(path + "\\SPARK00000000.csv", "w")
        for i in range(0, NR):
            file_s0.write(str(CCACYT[i]) + "," + str(CCAF[i]) + "," + str(CCAM[i]) + "," + str(CTRC[i]) + "," + str(
                CSRM[i]) + "," + str(CSLM[i]) + "\n")
        file_s0.write(str(average_ccacyt) + "," + str(average_ccaf) + "\n")
        file_s0.write(str(current_step) + "," + str(times) + "\n")
        file_s0.close()

    else:  # 不是第一步开始，先读取先前保存的Spark
        file_name = "SPARK" + str(START_STEP - 1).zfill(8) + ".csv"
        read_spark = open("DATA\\SPARK_solveEquations_dt_-8_saveInterval_1\\SPARK\\" + file_name, "r")
        I = 0
        for line in read_spark.readlines():
            current_line = line.strip("\n").split(",")
            if I < NR:
                CCACYT[I] = float(current_line[0])
                CCAF[I] = float(current_line[1])
                CCAM[I] = float(current_line[2])
                CTRC[I] = float(current_line[3])
                CSRM[I] = float(current_line[4])
                CSLM[I] = float(current_line[5])
            else:
                # if I3 == NR:
                #     before_avg_ccacyt = float(current_line[0])
                #     before_avg_ccaf = float(current_line[1])
                if I == NR + 1:
                    current_step = int(current_line[0])
                    times = float(current_line[1])
            I = I + 1
        read_spark.close()
        print("读取上一步SPARK文件：", file_name)
        print()

    # if (abs(times - RELEASE_TIME) <= DT / 2):
    if (times >= (RELEASE_TIME - DT / 2)):
        # DCARYR = 0
        KRyR2 = 0
        print("关闭RyR通道")

    # for J in range(START_STEP, END_STEP + 1):
    while (times <= END_TIME):

        current_step = current_step + 1
        times = times + DT
        print("计算第" + str(current_step) + "步")

        cal_dye()
        cal_buffers()
        # 因为上面两个函数已经计算出fn和gn方程需要的dye和buffers
        # 所以下面计算缓冲物的函数对后面两个fn和gn方程没有影响
        # 式子放在fn和gn前面是因为用到的浓度是n时刻旧的值
        cytosolic_buffers_equation()
        cytosolic_ca_equation()
        cytosolic_gn_equation()
        average_fn_gn()

        # if (current_step % SAVE_INTERVAL == 0) or (current_step == END_STEP):
        if (current_step % SAVE_INTERVAL == 0) or (times == END_TIME):
            save_name = "SPARK" + str(str(current_step).zfill(8)) + ".csv"
            file_spark = open(path + "\\" + save_name, "w")
            print("***** 数据开始写入 " + save_name + " *****")
            for i in range(0, NR):
                file_spark.write(
                    str(CCACYT[i]) + "," + str(CCAF[i]) + "," + str(CCAM[i]) + "," + str(CTRC[i]) + "," + str(
                        CSRM[i]) + "," + str(CSLM[i]) + "\n")
            file_spark.write(str(average_ccacyt) + "," + str(average_ccaf) + "\n")
            file_spark.write(str(current_step) + ",")
            file_spark.write(str(times) + "\n")
            file_spark.close()
            print("***** 数据写入完毕 *****")

        if (abs(times - RELEASE_TIME) <= DT / 2):
            # DCARYR = 0
            KRyR2 = 0
            print("关闭RyR通道")


if __name__ == '__main__':
    load_spark_gridinfo()
    spark_parameter()
    doloop_diffusion()
