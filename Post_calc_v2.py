from math import sqrt, log, exp, fmod, log10
from Parameter.t1 import PI, DR, NR, KDCAF, BCAF, NR2, CCAFREST, NPSF, MAX_STEP, SAVE_INTERVAL, RELEASE_STEP, DT
import numpy as np


def main():
    global step_time, n_step
    global CCACYT, CCAF, CCAM, CTRC, CSRM, CSLM
    global profile_t
    CCACYT = np.zeros(NR, float)
    CCAF = np.zeros(NR, float)
    CCAM = np.zeros(NR, float)
    CTRC = np.zeros(NR, float)
    CSRM = np.zeros(NR, float)
    CSLM = np.zeros(NR, float)

    print("计算psf")
    psf_parameter()
    print("计算完毕，写spark")

    # dut_spark.write(
    #     'zone i = ' + str(NR2 * 2 - 1) + ', j = ' + str(MAX_STEP / SAVE_INTERVAL + 1) + ', f = point' + "\n")

    # dut_blink = open("DATA_postcalc\\DUT_BLINK.csv", "w")
    # dut_blink.write('"X", "T", "Ca"' + "\n")

    # ca_current = open("DATA_postcalc\\CACURRENT.csv", "w")
    # ca_current.write('variables = "T", "Ica"' + '\n')
    # ca_current.write('zone I = ' + str(MAX_STEP / SAVE_INTERVAL + 1) + ', f = point')

    # profile_x.write('zone I = ' + str(NR2 * 2 - 1) + ', f = point')

    # profile_t.write('zone I = ' + str(MAX_STEP / SAVE_INTERVAL + 1) + ', f = point')
    profile_t = open("DATA_postcalc\\Fortran_dt_-8\\PROFILE_T\\PROFILE_T.csv", "w")
    profile_t.write('"T", "Ca", "F", "FPSF", "CaC", "CaPSFC"' + '\n')

    for n_step in range(0, int(MAX_STEP + 1)):
        if fmod(n_step, SAVE_INTERVAL) == 0:
            step_time = n_step * DT
            print(n_step)
            print(step_time)
            # 读取文件
            STRN = str(n_step).zfill(8)
            read_spark = open("DATA\\SPARK_Fortran_dt_-8_saveInterval_1000\\SPARK\\SPARK" + str(STRN) + ".csv", "r")
            i = 0
            for line in read_spark.readlines():
                current_line = line.strip("\n").split(",")
                if i < NR:
                    CCACYT[i] = float(current_line[0])
                    CCAF[i] = float(current_line[1])
                    CCAM[i] = float(current_line[2])
                    CTRC[i] = float(current_line[3])
                    CSRM[i] = float(current_line[4])
                    CSLM[i] = float(current_line[5])
                i = i + 1
            read_spark.close()

            dist = 0.0
            print("spark_space_distribution" + str(n_step))
            spark_space_distribution(dist)
            print("完毕")
    profile_t.close()


def spark_space_distribution(dist):
    c_ca_r = np.zeros(NR2)
    c_caf_r = np.zeros(NR2)
    c_caf_psf = np.zeros(NR2)
    c_cac = np.zeros(NR2)
    c_cac_psf = np.zeros(NR2)

    f_z_psf = np.zeros(NR)
    f_r_psf = np.zeros(NR)

    dist0 = (dist / DR) ** 2

    for i in range(0, NR2):
        ratio = sqrt(dist0 + i * i) + 1.0
        j = int(ratio)
        ratio = ratio - j
        c_ca_r[i] = CCACYT[j] * ratio + (1.0 - ratio) * CCACYT[j - 1]  # CCACYT
        c_caf_r[i] = CCAF[j] * ratio + (1.0 - ratio) * CCAF[j - 1]  # CCAF
        c_cac[i] = KDCAF * c_caf_r[i] / (BCAF - c_caf_r[i])

    for i in range(0, NR):
        f_z_psf[i] = 0
        for j in range(1, NPSF):
            ratio = sqrt(1.0 * i * i + j * j) + 1.0
            k = int(ratio)
            if k < NR:
                ratio = ratio - k
                f_psf_temp = CCAF[k] * ratio + (1.0 - ratio) * CCAF[k - 1]
                f_z_psf[i] = f_z_psf[i] + z_psf[j] * f_psf_temp
            else:
                f_z_psf[i] = f_z_psf[i] + z_psf[j] * CCAFREST
        f_z_psf[i] = f_z_psf[i] * 2 + z_psf[0] * CCAF[i]

    dist0 = dist / DR
    for i in range(0, NR):
        ratio = sqrt(dist0 ** 2 + i * i) + 1.0
        k = int(ratio)
        if (k < NR):
            ratio = ratio - k
            f_psf_temp = ratio * f_z_psf[k] + (1.0 - ratio) * f_z_psf[k - 1]
            f_r_psf[i] = xy_psf[0] * f_psf_temp
        else:
            f_r_psf[i] = xy_psf[1] * CCAFREST

        for j in range(1, NPSF):
            ratio = sqrt((dist0 - j) ** 2 + i * i) + 1.0
            k = int(ratio)
            if k < NR:
                ratio = ratio - k
                f_psf_temp = ratio * f_z_psf[k] + (1.0 - ratio) * f_z_psf[k - 1]
                f_r_psf[i] = f_r_psf[i] + xy_psf[j] * f_psf_temp
            else:
                f_r_psf[i] = f_r_psf[i] + xy_psf[j] * CCAFREST

            ratio = sqrt((dist0 + j) ** 2 + i * i) + 1.0
            k = int(ratio)
            if k < NR:
                ratio = ratio - k
                f_psf_temp = ratio * f_z_psf[k] + (1.0 - ratio) * f_z_psf[k - 1]
                f_r_psf[i] = f_r_psf[i] + xy_psf[j] * f_psf_temp
            else:
                f_r_psf[i] = f_r_psf[i] + xy_psf[j] * CCAFREST

    for i in range(0, NR2):
        c_caf_psf[i] = f_r_psf[i] * xy_psf[0]
        for j in range(1, NPSF):
            k = i + j + 1
            l = i - j + 1
            if k <= NR:
                c_caf_psf[i] = c_caf_psf[i] + xy_psf[j] * f_r_psf[k - 1]
            else:
                c_caf_psf[i] = c_caf_psf[i] + xy_psf[j] * CCAFREST
            if l <= 0:
                l = 2 - l
            c_caf_psf[i] = c_caf_psf[i] + xy_psf[j] * f_r_psf[l - 1]
        c_cac_psf[i] = KDCAF * c_caf_psf[i] / (BCAF - c_caf_psf[i])

    print("DUT_SPARK")
    dut_spark = open("DATA_postcalc\\Fortran_dt_-8\\DUT_SPARK\\DUT_SPARK_" + str(n_step).zfill(8) + ".csv", "w")
    dut_spark.write('"X", "Y", "Ca", "F", "FPSF", "CaC", "CaPSFC"' + '\n')
    for i in range(NR2 - 1, 0, -1):
        dut_spark.write(str(-i * DR) + ',' + str(step_time) + ',' + str(c_ca_r[i]) + ',' + str(c_caf_r[i]) + ',' + str(
            c_caf_psf[i]) + ',' + str(c_cac[i]) + ',' + str(c_cac_psf[i]) + '\n')
    for i in range(0, NR2):
        dut_spark.write(str(i * DR) + ',' + str(step_time) + ',' + str(c_ca_r[i]) + ',' + str(c_caf_r[i]) + ',' + str(
            c_caf_psf[i]) + ',' + str(c_cac[i]) + ',' + str(c_cac_psf[i]) + '\n')
    dut_spark.write(str(n_step) + "\n")
    dut_spark.close()

    if (n_step == RELEASE_STEP - 0.4):
        print("PROFILE_X")
        profile_x = open("DATA_postcalc\\Fortran_dt_-8\\PROFILE_X\\PROFILE_X_" + str(n_step).zfill(8) + ".csv", "w")
        profile_x.write('"X", "Ca", "F", "FPSF", "CaC", "CaPSFC"' + '\n')
        for i in range(NR2 - 1, 0, -1):
            profile_x.write(str(-i * DR) + ',' + str(log10(c_ca_r[i])) + ',' + str(
                log10(c_caf_r[i])) + ',' + str(
                log10(c_caf_psf[i])) + ',' + str(log10(c_cac[i])) + ',' + str(log10(c_cac_psf[i])) + '\n')
        for i in range(0, NR2):
            profile_x.write(str(i * DR) + ',' + str(log10(c_ca_r[i])) + ',' + str(
                log10(c_caf_r[i])) + ',' + str(
                log10(c_caf_psf[i])) + ',' + str(log10(c_cac[i])) + ',' + str(log10(c_cac_psf[i])) + '\n')
        # profile_x.write(str(n_step) + "\n")
        profile_x.close()

    print("PROFILE_T")
    profile_t.write(str(step_time) + ',' + str(log10(c_ca_r[0])) + ',' + str(
        log10(c_caf_r[0])) + ',' + str(
        log10(c_caf_psf[0])) + ',' + str(log10(c_cac[0])) + ',' + str(log10(c_cac_psf[0])) + '\n')


def psf_parameter():
    global xy_psf, z_psf
    xy_psf = np.zeros(NPSF)
    z_psf = np.zeros(NPSF)

    coe_f = 1.0 / (sqrt(0.25 * PI / log(2.0)) * (400 / DR))
    sigma = (400 / DR) ** 2 / log(2.0) * 0.25
    for i in range(0, NPSF):
        xy_psf[i] = coe_f * exp(-i ** 2 / sigma)

    coe_f = 1.0 / (sqrt(0.25 * PI / log(2.0)) * (800 / DR))
    sigma = (800 / DR) ** 2 / log(2.0) * 0.25
    for i in range(0, NPSF):
        z_psf[i] = coe_f * exp(-i ** 2 / sigma)


if __name__ == '__main__':
    main()
