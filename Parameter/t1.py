PI = 3.14159265358979
NR = 1001
NR0 = 2
NR1 = NR - 1
NR2 = 201
R0 = 0
R1 = 10000
R2 = 2000
X0 = 0
X1 = 300
Y0 = 0
Y1 = PI / 4
NX = 151
NY = 51
NX0 = 2
NX1 = NX - 1
NY0 = 2
NY1 = NY - 1
NPSF = 300

# control parameter
DT = 2 * 10 ** -8
MAX_STEP = 0.4 + 1.0 * 10 ** -1 / DT
SAVE_INTERVAL = 10000
DR = (R1 - R0) / (NR - 1)
RELEASE_TIME = 2 * 10 ** -2
RELEASE_STEP = (0.4 + RELEASE_TIME / DT)

# fluo-3 parameter
KPCAF = 80000
KMCAF = 90
BCAF = 0.05
KDCAF = KMCAF / KPCAF
# cytosolic spark
CCAFREST = 0.0001 * BCAF / (0.0001 + KDCAF)
