from math import asin
'''
SmithSpark Constant
'''

PI = 2 * asin(1.0)

# mesh
NR = 1001  # 点数
NR0 = 2
NR1 = NR - 1
R0 = 0
R1 = 10000
DR = (R1 - R0) / (NR - 1)  # dt=10   间距

# pump
MIUPUMP = 0.208
KPUMP = 0.184 * 10 ** -3
MPUMP = 3.9

# cytosolic spark
CCACYTREST = 1.0 * 10 ** -4
DCACYT = 3.5 * 10 ** 8
KRYR2 = 1.22522 * 10 ** 10

# Fluo-3 parameter
KPCAF = 80000
KMCAF = 90
BCAF = 0.05
KDCAF = KMCAF / KPCAF
DCAF = 2.0 * 10 ** 7

# Calmodulin parameter
KPCAM = 100000
KMCAM = 38
BCAM = 0.024
KDCAM = KMCAM / KPCAM

# Troponin C parameter
KPTRC = 39000
KMTRC = 20
BTRC = 0.07
KDTRC = KMTRC / KPTRC

# SR membrane parameter
KPSRM = 115000
KMSRM = 100
BSRM = 0.047
KDSRM = KMSRM / KPSRM

# SL membrane parameter
KPSLM = 115000
KMSLM = 1000
BSLM = 1.124
KDSLM = KMSLM / KPSLM