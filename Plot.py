import os
import linecache
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import csv
plt.rcParams["font.family"] = "Microsoft Yahei"
plt.rcParams['font.sans-serif'] = ['SimHei']
plt.rcParams['axes.unicode_minus'] = False
'''
Spark 可视化曲线图
'''

path = "DATA\\SPARK_Iteration_dt_-8_saveInterval_1000"
save_interval = 1000
file = "SPARK00010001.csv"


# ******************************* 画平均值曲线 ****************************
def exact_average():
    """
    提取每个文件平均值
    """
    fn = open(path + "\\AVG_CCACYT.csv", "w")
    gn = open(path + "\\AVG_CCAF.csv", "w")
    for filename in os.listdir(path + "\\SPARK"):
        print("读取" + filename)
        line = linecache.getline(path + "\\SPARK\\" + filename, 1002)
        data = line.strip("\n").split(",")
        fn.write(str(data[0]) + "\n")
        gn.write(str(data[1]) + "\n")
        print(filename)
    fn.close()
    gn.close()


def plot_average():
    """
    绘制平均值曲线图
    """
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    i = 0
    j = 0
    with open(path + "\\AVG_CCACYT.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x1.append(i * save_interval)
            y1.append(float(row[0]))
            i = i + 1
    plt.figure()
    plt.grid()
    plt.plot(x1, y1)
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.8f'))
    plt.gca().xaxis.get_major_formatter().set_scientific(False)
    plt.title("Fn")
    # plt.xlim((0,i*save_interval))
    plt.savefig(path + "\\avg_fn.jpg", bbox_inches='tight')
    plt.show()

    with open(path + "\\AVG_CCAF.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x2.append(j * save_interval)
            y2.append(float(row[0]))
            j = j + 1
    plt.figure()
    plt.grid()
    plt.plot(x2, y2)
    # plt.xlim((0, j * save_interval))
    plt.title("Gn")
    plt.gca().xaxis.get_major_formatter().set_scientific(False)
    plt.savefig(path + "\\avg_gn.jpg", bbox_inches='tight')
    plt.show()


def plot_average_multiple():
    x1 = []
    y1 = []
    x2 = []
    y2 = []

    i = 0
    j = 0
    path1 = "DATA\\SPARK_Iteration_dt_-8_saveInterval_1000"
    path2 = "DATA\\SPARK_Fortran_dt_-8_saveInterval_1000"
    save_interval = 1000
    with open(path1 + "\\AVG_CCACYT.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x1.append(i * save_interval)
            y1.append(float(row[0]))
            i = i + 1
    with open(path2 + "\\AVG_CCACYT.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x2.append(j * save_interval)
            y2.append(float(row[0]))
            j = j + 1
    plt.figure()
    plt.grid()
    l1, = plt.plot(x1, y1)
    l2, = plt.plot(x2, y2)
    plt.gca().yaxis.set_major_formatter(ticker.FormatStrFormatter('%.8f'))
    plt.gca().xaxis.get_major_formatter().set_scientific(False)
    plt.title("Fn")
    # plt.xlim((0,i*save_interval))
    plt.legend(handles=[l1, l2], labels=["迭代法", "Fortran"])
    plt.savefig("Figure\\avg_fn.jpg", bbox_inches='tight')
    plt.show()

    x3 = []
    y3 = []
    x4 = []
    y4 = []
    i = 0
    j = 0
    with open(path1 + "\\AVG_CCAF.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x3.append(i * save_interval)
            y3.append(float(row[0]))
            i = i + 1
    with open(path2 + "\\AVG_CCAF.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x4.append(j * save_interval)
            y4.append(float(row[0]))
            j = j + 1
    plt.figure()
    plt.grid()
    l3, = plt.plot(x3, y3)
    l4, = plt.plot(x4, y4)
    # plt.xlim((0, j * save_interval))
    plt.title("Gn")
    plt.gca().xaxis.get_major_formatter().set_scientific(False)
    plt.legend(handles=[l3, l4], labels=["迭代法", "Fortran"])
    plt.savefig("Figure\\avg_gn.jpg", bbox_inches='tight')
    plt.show()


# *****************************************************************

def plot_center():
    """
    提取F0
    """
    fn = open(path + "\\F0.csv", "w")
    # gn = open(path + "\\AVG_G0.csv", "w")
    for filename in os.listdir(path + "\\SPARK"):
        print("读取" + filename)
        line = linecache.getline(path + "\\SPARK\\" + filename, 1)
        data = line.strip("\n").split(",")
        fn.write(str(data[0]) + "\n")
        # gn.write(str(data[1]) + "\n")
    fn.close()
    # gn.close()

    x = []
    y = []
    i = 1
    with open(path + "\\F0.csv", 'r') as ca_csvfile:
        plots = csv.reader(ca_csvfile)
        for row in plots:
            x.append(i * save_interval)
            y.append(float(row[0]))
            i = i + 1
    plt.figure()
    plt.grid()
    plt.plot(x, y)
    plt.savefig(path + "\\f0.jpg", bbox_inches='tight')
    plt.show()


# **********************************************************************
def plot_one_file():
    """
    根据某一个文件，绘制所有点值的曲线图
    """
    x1 = []
    y1 = []
    x2 = []
    y2 = []
    i = 0
    j = 0
    with open(path + "\\SPARK\\" + file, 'r') as ca_ccycto:
        data = csv.reader(ca_ccycto)
        for row in data:
            if i < 1000:
                x1.append(i)
                y1.append(float(row[0]))
                i = i + 1
    plt.figure()
    plt.grid()
    plt.plot(x1, y1)
    plt.xlim((0, i))
    plt.title(file + " CCACYT")
    plt.show()

    with open(path + "\\SPARK\\" + file, 'r') as ca_ccycto:
        data = csv.reader(ca_ccycto)
        for row in data:
            if j < 1000:
                x2.append(j)
                y2.append(float(row[1]))
                j = j + 1
    plt.figure()
    plt.grid()
    plt.plot(x2, y2)
    plt.xlim((0, j))
    plt.title(file + " CCAF")
    plt.show()


if __name__ == '__main__':
    # exact_average()
    plot_average()
    #

    # plot_one_file()

    # plot_center()

    # plot_average_multiple()
