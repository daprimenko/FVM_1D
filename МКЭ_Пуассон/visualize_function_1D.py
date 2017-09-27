import matplotlib.pyplot as plt
import numpy as np
import struct


def init():
    # Инициализирует необходимые переменные
    global fmt, length, NPoints
    fmt = 'd'
    length = struct.calcsize(fmt)
    length_int = struct.calcsize('i')

    # Открываем информационный файл
    with open('partition_info.dat', 'rb') as file:
        inf_dat = file.read(length_int)
        NPoints = struct.unpack('i', inf_dat)[0]


def plot_function():
    # Строит узловые точки
    x = np.zeros(NPoints)
    with open('partition.dat', 'rb') as file:
        for i in range(NPoints):
            x_dat = file.read(length)
            x[i] = struct.unpack(fmt, x_dat)[0]

    f = np.zeros(NPoints)
    with open('u.dat', 'rb') as file:
        for i in range(NPoints):
            f_dat = file.read(length)
            f[i] = struct.unpack(fmt, f_dat)[0]

    plt.plot(x, f, 'o', markersize=2, color='black')
    plt.show()


init()
plot_function()