import pyvisa
import numpy as np
import time
import matplotlib.pyplot as plt
import decimal
import LS335, KL6517B
import shimizu_file as smz
import keyboard  # 必要な場合はpipでインストール

lake = LS335.LS335(12)
elec = KL6517B.KL6517B(1)

elec.set_dc_range(dcrange=2e-3)
elec.operate(state=1)
vn = 0
elec.set_voltage(voltage=vn)
curr = elec.get_current()

# 条件設定
target_temp = 300
vmax = 700
filename = smz.genfilename(vmax=vmax)
smz.writeheader(filename)

# vnの掃引用変数
increment = 2  # 例として10Vずつ変動
direction = 1  # 1 = 増加, -1 = 減少
n_swings = 0  # vnの掃引の回数をカウント
max_swings = 1  # 往復の最大回数

# 温度設定用の変数
target_temp1 = 415
target_temp2 = 310

# 測定開始
t0 = time.time()
temp = lake.getTemp()
while True:
    try:
        systime = time.time()
        reltime = time.time() - t0

        vn += increment * direction
        if vn >= vmax:
            direction = -1
            n_swings += 1
        elif vn <= 0:
            direction = 1
            n_swings += 1

        elec.set_voltage(voltage=vn)
        curr = elec.get_current()
        temp = lake.getTemp()
        data = [temp, curr, systime, reltime, vn]
        smz.writedata(filename, data)

        if temp > 450 or n_swings >= 2 * max_swings:
            elec.operate(state=0)
            break

    except KeyboardInterrupt:
        print("Ctrl+C detected")
        #elec.operate(state=0)
        break
while True:
    try:
        lake.setTemp(1, target_temp1, 30)

        systime = time.time()
        reltime = time.time() - t0

        vn += increment * direction
        if vn >= vmax:
            direction = -1
            n_swings += 1
        elif vn <= 0:
            direction = 1
            n_swings += 1

        elec.set_voltage(voltage=vn)
        curr = elec.get_current()
        temp = lake.getTemp()
        data = [temp, curr, systime, reltime, vn]
        smz.writedata(filename, data)

        if temp > 450 or n_swings >= 2 * max_swings:
            elec.operate(state=0)
            break

    except KeyboardInterrupt:
        print("Ctrl+C detected")
        #elec.operate(state=0)
        break
 
while True:
    try:
        lake.setTemp(1, target_temp2, 30)

        systime = time.time()
        reltime = time.time() - t0

        vn += increment * direction
        if vn >= vmax:
            direction = -1
            n_swings += 1
        elif vn <= 0:
            direction = 1
            n_swings += 1

        elec.set_voltage(voltage=vn)
        curr = elec.get_current()
        temp = lake.getTemp()
        data = [temp, curr, systime, reltime, vn]
        smz.writedata(filename, data)

        if temp > 450 or n_swings >= 2 * max_swings:
            elec.operate(state=0)
            break

    except KeyboardInterrupt:
        print("Ctrl+C detected")
        elec.operate(state=0)
        break
