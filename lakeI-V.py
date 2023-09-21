import sys
sys.path.append("C:/Users/Keisuke_Matsuura/Desktop/Programming/python/shimizu_TCNQ/common")
import pyvisa
import numpy as np
import time
import matplotlib.pyplot as plt
import decimal
import LS335,KL6517B
import shimizu_file as smz

lake = LS335.LS335(12)
elec = KL6517B.KL6517B(1)

elec.operate(state=1)
vn = 0
elec.set_voltage(voltage=vn)
curr = elec.get_current()

# 条件設定
target_temp = 315
vmax = 500
filename=smz.genfilename(setdir="../data/",vmax=vmax)
smz.writeheader(filename)

# vnの掃引用変数
increment = 2  # 例として10Vずつ変動
direction = 1  # 1 = 増加, -1 = 減少
n_swings = 0  # vnの掃引の回数をカウント
max_swings = 1  # 往復の最大回数

#温度安定
#lake.wait_temp(target_temp=target_temp,interval=1)

#測定開始
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
        data = [temp,curr,systime,reltime,vn]
        smz.writedata(filename,data)

        if(temp > 450 or  n_swings >= 2*max_swings):
            elec.operate(state=0)
            break

    except KeyboardInterrupt:
        print("Ctrl+C detected")
        elec.operate(state=0)
        break
