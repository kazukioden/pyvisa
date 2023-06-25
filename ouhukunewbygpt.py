import pyvisa
import numpy as np
import time
import matplotlib.pyplot as plt
import decimal
from datetime import datetime
import subprocess
gnuplot_processes = []











rm = pyvisa.ResourceManager()
inst = rm.open_resource("GPIB1::1::INSTR")
inst.write('*IDN?') 
print(inst.read()) 
vs=0.0
vl=600
dv=1
nall=(vl-vs)/dv
nalli=round(nall)
num=0
vn=0.0
nowvr=[]

inst.write("*CLS")
inst.write("FUNC 'CURR:DC'")
inst.write("SOUR:VOLT:MCON on")
inst.write("SENS1:CURR:DC:NPLC 1")
inst.write("SENS1:CURR:DC:RANG 20e-3")
inst.write("SOUR:VOLT:RANG 1000")
inst.write("SENS1:RES:MAN:VSO:OPER on")

i=0
while i<1:
    num=0
    time0=time.time()
    while  num < 2*nall+2:
        if abs(vn) >= 300 and abs(value) <= 1e-7:
            print("error air discharge")
            inst.write("SENS1:RES:MAN:VSO:OPER off")
            break
        if num <= nall:
            vn=vs+(-1)**i*num*dv
            inst.write("SOUR:VOLT {}".format(vn))
            nowv=datetime.now()
            inst.write("FORM:ELEM READ")
            inst.write("FETC?")
            value=str(inst.read("FETC?"))
            t=time.time()
            t=t-time0
            value=float(value.replace("E","e"))
            decimal.getcontext().prec=4
            vn=decimal.Decimal(vn)
            value=decimal.Decimal(value)
            print("Voltage "+ str(+vn) +" V"+"      Current "+str(1000*value)+ " mA")
            nowvr.append(str(nowv))
            if num==0:
              now=datetime.now()
              now=now.replace(microsecond=0)
              now=str(now)
              now=now.replace(":"," ")
              cap=(-1)**i*700
              xytfilename=now+'xytdata 0-'+str(cap)+ 'v '+str(i) +'th.txt'
              xytfilename=xytfilename.replace("'","")
              import os
              xytpath=r'C:/Users/Kagawa_lab/Desktop/python_lesson/shimizu_data/'
              xytpath=os.path.join(xytpath,xytfilename)
            with open(xytpath, 'a') as f:
                if f.tell() == 0:  # ファイルが空の場合（ファイルポインタが先頭の場合）
                 f.write("Voltage(V),Current(mA),UNIX_Time\n")  # ヘッダ行を書き込む
                f.write("{},{},{}\n".format(vn, value*1000, t))
            if num==0:
             datafile1=xytpath
             # gnuplotにプロットコマンドを送信
             # テンプレートを読み込む
             with open('template1.gnuplot') as f:
                template1 = f.read()

             # データファイル名を挿入して新しいスクリプトを作成
                script = template1.format(datafile1=os.path.abspath(datafile1))

             # 新しいスクリプトを一時的なファイルに保存
             with open('temp1.gnuplot', 'w') as f:
                f.write(script)

             # gnuplotでスクリプトを実行
             p1 = subprocess.Popen(['gnuplot', 'temp1.gnuplot'])
            
          
        else:
            nums=2*nall+1-num
            vn=vs+(-1)**i*nums*dv
            inst.write("SOUR:VOLT {}".format(vn))
            nowv=datetime.now()
            inst.write("FORM:ELEM READ")
            inst.write("FETC?")
            value=str(inst.read("FETC?"))
            t2=time.time()
            t2=t2-time0
            value=float(value.replace("E","e"))
            decimal.getcontext().prec=4
            vn=decimal.Decimal(vn)
            value=decimal.Decimal(value)
            print("Voltage "+ str(+vn) +" V"+"      Current "+str(1000*value)+ " mA")
            nowvr.append(str(nowv))
            if num==nall+1:
              now=datetime.now()
              now=now.replace(microsecond=0)
              now=str(now)
              now=now.replace(":"," ")
              cap=(-1)**i*700
              xytfilename=now+'xytdata'+str(cap)+'-0v'+ ' '+str(i) +'th.txt'
              xytfilename=xytfilename.replace("'","")
              xytpath=r'C:/Users/Kagawa_lab/Desktop/python_lesson/shimizu_data/'
              xytpath2=os.path.join(xytpath,xytfilename)
            with open(xytpath2, 'a') as f:
                if f.tell() == 0:  # ファイルが空の場合（ファイルポインタが先頭の場合）
                 f.write("Voltage(V),Current(mA),UNIX_Time\n")  # ヘッダ行を書き込む
                f.write("{},{},{}\n".format(vn, value*1000, t2))
            if num==nall+1:
             datafile2=xytpath2
             # gnuplotにプロットコマンドを送信
             # テンプレートを読み込む
             with open('template2.gnuplot') as f:
                template2 = f.read()

             #データファイル名を挿入して新しいスクリプトを作成
                script = template2.format(datafile1=os.path.abspath(datafile1),datafile2=os.path.abspath(datafile2))
             # 新しいスクリプトを一時的なファイルに保存
             with open('temp2.gnuplot', 'w') as f:
                f.write(script)

            # gnuplotでスクリプトを実行
             p2 = subprocess.Popen(['gnuplot', 'temp2.gnuplot'])
             
             #if num==2*nall+1:
                #p.kill()
        num=num+1
    i=i+1

inst.write("SENS1:RES:MAN:VSO:OPER off")
now=datetime.now()
now=now.replace(microsecond=0)
now=str(now)
now=now.replace(":"," ")
timefilename=now+'timefile -700kara700madeouhuku.txt'
timefilename=timefilename.replace("'","")
timefilepath=os.path.join(xytpath,timefilename)
timef=open(timefilepath,mode='w')
timef.write("Timestamp\n")  # ヘッダ行を書き込む
timef.write("\n".join(nowvr))
