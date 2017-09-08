#! .env/bin/python
# -*- encoding: utf-8 -*-

from signalprocessing import duhamel, dmaclin, volts_to_gales, desplin
import time
import csv

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file = open('GO003.001','r')
    lines = file.readlines()

    data = list(map(lambda x: float(x),lines))
    data = volts_to_gales(data,1.25)

    cstart = time.time()
    pi = 3.14159265359
    #t0,d0 = duhamel(p=data,m=1.0,w=2.0*pi,xi=0.05,dt=0.01)
    #t1,d1,v1,a1 = dmaclin(p=data,m=1.0,w=2.0*pi,xi=0.05,dt=0.01)
    t3,sd3,sv3,sa3 = desplin(acc=data,tmin=0.0,tmax=4.0,dt_period=0.02,xi=0.05,dt_accelerogram=0.01)
    cend = time.time()

    # wstart = time.time()
    # with open('result.csv','wb') as csvfile:
    #   spamwriter = csv.writer(csvfile,delimiter=';',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #   write = spamwriter.writerow
    #   for t,d in zip(ts,ds):
    #       write([t,d])
    # wend = time.time()

    plt.plot(t3, sd3)
    plt.show()

    print "Calculation time", (cend-cstart)
    # print "Writting time" , (wend-wstart)
