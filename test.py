#!/usr/bin/python

from signalprocessing import duhamel, dmaclin, volts_to_gales
import time
import csv

import numpy as np
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file = open('KD003.002','r')
    lines = file.readlines()

    data = list(map(lambda x: float(x),lines))
    data = volts_to_gales(data,1.25)

    cstart = time.time()
    pi = 3.14159265359
    t0s,d0s = duhamel(p=data,m=1.0,w=2.0*pi,xi=0.05,dt=0.01)
    t1s,d1s,v1s,a1s = dmaclin(p=data,m=1.0,w=2.0*pi,xi=0.05,dt=0.01)
    cend = time.time()

    # wstart = time.time()
    # with open('result.csv','wb') as csvfile:
    #   spamwriter = csv.writer(csvfile,delimiter=';',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #   write = spamwriter.writerow
    #   for t,d in zip(ts,ds):
    #       write([t,d])
    # wend = time.time()

    plt.plot(t1s, d1s)
    plt.show()

    # print "Calculation time", (cend-cstart)
    # print "Writting time" , (wend-wstart)
