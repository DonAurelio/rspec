# -*- encoding: utf-8 -*-

from signalprocessing import *
import time
import csv

import numpy
import matplotlib.pyplot as plt

if __name__ == '__main__':
    file = open('KD003.003','r')
    # file = open('GO003.001','r')
    lines = file.readlines()

    DATA = list(map(lambda x: float(x),lines))
    #DATA = volts_to_gales(DATA,1.25)
    DAMPING = 0.05
    SENSITIVITY = 1.25
    SAMPLE_RATE = 0.01

    cstart = time.time()
    # Response Espectrum Diego
    # signal: vector de carga externa
    # damping    
    # sv, device sensitivity
    # fm, sample rate

    # T: oscilation periods
    # Pd: Pseudo displacement response spectrum
    # Pv: Pseudo displacement response spectrum
    # Pa: Pseudo displacement response spectrum
    # T, Pd, Pv, Pa = diego(signal=DATA,damping=DAMPING,sv=SENSITIVITY,fm=SAMPLE_RATE)
    # plt.plot(T,Pd)

    # Response of a system with Duhamel Integral
    # p: vector de carga externa
    # m: masa del sistema
    # w: frecuencia natural del sistema
    # xi: fracción de amortiguamiento viscoso
    # dt: paso de timepo ¿o frecuencia de muestreo?

    # t: vector de timepo
    # d: desplazamiento de respuesta
    T,D = duhamel(P=DATA,m=1.0,w=2.0*numpy.pi,xi=DAMPING,dt=SAMPLE_RATE)
    plt.plot(T,D)
    
    # Response of a system with the Lineal Acceleration Method
    # p: vector de carga externa
    # m: masa del sistema
    # w: frecuencia natural del sistema
    # xi: fracción de amortiguamiento viscoso
    # dt: paso de timepo

    # t: vector de timepo
    # d: desplazamiento de respuesta
    # v: velocidad en respuesta
    # a: acceleración en respuesta
    # t, d, v, a = dmaclin(p=data,m=1.0,w=2.0*numpy.pi,xi=0.05,dt=0.01)

    # Response espetrum calculation with Linear Acceleration Method
    # acc: vector de la acceleración del suelo
    # tmin: periodo mínimo del calculo
    # tmax: periodo máximo del calculo
    # dt_period: incremento del periodo
    # vxi: fraccion de amortiguamiento
    # viscoso para las cuales se han de calcular los espectros
    # dt_accelerogram: paso del timepo del accelerograma

    # Sd: espectro de desplzamiento
    # Sv: espectro de velocidad
    # Sa: espectro de acceleración
    # t,d= desplin(acc=data,tmin=0.0,tmax=4.0,dt_period=0.02,xi=0.05,dt_accelerogram=0.01)

    # Response spectrum with transferFunction and lsim
    # A: Señal de acceleraciones
    # Z: amortiguamiento
    # Fs: Frecuencia de adqusición
    # Tn: Periodo

    # d: Pseudo acceleración
    # t, d = rspect(A=data,z=0.05,fs=0.01)

    # Response espectrum with simplified lineal acceleration method
    # p: es el accelerograma
    # m: es la masa del sistema de 1 gdl
    # c: es el amortiguamiento del sistema de 1 gdl
    # k: es la rigidez del sistema de 1 gdl
    # dt: es el incremento de tiempo con el cual se desea hallar
    # la respuesta. El mismo que tiene que ser igual al incremento
    # de tiempo con el cual se obtuvo el accelerograma.

    # d: desplazamiento
    # v: velocidad
    # a: acceleración del sistema
    # t, d, v, a = lineal(p=data,m=1.0,c=0.05,k=1.0,dt=0.01)


    cend = time.time()

    # wstart = time.time()
    # with open('result.csv','wb') as csvfile:
    #   spamwriter = csv.writer(csvfile,delimiter=';',quotechar='|', quoting=csv.QUOTE_MINIMAL)
    #   write = spamwriter.writerow
    #   for t,d in zip(ts,ds):
    #       write([t,d])
    # wend = time.time()

    plt.show()

    print "Calculation time", (cend-cstart)
    # print "Writting time" , (wend-wstart)
