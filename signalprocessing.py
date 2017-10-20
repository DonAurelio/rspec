# -*- encoding: utf-8 -*-

from scipy import signal as sg
from signalfilters import resampling
from signalfilters import bandpassfilter
import numpy

GALES = 981.0 # Gravedad (cm/s2 o gales)

def volts_to_gales(signal,sv):
    """
        signal: Señal en voltios 
        sv: Sensibilidad del dispositivo
    """

    # Paso de señal de voltios a cm/s2
    FE = GALES/sv # Conversión unidades fisicas
    signal = map(lambda x:float(x) * FE ,signal[:])

    # Corrección de linea base orden 0
    # signal = map(lambda x: x - mean(signal),signal)
    signal = sg.detrend(signal).tolist()

    return signal

def max_peak(ts,accls):
    return max(zip(ts,accls), key=lambda item: abs(item[1]))

def diego(signal,damping,sv,fm):
    """
        signal: Señal en voltios
        damping: Razón de amortiguamiento
        sv: Sensibilidad del sensor
        fm: Frecuencia de muestreo
    """

    # DATOS INICIALES
    # Periodo de cada oscilador que seran sometidos a la señal de acceleraciones
    ts = numpy.arange(start=0.02,stop=(4.0+0.02),step=0.02) 

    T = 1.0 # Periodo
    z = damping # Razón de amortiguamiento
    M = 1.0 # Masa
    DT = 1.0/fm # Delta tiempo

    ys = signal[:]

    # FILTRADO DE SEÑAL
    # ys = bandpassfilter(signal=ys,densidad=(fm/2.0)) # Pasabanda
    filtered_signal = ys[:]
    
    # Tabla datos espectro de respuesta en desplzamiento
    Dmaxs = []
    # Tabla datos espectro de respuesta en Pseudovelocidad
    Vmaxs = []
    # Tabla datos espectro de respuesta en Pseudoaceleración
    Amaxs = []
    # Periodos de oscilación (es el mismo para todas tablas)
    Ts = []

    # Coefcioentes de la ecuación diferencial
    uis = []
    # Derivadas de los coeficientes de la ecuación diferencial
    upis = []

    # Por cada periodo de oscilación
    for t in ts:

        Wn=(2.0*numpy.pi)/t
        k=M*(Wn**2.0)
        d=DT
        a1=(2.71828)**(-z*Wn*d)
        Wd=Wn*numpy.sqrt(1.0-(z**2.0))
        a2=numpy.sin(Wd*d)
        a3=numpy.cos(Wd*d)
        
        A=a1*(((z/(numpy.sqrt(1.0-z**2.0)))*a2)+a3)
        B=a1*((1.0/Wd)*a2)
        C=(1.0/k)*((2.0*z/(Wn*d))+a1*((((1.0-2.0*(z**2))/(Wd*d))-(z/(numpy.sqrt(1.0-z**2.0))))*a2-(1.0+(2.0*z/(Wn*d)))*a3))
        D=(1.0/k)*(1.0-(2.0*z/(Wn*d))+a1*(((2.0*(z**2.0)-1.0)/(Wd*d))*a2+(2.0*z/(Wn*d))*a3))
        A1=-a1*(((Wn/(numpy.sqrt(1.0-z**2.0)))*a2))
        B1=a1*(a3-((z/(numpy.sqrt(1.0-z**2.0)))*a2))
        C1=(1.0/k)*((-1.0/d)+a1*(((Wn/(numpy.sqrt(1.0-(z**2.0))))+(z/(d*numpy.sqrt(1.0-(z**2.0)))))*a2+(a3/d)))
        D1=(1.0/(k*d))*(1.0-a1*((z/(numpy.sqrt(1.0-(z**2.0))))*a2+a3))

        uis = [0.0]
        upis = [0.0] 

        # Por cada una de las acceleraciones (Columna 2 del archivo)
        for i in range(len(ys)-1): # Revisar por que ? el - 1
            uis.append( A*uis[i]+ B*upis[i] + C*ys[i] + D*ys[i+1] )
            upis.append( A1*uis[i] + B1*upis[i] + C1*ys[i] + D1*ys[i+1]  )
            
        # Tabla Dmaxs 
        dmax = max(map(lambda x: abs(x),uis))
        Dmaxs.append(dmax)
        # Tabla Vmaxs 
        vmax = max(map(lambda x: abs(x*Wn),uis))
        Vmaxs.append(vmax)
        # Tabla Amaxs
        amax = max(map(lambda x: abs( (x*(Wn**2)) ),uis))
        Amaxs.append(amax)
        # Tabla con valores de periodos de oscilación
        Ts.append(t)

    return Ts,Dmaxs,Vmaxs,Amaxs

def duhamel(P,m,w,xi,dt):
    """ Duhamel integral calculation.

    Version: 3.0
    Calculation time 0.349811077118
    To determine the general response of a 
    simple lineal system. 
    
    by: Jorge E. Hurtado G.
    Universidad Nacional de Colombia

    P: vector de carga externa
    m: masa del sistema
    w: frecuencia natural del sistema
    xi: fracción de amortiguamiento viscoso
    dt: paso de timepo

    t: vector de timepo
    d: desplazamiento de respuesta
        
    """

    n = len(P)
    tmax = dt*n
    T = numpy.linspace(start=0.0,stop=tmax,num=n)
    # T = numpy.arange(start=0.0,stop=(tmax),step=dt)

    wa = w*numpy.sqrt(1-xi**2.0)

    # Alternatives scalar and numpy array multiplication
    # wa_x_T = map((lambda x: x*wa),T)
    # wa_x_T = numpy.multiply(T,wa)
    wa_x_T = T*wa

    # Aternative python list by numpy array multiplication (all has the same behavior)
    # F = map((lambda x,y: x*y),P,numpy.cos(wa_x_T))
    # F = numpy.multiply(P,numpy.cos(wa_x_T))
    F = P * numpy.cos(wa_x_T)

    # Aternative python list by numpy array multiplication (all has the same behavior)
    # G = map((lambda x,y: x*y),P,numpy.sin(wa_x_T))
    G = numpy.multiply(P,numpy.sin(wa_x_T))
    # G = P * numpy.sin(wa_x_T)

    #  the first value differs 
    # for f in G:
    #     print round(f,4)
    
    # Insert, (array), (index), (value)
    F1 = numpy.insert(F,0,0.0)[:-1]
    G1 = numpy.insert(G,0,0.0)[:-1]

    # the second value differs
    # for value in G1:
    #     print round(value,4)

    # Pc = map((lambda x,y: x*numpy.exp(-xi*w*dt)+y),F1,F)
    Pc = F1 * numpy.exp(-xi*w*dt) + F

    # Ps = map((lambda x,y: x*numpy.exp(-xi*w*dt)+y),G1,G)
    Ps = G1 * numpy.exp(-xi*w*dt) + G

    # Pc = map((lambda x: ((((x*dt)/m)/wa)/2.0) ),Pc)
    Pc = Pc  * dt 
    # / m / wa / 2.0

    # ERROR NO COINCIDE
    for value in Pc:
        print round(value,5)
    
    # Ps = map((lambda x: ((((x*dt)/m)/wa)/2.0) ),Ps)
    Ps = Ps * dt / m / wa / 2.0

    C = [Pc[0]]
    S = [Ps[0]]

    for i in range(1,n):
        C += [ C[i-1]*numpy.exp(-xi*w*dt)+Pc[i] ] 
        S += [ S[i-1]*numpy.exp(-xi*w*dt)+Ps[i] ]

    # Both answers are different ¿?
    # d = map((lambda x,y,z,w: x*y-z*w), C, numpy.sin(wa_x_T), S, numpy.cos(wa_x_T))
    d = C * numpy.sin(wa_x_T) + S * numpy.cos(wa_x_T) # I beleive this is the correct

    return T.tolist(), d

def dmaclin(p,m,w,xi,dt):
    """ Lineal Acceleration Method.

    To determine the general response of a 
    simple lineal system by the lineal 
    acceleration method, it is an alternative 
    to the Duhamel method.
    
    by: Jorge E. Hurtado G.
    Universidad Nacional de Colombia

    p: vector de carga externa
    m: masa del sistema
    w: frecuencia natural del sistema
    xi: fracción de amortiguamiento viscoso
    dt: paso de timepo

    t: vector de timepo
    d: desplazamiento de respuesta
    v: velocidad en respuesta
    a: acceleración en respuesta
    """

    n = len(p)
    tmax = dt*n
    t = numpy.linspace(start=0.0,stop=tmax,num=n)
    d0 = 0.0
    v0 = 0.0
    a0 = 0.0

    k = m*w**2.0
    c = 2.0*m*w+xi
    kbar = (k+((3.0*c)/dt)) + ((6.0*m)/(dt**2.0))
    ikbar = 1.0/kbar

    d = []
    v = []
    a = []

    for i in range(n):
        p1 = p[i]
        dp = m*( ((6.0*d0)/(dt**2.0))+((6.0*v0)/dt)+(2.0*a0) )
        dp = dp + ( c*((3.0*d0)/dt)+(2.0*v0)+((dt*a0)/2.0) )
        pbar = p1 + dp
        d1 = ikbar * pbar
        v1 = ((3.0*(d1-d0))/dt) - (2.0*v0) - ((dt*a0)/2.0)
        a1 = ((6.0*(d1-d0))/dt**2.0) - ((6.0*v0)/dt) - (2.0*a0)
        d.append(d1)
        v.append(v1)
        a.append(a1)
        d0 = d1
        v0 = v1
        a0 = a1

    return t.tolist(), d, v, a

def desplin(acc,tmin,tmax,dt_period,xi,dt_accelerogram):
    """Calcula los espectros de respuesta de un sistema 
    sencillo lineal por el método de la acceleración lineal.

    Por: Jorge E. Hurtado G.
        Universidad Nacional de Colombia

    acc: vector de la acceleración del suelo
    tmin: periodo mínimo del calculo
    tmax: periodo máximo del calculo
    dt_period: incremento del periodo
    vxi: fraccion de amortiguamiento
    viscoso para las cuales se han de calcular los espectros
    dt_accelerogram: paso del timepo del accelerograma

    Sd: espectro de desplzamiento
    Sv: espectro de velocidad
    Sa: espectro de acceleración
    """
    m = (tmax-tmin) / dt_period+1.0
    print m
    T = numpy.linspace(start=tmin,stop=tmax,num=m)
    # error por division por cero, el primer valor de la lista es inf
    for i in T:
        print (2.0*numpy.pi)/i
    W = map((lambda x: (2.0*numpy.pi)/x ), T)

    Sd = []
    Sv = []
    Sa = []

    for j in range(int(m)):
        w = W[j]
        # t, d,v,a = dmaclin(map( (lambda x: x*(-1.0)) ,acc), 1.0,w,xi,dt_accelerogram)
        # Sd.append(max(map(lambda x: abs(x),d)))
        # Sv.append(max(map(lambda x: abs(x),v)))
        # Sa.append(max(map(lambda x: abs(x),a)))

     # return T.tolist(), Sd, Sv, Sa

        t, d = duhamel(map( (lambda x: x*(-1.0)) ,acc), 1.0,w,xi,dt_accelerogram) 
        Sd.append(max(map(lambda x: abs(x),d)))

    return T.tolist(), Sd

def rspect(A,z,fs):
    """
        by Omar

        A: Señal de acceleraciones
        z: amortiguamiento
        fs: Frecuencia de adqusición o muestreo
    """

    # ¿Qué es T?
    T = [ (x*1/fs) for x in range(0,len(A))]
    # Pseudo acceleration response spectrum
    Pa = []
    #  Oscilator periods for start=0.0,stop=4.0,step=0.02
    Ts = spec_period_serie()
 
    for i in range(len(Ts)):
        t = Ts[i]
        wn = (2.0*numpy.pi)/t

        # H = sg.TransferFunction(num=[0.0, 0.0, -1.0], den=[1.0, 2.0*z*wn, wn**2])
        # H2 = sg.TransferFunction(num=[-1.0, 0.0, 0.0], den=[1.0, 2.0*z*wn, wn**2])
        H = sg.TransferFunction([-1.0, 0.0, 0.0],[wn, 2.0*z*wn, 1.0])
        H2 = sg.TransferFunction([0.0, 0.0, -1.0],[wn, 2.0*z*wn, 1.0 ])

        # tout, a 1D array, time values for the output
        # yout, 1D array system response
        # xout, time evolution of the state vector
        # displacement response
        Tout, Yout, Xout = sg.lsim(system=H,U=A,T=T)
        # acceleration response
        Tout1, Yout1, Xout1 = sg.lsim(system=H2, U=A,T=T)
        D = max(map(lambda x: abs(x),Yout)) 
        EA = max(map(lambda x: abs(x),Yout1)) # No sé donde se usa esta variable
        
        Pa.append( D * (2*numpy.pi/Ts[i]) ** 2/9.81 ) # ¿Es acceleración o desplazamiento?

    return Ts, Pa

def lineal(p,m,c,k,dt):
    """Respuesta en tiempo de un sistema de un grado de
    de libertad por método de la acceleración lineal

    by Dinámica de Estructuras con MATLAB

    p: es el accelerograma
    m: es la masa del sistema de 1 gdl
    c: es el amortiguamiento del sistema de 1 gdl
    k: es la rigidez del sistema de 1 gdl
    dt es el incremento de tiempo con el cual se desea hallar
    la respuesta. El mismo que tiene que ser igual al incremento
    de tiempo con el cual se obtuvo el accelerograma.

    d: desplazamiento
    v: velocidad
    a: acceleración del sistema
    """

    n = len(p)
    tmax = dt*n
    t = numpy.linspace(start=0,stop=tmax,num=n)
    ma = m + (c*dt/2) + (k*dt*(dt/6))

    d = [0]
    v = [0]
    a = [0]

    # range exclide the last value (n)
    for i in range(0,n-1):
        dq = -m*(p[i+1]-p[i])
        dqa = dq - a[i]*(c*dt+k*dt*(dt/2)) - v[i]*k*dt
        inca = dqa/ma
        incv = a[i]*dt+inca*dt/2
        incd = v[i]*dt+a[i]*dt*dt/2+inca*dt*dt/6
        d.append(d[i]+incd)
        v.append(v[i]+incv)
        a.append(a[i]+inca)
        d[i] = d[i+1]
        v[i] = v[i+1]
        a[i] = a[i+1]

    return t.tolist(), d, v, a
