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

def acc_time_serie(signal,fm):
    dt = 1.0/fm;
    signal_duration = dt * len(signal)
    return numpy.arange(start=0.00,stop=(signal_duration+dt),step=dt).tolist()

def spec_period_serie():
    return numpy.arange(start=0.02,stop=(4.0+0.02),step=0.02).tolist()

def rspectrum(signal,damping,sv,fm):
    """
        signal: Señal en voltios
        damping: Razón de amortiguamiento
        sv: Sensibilidad del sensor
        fm: Frecuencia de muestreo
    """

    # DATOS INICIALES
    # Periodo de cada oscilador que seran sometidos a la señal de acceleraciones
    ts = arange(start=0.02,stop=(4.0+0.02),step=0.02) 

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

        Wn=(2.0*pi)/t
        k=M*(Wn**2.0)
        d=DT
        a1=(2.71828)**(-z*Wn*d)
        Wd=Wn*sqrt(1.0-(z**2.0))
        a2=sin(Wd*d)
        a3=cos(Wd*d)
        
        A=a1*(((z/(sqrt(1.0-z**2.0)))*a2)+a3)
        B=a1*((1.0/Wd)*a2)
        C=(1.0/k)*((2.0*z/(Wn*d))+a1*((((1.0-2.0*(z**2))/(Wd*d))-(z/(sqrt(1.0-z**2.0))))*a2-(1.0+(2.0*z/(Wn*d)))*a3))
        D=(1.0/k)*(1.0-(2.0*z/(Wn*d))+a1*(((2.0*(z**2.0)-1.0)/(Wd*d))*a2+(2.0*z/(Wn*d))*a3))
        A1=-a1*(((Wn/(sqrt(1.0-z**2.0)))*a2))
        B1=a1*(a3-((z/(sqrt(1.0-z**2.0)))*a2))
        C1=(1.0/k)*((-1.0/d)+a1*(((Wn/(sqrt(1.0-(z**2.0))))+(z/(d*sqrt(1.0-(z**2.0)))))*a2+(a3/d)))
        D1=(1.0/(k*d))*(1.0-a1*((z/(sqrt(1.0-(z**2.0))))*a2+a3))

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

def duhamel(p,m,w,xi,dt):
    """ Duhamel integral calculation.

    To determine the general response of a 
    simple lineal system. 
    
    by: Jorge E. Hurtado G.
    Universidad Nacional de Colombia

    p: vector de carga externa
    m: masa del sistema
    w: frecuencia natural del sistema
    xi: fracción de amortiguamiento viscoso
    dt: paso de timepo

    t: vector de timepo
    d: desplazamiento de respuesta
        
    """

    n = len(p)
    tmax = dt*n
    t = numpy.linspace(start=0.0,stop=tmax,num=n)
    # t = numpy.arange(start=0.0,stop=(tmax),step=dt)

    wa = w*numpy.sqrt(1-xi**2.0)

    f = map((lambda x,y: x*y),p,numpy.cos(wa*t))
    g = map((lambda x,y: x*y),p,numpy.sin(wa*t))
    
    f1 = [0.0] + f[:-1]
    g1 = [0.0] + g[:-1]

    pc = map((lambda x,y: x*numpy.exp(-xi*w*dt)+y),f1,f)
    ps = map((lambda x,y: x*numpy.exp(-xi*w*dt)+y),g1,g)

    pc = map((lambda x: ((((x*dt)/m)/wa)/2.0) ),pc)
    ps = map((lambda x: ((((x*dt)/m)/wa)/2.0) ),ps)

    c = [pc[0]]
    s = [ps[0]]

    for i in range(1,n):
        c += [ c[i-1]*numpy.exp(-xi*w*dt)+pc[i] ] 
        s += [ c[i-1]*numpy.exp(-xi*w*dt)+pc[i] ]

    d = map((lambda x,y,z,w: x*y-z*w),c,numpy.sin(wa*t),s,numpy.cos(wa*t))

    return t.tolist(), d

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
        t, d,v,a = dmaclin(map( (lambda x: x*(-1.0)) ,acc), 1.0,w,xi,dt_accelerogram) 
        Sd.append(max(map(lambda x: abs(x),d)))
        Sv.append(max(map(lambda x: abs(x),v)))
        Sa.append(max(map(lambda x: abs(x),a)))

    return T.tolist(), Sd, Sv, Sa

# Arrays in mayus
def rspect(a,z,fs,tn):
    """
        by Omar

        A: Señal de acceleraciones
        Z: amortiguamiento
        Fs: Frecuencia de adqusición
        Tn: Periodo
    """

    t = [ (x*1/fs) for x in range(0,len(a))]

    PA = []

    for i in range(tn):
        T = tn[i]
        wn = (2.0*numpy.pi)/T

        H = sg.TransferFunction(num=[0.0, 0.0, -1.0], den=[1.0, 2.0*z*wn, wn**2])
        H2 = sg.TransferFunction(num=[-1.0, 0.0, 0.0], den=[1.0, 2.0*z*wn, wn**2])

        De = sg.lsim(system=H,U=a,T=T)
        Ac = sg.lsim(system=H2, U=a,T=T)

        D = max(map(lambda x: abs(x),De))
        EA = max(map(lambda x: abs(x),AC)) # No sé donde se usa esta variable

        PA.append( D * (2*numpy.pi/tn[i]) ** 2/9.81 )

        return t, PA


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
