# -*- encoding: utf-8 -*-

import scipy
import numpy

def resampling(signal=[]):
	"""
		Remuestrea la señal entrante.

		Dada una señal con un paso de tiempo definido dt (para este caso
		0.005 ) permite interpolar la señal dado un paso de tiempo n_dt distinto.

		Cómo para el observatorio todas las señales  recibidas estan estandarizadas 
		a una dt de 0.005 entonces no es necesario hacer remuestreo.

			signal = señal o acceleraciones del sismo
	"""
	signal_array = signal
	dt = 0.005 # Paso del tiempo estandar en todas las señales de acceleraciones 0.005
	dur = len(signal_array) * dt # Duración del sismo (Número de datos por paso del tiempo)
	tiempo = numpy.linspace(0.0, dur, num=len(signal_array))

	n_dt = 0.5
	n_tiempo = numpy.linspace(0.0, dur, num=dur/n_dt)
	r_signal = numpy.interp(x=n_tiempo, xp=tiempo, fp=signal_array)

	return  r_signal


def bandpassfilter(signal,densidad):
	"""
		Butterworth filter

		Un filtro pasa banda es un tipo de filtro que 
		deja pasar un determinado rango de freceuncias 
		de una señal y atenual el paso del resto
			
			signal = Señal		
			fcinf, Frecuencia de corte inferior 
			fcsup, Frecuencia de corte superior
			densidad, mitad de cantidad de datos por segundo
	"""

	fcinf = 0.1 # Frecuencia de corte inferior 
	fcsup = 50.0 # Frecuencia de corte superior
	
	fbp = [fcinf/densidad,fcsup/densidad] # Vector de constantes para modificar señal
	
	# Reference http://docs.scipy.org/doc/scipy-0.14.0/reference/generated/scipy.signal.butter.html
	b,a = scipy.signal.butter(N=2,Wn=fbp,btype='bandpass') # Filtro pasabanda
	yT = scipy.signal.filtfilt(b=b,a=a,x=signal) # Señal ya filtrada entre 0.1 y 50 Hz

	return yT