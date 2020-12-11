import scipy.signal as signal
from scipy.io.wavfile import write, read
import numpy as np
import xlrd
from xlutils.copy import copy
import scipy.signal as signal
import copy
import math

def invoke(sig, mov, sen, vel):
	"""
	Importa archivos .wav de prueba de acuerdo a sus caracteristicas
	"""
	
	fs, x = read('../Sonidos/viet/'+sig+mov+vel+sen+'m.wav')
	fs, y = read('../Sonidos/viet/'+sig+mov+vel+sen+'f.wav')
	fs = float(fs)
	return fs, x, y

def trim(X1, X2, fs, t):
	
	if len(X1)>t*fs and len(X2)>t*fs:
		X1 = X1[:int(t*fs)]
		X2 = X2[:int(t*fs)]
	
	elif len(X1)>len(X2):
		X1 = X1[:len(X2)]
	
	else:
		X2 = X2[:len(X1)]

	return X1, X2

def norm(x):
	"""
	Normaliza los archivos de audio poniendo el pico en el valor 1
	"""
	m = float(max(abs(x)))
	X = x/m
	return X
	
def filtroPB(n,f1,f2,fs,x,lv):
	"""
	Filtro pasa banda para archivos .wav mono.
	Input:
		n: Numero de muestras de la respuesta al impulso
		f1: Frecuencia de corte inferior
		f2: Frecuencia de corte superior
		fs: Frecuencia de muestreo
		x: Archivo .wav a filrar
	Output:
		array del archivo filtrado
	"""
	fs = float(fs)
	# f1 = fs/float((lv/2.0))  # Define la frecuencia de corte inferior de forma que entren al menos 2 ciclos por ventana
	# if f1 < 100:
	# 	f1 = 100
	a = signal.firwin(n, cutoff = f1/(fs/2.0), window = 'blackman')
	b = - signal.firwin(n, cutoff = f2/(fs/2.0), window = 'blackman')
	b[n/2] = b[n/2] + 1
	d = - (a+b)
	d[n/2] = d[n/2] + 1
		
	X = signal.fftconvolve(x, d, 'same')
	return X

def vinicial(X1, X2, lv, y1, y2):

	X11 = copy.deepcopy(X1)
	X22 = copy.deepcopy(X2)
	X11 = np.append(X11, np.zeros(lv-len(X11)%lv)) 
	X22 = np.append(X22, np.zeros(lv-len(X22)%lv)) 
	l   = len(X11)/lv

	rms1 = float(np.sqrt(np.mean(X11**2)))
	rms2 = float(np.sqrt(np.mean(X22**2)))
	X11  = X11/rms1
	X22  = X22/rms2

	X11  = X11[:l*lv]
	X111 = copy.deepcopy(X11)
	X22  = X22[:l*lv]
	nl   = np.floor(np.log2(len(X11)/lv)) #Numero de loops

	cmax  = np.zeros(2)
	cpos  = np.zeros(2)
	Cposa = 0
	X2pos = 0 # Posicion inicial del segmento final de X2
	j     = 0

	while j<nl:
		l2 = np.round(len(X111)/2)

		for i in range(2):
			u1      = l2*i
			u2      = u1 + l2
			X2C     = X22[u1:u2] # Corta la senal modificada en partes igualles
			Y       = signal.correlate(X111,X2C,'full','fft') # Correlacion entre senal modificada cortada y de referencia
			cmax[i] = float(max(Y)/float(len(Y))) # Vector de valores maximos de correlacion
			s       = len(Y)/2
			cpos[i] = np.argmax(Y) - s # Vector de posicion de valores maximos de correlacion
		
		Cposr = np.argmax(cmax) # Maxima posicion relativa
		X2pos = Cposr*l2 + X2pos
		Cposa = int(cpos[Cposr]) + Cposa   # Maxima posicion absoluta
		X22   = X22[l2*Cposr:l2*(Cposr+1)]
		# print Cposa

		mp = len(X11)/2+Cposa-l2/2  # Posicion inicial del vector X111 con respecto a X11
		
		if mp>=0:
			X111 = X11[mp:mp+l2]
		else:
			X111 = np.append(np.zeros(abs(mp)),X11[:l2+mp])
		
		# c = signal.correlate(X111,X22,'full','fft')	
		# print(np.argmax(c)-len(c)/2)
		
		j+=1

	n = mp-X2pos

	if n<0:
		z = np.zeros(abs(n))
		X1 = np.append(z,X1)
		y1 = np.append(z,y1)
		X2 = np.append(X2,z)
		y2 = np.append(y2,z)
		r = X2pos
	elif n>0:
		z = np.zeros(abs(n))
		X2 = np.append(z,X2)
		y2 = np.append(z,y2)
		X1 = np.append(X1,z)
		y1 = np.append(y1,z)
		r = mp

	return X1, X2, r, y1, y2

def cor2(X1, X2, lv, p, g, mv, fs, CCM, u):
	"""
	Correlacion cruzada de audios.
	lv es la longitud de la ventana temporal a utilizar. 
	p es el porcentaje de lv donde transiciona entre 0 y 1.
	g es el umbral del gate
	mv es la maxima velocidad de desplazamiento entre microfonos
	fs es la frecuencia de muestreo
	CCM es el coeficiente de correlacion minima
	"""

	X11 = copy.deepcopy(X1)
	X12 = copy.deepcopy(X2)	
	X11 = np.append(X11,np.zeros(lv-len(X11)%lv)) 
	X12 = np.append(X12,np.zeros(lv-len(X12)%lv)) 

	
	################## Ventana temporal #########################
	w  = np.blackman(round(lv*p)) 	 # Genera la ventana blackman
	w1 = np.append(w[:len(w)/2],np.ones(int(round(lv*(1-p))))) 
	w1 = np.append(w1, w[len(w)/2:]) # Crea la ventana temporal de longitud lv

	################## Generacion de arrays del bucle for #######
	l  = int(len(X11)/float(lv))  # Cantidad de ventanas temporales/ cilos for
	ta = np.zeros(l) 			  # Vector de corrimentos absolutos
	CC = np.zeros(l) 			  # Vector de maximas correlaciones por ventana
	M  = np.ones(l)   			  # Vector de coeficientes de correccion por umbral de amplitud
	N  = np.ones(l)   			  # Vector de coeficientes de correccion por umbral de correlacion

	################## Maximo Corrimiento Posible ###############
	h    = 2.0 							# Multiplicador del maximo corrimiento
	mcs  = math.ceil(fs*mv/343.0) 		# Maximo corrimiento en muestras por segundo
	mcv  = int(math.ceil(lv*mcs/fs)*h)  # Maximo corrimiento en muestras por ventana
	mcvl = int(math.ceil(mcv/2.0)) 		# Maximo corrimiento lateral

	################## Parametros del bucle #####################
	m    = 3									# Coeficiente de correccion por umbral de amplitud 
	mm   = 5									# Maximo valor de m
	n    = 10 									# Coeficiente de correccion por umbral de correlacion
	nm   = 5									# Maximo valor de n
	tt   = 0 							    	# Parametro de ajuste para la correlacion
	rms1 = float(np.sqrt(np.mean(X11**2))) 	# Valor rms de la senal X11
	rms2 = float(np.sqrt(np.mean(X12**2))) 	# Valor rms de la senal X12

	################## Bucle de correlacion #####################

	z  = np.arange(u/lv-1,-1,-1)					
	t0 = np.zeros(len(z)) 			            # Vector de corrimientos relativos a i-1
	for i in z:									# Calculo de corrimientos a la izquierda de la ventana inicial
		lv1  = int(lv*i) 						# Limite inferior
		lv2  = int(lv*(i+1))					# Limite superior 
		Y1   = X11[lv1:lv2]*w1 					# Senal de referencia ventaneada en tiempo
		RMS1 = float(np.sqrt(np.mean(Y1**2))) 	# Valor RMS de la ventana de referencia
		Y1   = Y1/RMS1 							# Senal de referencia normalizada a su valor RMS
		
		if RMS1/rms1>g:							# Gate de la senal de referencia
			r = int(np.ceil(mcvl*m*n)) 			# Cantidad de muestras del corrimiento lateral por ventana
			c = np.zeros(int(np.ceil(mcv*m*n)+1))	# Vector de correlaciones por ciclo del bucle
			
			for j in range(-r, r+1): # Ventaneo de la senal desplazada
			
				if i == 0 and j < 0: 		# Condicion inicial de corrimiento hacia la izquierda																					 
					Y2   = np.append(np.zeros(abs(j)),X12[:int(lv2+j)])*w1 # Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2)))				   # Valor RMS de Y2
					Y2   = Y2/RMS2 										   # Y2 Normalizada a su valor RMS

				elif i == 0 and j >= 0:		# Condicion inicial de corrimiento hacia la derecha
					Y2   = X12[j:int(lv2+j)]*w1  						   # Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2)))				   # Valor RMS de Y2
					Y2   = Y2/RMS2										   # Y2 Normalizada a su valor RMS

				elif len(X12[int(lv1-tt-j):int(lv2-tt-j)]) == lv: # Condicion normal
					Y2   = X12[int(lv1-tt-j):int(lv2-tt-j)]*w1 	# Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2))) 		# Valor cuadratico medio de Y2
					Y2   = Y2/RMS2								# Y2 Normalizada a su valor RMS

				else:
					T    = len(X12[int(lv1-tt-j):int(lv2-tt-j)])	 	                # Ultima ventana
					TT   = np.append(X12[int(lv1-tt-j):int(lv2-tt-j)],np.zeros(lv-T))	# Ventana Y2
					Y2   = TT*w1 														
					RMS2 = float(np.sqrt(np.mean(Y2**2)))   					    	# Valor cuadratico medio de Y2
					Y2   = Y2/float(RMS2)						 					    # Y2 Normalizada a su valor RMS

				if RMS2/rms2>g: 								    	# Gate de la senal desplazada
					m      = 1 											# Reset del coeficiente de corrimiento maximo por amplitud
					c[j+r] = signal.correlate(Y1, Y2,'valid','fft')     # Correlacion de senales
					C      = max(c)/float(lv) 							# Valor de correlacion maxima

					if C >= CCM: 			# Condicion para valores de correlacion por debajo del umbral CCM
						cm    = np.argmax(c)-r 	# Posicion del valor de maxima correlacion
						t0[i] = cm
						n     = 1				# Reset del coeficiente de corrimiento maximo por correlacion
						CC[i] = C				# Vector de valores maximos de correlacion

					else:
						cm     = 0				 
						CC [i] = CCM
					
						if n<nm:			
							n += 1
						
					N[i] = n # Vector de coeficientes de correccion por umbral de correlacion
				
				else:
					t0[i] = 0
					
					if m<mm:
						m += 1
		else:
			t0[i] = 0

			if m<mm:
				m += 1

		M[i] = m # Vector de coeficientes de correccion por umbral de amplitud
		
		####### Vectores de corrimientos tt y ta ########
		tt = sum(t0)
		ta[i] = tt

	t1 = np.zeros(2+l-u/lv) 			    # Vector de corrimientos relativos a i-1
	m  = 1									# Coeficiente de correccion por umbral de amplitud 
	n  = 10 									# Coeficiente de correccion por umbral de correlacion
	tt = 0 									# Parametro de ajuste para la correlacion
	z  = np.arange(u/lv,l)

	for i in z:
		lv1  = int(lv*i) 						# Limite inferior
		lv2  = int(lv*(i+1))						# Limite superior 
		Y1   = X11[lv1:lv2]*w1 					# Senal de referencia ventaneada en tiempo
		RMS1 = float(np.sqrt(np.mean(Y1**2))) 	# Valor RMS de la ventana de referencia
		Y1   = Y1/RMS1 							# Senal de referencia normalizada a su valor RMS
		
		if RMS1/rms1>g:							# Gate de la senal de referencia
			r = int(np.ceil(mcvl*m*n)) 			# Cantidad de muestras del corrimiento lateral por ventana
			c = np.zeros(int(np.ceil(mcv*m*n)+1))	# Vector de correlaciones por ciclo del bucle
	
			for j in range(-r, r+1): 			# Ventaneo de la senal desplazada
				if i == 0 and j < 0: 		# Condicion inicial de corrimiento hacia la izquierda																					 
					Y2   = np.append(np.zeros(abs(j)),X12[:int(lv2+j)])*w1 # Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2)))				 # Valor RMS de Y2
					Y2   = Y2/RMS2 										 # Y2 Normalizada a su valor RMS

				elif i == 0 and j >= 0:		# Condicion inicial de corrimiento hacia la derecha
					Y2   = X12[j:int(lv2+j)]*w1  							 # Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2)))				 # Valor RMS de Y2
					Y2   = Y2/RMS2										 # Y2 Normalizada a su valor RMS

				elif len(X12[int(lv1-tt-j):int(lv2-tt-j)]) == lv: # Condicion normal
					Y2   = X12[int(lv1-tt-j):int(lv2-tt-j)]*w1 	  # Ventana Y2
					RMS2 = float(np.sqrt(np.mean(Y2**2))) 		  # Valor cuadratico medio de Y2
					Y2   = Y2/RMS2								  # Y2 Normalizada a su valor RMS

				else:
					T    = len(X12[int(lv1-tt-j):int(lv2-tt-j)])	  					# Ultima ventana
					TT   = np.append(X12[int(lv1-tt-j):int(lv2-tt-j)],np.zeros(lv-T))	# Ventana Y2
					Y2   = TT*w1 														
					RMS2 = float(np.sqrt(np.mean(Y2**2)))   						# Valor cuadratico medio de Y2
					Y2   = Y2/float(RMS2)						 						# Y2 Normalizada a su valor RMS

				if RMS2/rms2>g: 									# Gate de la senal desplazada
					m      = 1 											# Reset del coeficiente de corrimiento maximo por amplitud
					c[j+r] = signal.correlate(Y1, Y2,'valid','fft') # Correlacion de senales
					C      = max(c)/float(lv) 							# Valor de correlacion maxima
					
					if C >= CCM: 				# Condicion para valores de correlacion por debajo del umbral CCM
						cm         = np.argmax(c)-r 	# Posicion del valor de maxima correlacion
						t1[i-u/lv] = cm
						n          = 1					# Reset del coeficiente de corrimiento maximo por correlacion
						CC[i]      = C				# Vector de valores maximos de correlacion

					else:
						cm     = 0				 
						CC [i] = CCM
					
						if n<nm:			
							n += 1
						
					N[i] = n # Vector de coeficientes de correccion por umbral de correlacion
				
				else:
					t1[i-u/lv] = 0
					
					if m<mm:
						m += 1

		else:
			t1[i-u/lv] = 0

			if m<mm:
				m += 1

		M[i] = m # Vector de coeficientes de correccion por umbral de amplitud
		
		####### Vectores de corrimientos tt y ta ########
		tt    = sum(t1)
		ta[i] = tt

	return t0, t1, ta, l, CC, N, M
	
def aligne(ta, X2, lv, pcf):
	"""
	Alinea temporalmente X2 de acuerdo a los coeficientes ta de a lv muestras
	"""

	############ Crossover ###############
	lon = lv*pcf                        # Longitud del CO
	t = np.arange(0, np.pi, np.pi/lon)  # Variable independiente del CO
	fad = 0.5 + 0.5*np.sin(t+np.pi/2.0) # Funcion del CO
	
	############ Generacion de arrays para el ciclo #################
	X2 = np.append(X2, np.zeros(abs(int(ta[-1]))))  # Extiende X2 en caso que el resultado final sea mas largo que el original
	X2 = np.append(X2, np.zeros(lv-len(X2)%lv)) 	# Extiende X2 para ser multiplo de lv
	Xf = np.zeros(len(X2)) 							# Genera el array Xf donde se guardan las muestras del archivo final

	############# Alineacion de la ultima ventana ##################
	if ta[-1]>0:
		Xf[len(Xf)-lv:len(Xf)-lv+len(X2[int((len(ta)-1)*lv):int((len(ta)-1)*lv)+lv])] = X2[int((len(ta)-1)*lv):int((len(ta)-1)*lv)+lv]
	
	else:
		Xf[int((len(ta)-1)*lv+ta[-1]):int((len(ta)-1)*lv+ta[-1])+len(X2[int((len(ta)-1)*lv):])] = X2[int((len(ta)-1)*lv):]
	
	############# Alineacion de las ventanas centrales #############
	for i in range(1,len(ta)-2):		
		# print(ta[i])
		# Evaluar caso donde ta[1]>lv
		if ta[i]>=0:
			Xf[int(i*lv-lon/2.0):int(i*lv+lon/2.0)] = X2[int(i*lv-ta[i-1]-lon/2.0):int(i*lv-ta[i-1]+lon/2.0)]*fad+X2[int(i*lv-ta[i]-lon/2.0):int(i*lv-ta[i]+lon/2.0)]*(1-fad) #Banda de transicion
			Xf[int(i*lv+lon/2.0):int(i*lv+lv+ta[i]-lon/2.0)] = X2[int(i*lv-ta[i]+lon/2.0):int(i*lv+lv-lon/2.0)]

		elif all([ta[i]<0,i*lv+ta[i]+lon/2.0>0]):
			Xf[int(i*lv-lon/2.0):int(i*lv+lon/2.0)] = X2[int(i*lv-ta[i-1]-lon/2.0):int(i*lv-ta[i-1]+lon/2.0)]*fad+X2[int(i*lv-ta[i]-lon/2.0):int(i*lv-ta[i]+lon/2.0)]*(1-fad) # Banda de transicion
			Xf[int(i*lv+ta[i]+lon/2.0):int(i*lv+lv-lon/2.0)] = X2[int(i*lv+lon/2.0):int(i*lv+lv-ta[i]-lon/2.0)]
	
	############# Alineacion de la primer ventana ##################
	if ta[0]>0: # mover senal desplazada hacia la derecha
		Xf[:int(ta[0]+lv)] = np.append(np.zeros(int(abs(ta[0]))), X2[:lv])
	
	elif all([ta[0]<0,abs(ta[0])<lv]): 		# mover senal desplazada hacia la izquierda
		Xf[:int(lv+ta[0])] = X2[int(abs(ta[0])):lv]

	elif all([ta[0]<0,abs(ta[0])>lv]):
		Xf = np.append(X2[:int(abs(t[0]))],Xf)
	
	return Xf

def expExcel(fileName, sheetNum, row0, col0, x):
	"""
	Exporta un array a un archivo .xls existente.
	Input:
		fileName: ubicacion del archivo
		sheetNum: numero de la hoja a modificar
		row0: fila inicial
		col0: columna inicial
		x: array 
	"""
	rb = xlrd.open_workbook(fileName)
	wb = copy(rb)
	w_sheet = wb.get_sheet(sheetNum)
		
	r, c = np.shape(x)
	
	for i in np.arange(r):
		for j in np.arange(c):
			w_sheet.write(i + row0, j + col0, x[i][j])
		
	wb.save(fileName)
