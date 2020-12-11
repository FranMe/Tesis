import numpy as np
import matplotlib.pyplot as plt
from scipy.io.wavfile import write, read
import resampy

n = 1		 # numero de velocidades
tr = 2		 # numero de audio: 1.sine-200 2. vozMm  
A = 1        # atenuador
a = 0.5		 # coeficiente de atenuacion
c = 343.0    # Velocidad del sonido en metros por segundo

if tr == 1:                                   # Senal de entrada
	t = 1 
	fs = 92000.0
	f = 300
	x = np.arange(0, t, 1/fs)
	y1 = np.sin(2*np.pi*f*x)
elif tr == 2:
	t = 7.0     # segundos		
	y1 = ('../Sonidos/Prueba/vozPruebaFinal.wav')
	fs, y1 = read(y1)      # Importa audios
	fs = float(fs)          
	y1 = y1/float(max(y1)) # Normalizacion de senal de entrada

if len(y1)>fs*t: 
	y1 = y1[0:int(fs*t)]	   # Acorta la longitud de los audios

# Resampleo
if n == 1:         
	vs = -1.029                           # metros por segundo. vs=0,343 para un desplazamiento de 1 muestra por intervalo de 1000 muestras 
	Nm = np.round(fs*vs/c)
	y2 = resampy.resample(y1, fs, fs+Nm)

else:
	vs1 = 0.686
	vs2 = -0.686
	Nm1 = np.round(fs*vs1/c)
	Nm2 = np.round(fs*vs2/c)
	y11 = resampy.resample(y1, fs, fs+Nm1)
	y12 = resampy.resample(y1, fs, fs+Nm2)
	y11 = y11[:int(round(len(y11)/2))]
	y12 = y12[int(round(len(y12)/2)):len(y12)]
	y2 = np.append(y11, y12)

# atenuacion
if A==1:
	At1 = np.arange(1, a, -(1-a)/(float(len(y1))))
	y1 = At1*y1
	At2 = np.arange(1, a, -(1-a)/(float(len(y2))))
	y2 = At2*y2
	
# generacion de archivos .wav
if tr==1: 
	write('../Sonidos/Prueba/sine-1000m.wav', fs, y2)
	if A==1:
		write('../Sonidos/Prueba/sine-200ma.wav', fs, y2)
elif tr==2:
	write('../Sonidos/Prueba/vozMmm.wav', fs, y2)
	if A==1:
		 write('../Sonidos/Prueba/vozMmma.wav', fs, y2)

# lp = 0.05
# pi = 0.545
# y1 = y1[int(pi*fs):int((pi+lp)*fs)]
# y2 = y2[int(pi*fs):int((pi+lp)*fs)]
# # ts = np.arange(pi,pi+lp-1/fs,1/fs)
# # ts = np.arange(pi,pi+lp,1/fs)
# ts = np.arange(pi,pi+lp+1/fs,1/fs)
# plt.figure(1)
# plt.title('f = 300 Hz')
# plt.plot(ts[:len(y1)],y1, 'blue')
# plt.plot(ts[:len(y2)],y2, 'red')
plt.plot(y1)
plt.plot(y2)
plt.xlabel('Segundos')
plt.ylabel('Amplitud')

plt.show()