from scipy.io.wavfile import write, read
import numpy as np
import matplotlib.pyplot as plt
import scipy.signal as signal
import tesisFun4_4


#################### Variables ####################

t    = 11 		# longitud de archivos a comparar
n    = 1001    # numero de muestras rta al impulso
f1   = 100
f2   = 10000  # fc superior
lv   = 200    # longitud de la ventana
g    = 0.08    # umbral del gate
p    = 0.1    # longitud de ventana temporal como porcentaje de lv
pcf  = 0.05 	# porcentaje de crossfade
mv   = 2		# velocidad maxima 
CCM  = 0.06  # Coeficiente de correlacion minima
v    = np.arange(-0.05, 0.05, 0.02) 
vCCM = v+CCM
vG   = v+g
vCCMM, vGM = np.meshgrid(vCCM, vG)
# fig = plt.figure()
# ax = fig.gca(projection='3d')

#################### Seleccion de archivos ####################

Senal    = "B"
nDeSenal = "2"

fs, x1 = read('AudiosFinales/Senales'+Senal+'/'+Senal+'des'+nDeSenal+'.wav')
fs, x2 = read('AudiosFinales/Senales'+Senal+'/'+Senal+'ref'+nDeSenal+'.wav')

#################### Recortar archivos ####################

y1, y2 = tesisFun4_4.trim(x1, x2, fs, t)

#################### Normalizacion ####################

X1 = tesisFun4_4.norm(y1)
X2 = tesisFun4_4.norm(y2)

#################### Filtrado ####################

X1 = tesisFun4_4.filtroPB(n,f1,f2,fs,X1,lv)
X2 = tesisFun4_4.filtroPB(n,f1,f2,fs,X2,lv)

################### Eleccion de la ventana inicial##############

X1,X2,u,y1,y2 = tesisFun4_4.vinicial(X1,X2,lv,y1,y2)
print(u/200.0)

################### Correlacion ####################
c = 0
num = 0
y1 = y1/float(max(y1))
corF = np.zeros(len(vCCM)*len(vG))
corF = np.resize(corF,(len(vCCM),len(vG)))
print(corF)
for i in range(len(vCCM)):
	for j in range(len(vG)):
		t01, t11, ta1, l, CC, N, M = tesisFun4_4.cor2(X1, X2, lv, p, round(vG[j],2), mv, fs, round(vCCM[i],2), u)
		Xf1 = tesisFun4_4.aligne(ta1, y2, lv, pcf)
		xm2 = Xf1/float(abs(max(Xf1)))
		xm1 = y1/float(abs(max(y1)))
		
		if len(xm1)>len(xm2):
			xm1=xm1[:len(xm2)]
		elif len(xm2)>len(xm1):
			xm2=xm2[:len(xm1)]

		RMS1a = float(np.sqrt(np.mean(xm1**2)))
		RMS2a = float(np.sqrt(np.mean(xm2**2)))
		x1  = xm1/RMS1a
		x2  = xm2/RMS2a

		c1 = signal.correlate(x1, x2, 'full', 'fft')
		c2 = max(c1)/len(x1)
		num += 1
		print num
		corF[i,j]=c2
		if c2>c:
			c   = c2
			ta  = ta1
			Xf  = Xf1
			t0  = t01
			t1  = t11
			g   = round(vG[j],2)
			CCM = round(vCCM[i],2)


################### Alineacion ##################

Xf = tesisFun4_4.aligne(ta, y2, lv, pcf)

################### Generacion de archivo .wav ###################

Xf = Xf/float(max(Xf))
y1 = y1/float(max(y1))
write('AudiosFinales/Senales'+Senal+'/'+Senal+',ref'+nDeSenal+'.wav', fs, Xf)
write('AudiosFinales/Senales'+Senal+'/'+Senal+',cor'+nDeSenal+'.wav', fs, y1)
