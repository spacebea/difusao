#!/usr/bin/python

from math import erf, exp, sqrt, pi, cos
import matplotlib.pyplot as plt
from numpy import zeros

Co =      0.25
Cs =      0.60 
tempo1  =  1.0*3600.00
tempo2  =  2.0*3600.00
tempo3  = 10.0*3600.00
tempo4  = 20.0*3600.00
Temp  =    1123.15
R  =    8.314
Do =    6.2e-07
Qd =    80000.0
L  =    0.020
n  =    41

axis_x1= zeros(n+1)
axis_y1 = zeros(n+1)
axis_x2= zeros(n+1)
axis_y2 = zeros(n+1)
axis_x3= zeros(n+1)
axis_y3 = zeros(n+1)
axis_x4= zeros(n+1)
axis_y4 = zeros(n+1)

D = Do * exp(-Qd/(R*Temp))

print('Calculando do coeficiente de difusao D = %g ' % (D))

dx = L/float(n-1)

print('Co = %g ' %(Co))
print('Cs = %g ' %(Cs))
print('tempo1 = %10.2f '  %(tempo1))
print('tempo2 = %10.2f '  %(tempo2))
print('tempo3 = %10.2f '  %(tempo3))
print('tempo4 = %10.2f '  %(tempo4))
print('Temp = %g ' %(Temp))
print('R = %g ' %(R))
print('Do = %.2e ' %(Do))
print('Qd = %.2e ' %(Qd))
print('L = %.4f ' %(L))
print('n = %i ' %(n)) 
print('dx = %g m' %(dx))


for i in range(1, n+1, 1):
    X = (i-1)*dx
    Arg = X/(2.0*sqrt(D*tempo1))
    Cx = Cs - (Cs-Co)*erf(Arg)
#print ('%.4f     %.4f' %(X, Cx))
    axis_x1[i] = X
    axis_y1[i] = Cx


for j in range(1, n+1, 1):
    X = (j-1)*dx
    Arg = X/(2.0*sqrt(D*tempo2))
    Cx = Cs - (Cs-Co)*erf(Arg)
#print ('%.4f     %.4f' %(X, Cx))    
    axis_x2[j] = X
    axis_y2[j] = Cx


for k in range(1, n+1, 1):
    X = (k-1)*dx
    Arg = X/(2.0*sqrt(D*tempo3))
    Cx = Cs - (Cs-Co)*erf(Arg)
#print ('%.4f     %.4f' %(X, Cx))
    axis_x3[k] = X
    axis_y3[k] = Cx

for w in range(1, n+1, 1):
    X = (w-1)*dx
    Arg = X/(2.0*sqrt(D*tempo4))
    Cx = Cs - (Cs-Co)*erf(Arg)
#print ('%.4f     %.4f' %(X, Cx))
    axis_x4[w] = X
    axis_y4[w] = Cx


plt.title('Difusao do C no Aço ASTM 1025')
plt.plot(axis_x1, axis_y1, 'b-')
plt.plot(axis_x2, axis_y2, 'r-')
plt.plot(axis_x3, axis_y3, 'g-')
plt.plot(axis_x4, axis_y4, 'y-')
#plt.ylim(ymin=0)
#plt.xlim(xmin=0)
plt.axis([0.0,0.02,0.0,0.6])
plt.legend(('1 hora erf', '2 horas erf', '10 horas erf', '20 horas erf'))
plt.ylabel('Concentracao em Percentagem Peso C, [%p C]')
plt.xlabel('Distancia X [m]')
plt.show();

