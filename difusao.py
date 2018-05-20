#!/usr/bin/python

from math import erf, exp, sqrt, pi, cos
import matplotlib.pyplot as plt
from numpy import zeros


def tau(d, t, l):
    return (d * t) / (l ** 2)


C0 = 0.25
Cs = 0.60
tempo1 = 1.0 * 3600.00
tempo2 = 2.0 * 3600.00
tempo3 = 10.0 * 3600.00
tempo4 = 20.0 * 3600.00
Temp = 1123.15
R = 8.314
D0 = 6.2e-07
Qd = 80000.0
L = 0.020
n = 41

axis_x1 = zeros(n + 1)
axis_y1 = zeros(n + 1)
axis_x2 = zeros(n + 1)
axis_y2 = zeros(n + 1)
axis_x3 = zeros(n + 1)
axis_y3 = zeros(n + 1)
axis_x4 = zeros(n + 1)
axis_y4 = zeros(n + 1)

D = D0 * exp(-Qd / (R * Temp))
tau1 = tau(D, tempo1, L)
tau2 = tau(D, tempo2, L)
tau3 = tau(D, tempo3, L)
tau4 = tau(D, tempo4, L)

csv1 = ''
csv2 = ''
csv3 = ''
csv4 = ''

print('Coeficiente de difusao D = %g ' % D)

dx = L / float(n - 1)

print('C0 = %g ' % C0)
print('Cs = %g ' % Cs)
print('Tempo 1 = %10.2f ' % tempo1)
print('Tempo 2 = %10.2f ' % tempo2)
print('Tempo 3 = %10.2f ' % tempo3)
print('Tempo 4 = %10.2f ' % tempo4)
print('Temperatura em Kelvin = %g ' % Temp)
print('R = %g ' % R)
print('D0 = %.2e ' % D0)
print('Qd = %.2e ' % Qd)
print('L (comprimento) = %.4f ' % L)
print('Número de pontos = %i ' % n)
print('dx = %g m' % dx)

print('\n1 hora')
print('Posição (mm)\tC_Series\t\tC_Erf')
csv1 += 'Posição (mm);C_Series;C_Erf\n'
for i in range(1, n + 1, 1):
    X = (i - 1) * dx
    Arg = X / (2.0 * sqrt(D * tempo1))
    Cerf = Cs - (Cs - C0) * erf(Arg)
    Xmm = X * 10 ** 3

    p = (L - X) / L
    soma = 0
    for r in range(1, 5000, 1):
        mi = (2 * r - 1) * (pi / 2)
        soma = soma + ((-1) ** (r + 1)) * ((cos(mi * p) * exp((-mi ** 2) * tau1)) / mi)
    Cseries = Cs - (Cs - C0) * 2 * soma

    print('%.4f\t\t\t%.4f\t\t\t%.4f' % (Xmm, Cseries, Cerf))
    csv1 += '%.4f;%.4f;%.4f\n' % (Xmm, Cseries, Cerf)
    axis_x1[i] = X
    axis_y1[i] = Cerf

print('\n2 horas')
print('Posição (mm)\tC_Series\t\tC_Erf')
csv2 += 'Posição (mm);C_Series;C_Erf\n'
for j in range(1, n + 1, 1):
    X = (j - 1) * dx
    Arg = X / (2.0 * sqrt(D * tempo2))
    Cerf = Cs - (Cs - C0) * erf(Arg)
    Xmm = X * 10 ** 3

    p = (L - X) / L
    soma = 0
    for r in range(1, 5000, 1):
        mi = (2 * r - 1) * (pi / 2)
        soma = soma + ((-1) ** (r + 1)) * ((cos(mi * p) * exp((-mi ** 2) * tau2)) / mi)
    Cseries = Cs - (Cs - C0) * 2 * soma

    print('%.4f\t\t\t%.4f\t\t\t%.4f' % (Xmm, Cseries, Cerf))
    csv2 += '%.4f;%.4f;%.4f\n' % (Xmm, Cseries, Cerf)
    axis_x2[j] = X
    axis_y2[j] = Cerf

print('\n10 horas')
print('Posição (mm)\tC_Series\t\tC_Erf')
csv3 += 'Posição (mm);C_Series;C_Erf\n'
for k in range(1, n + 1, 1):
    X = (k - 1) * dx
    Arg = X / (2.0 * sqrt(D * tempo3))
    Cerf = Cs - (Cs - C0) * erf(Arg)
    Xmm = X * 10 ** 3

    p = (L - X) / L
    soma = 0
    for r in range(1, 5000, 1):
        mi = (2 * r - 1) * (pi / 2)
        soma = soma + ((-1) ** (r + 1)) * ((cos(mi * p) * exp((-mi ** 2) * tau3)) / mi)
    Cseries = Cs - (Cs - C0) * 2 * soma

    print('%.4f\t\t\t%.4f\t\t\t%.4f' % (Xmm, Cseries, Cerf))
    csv3 += '%.4f;%.4f;%.4f\n' % (Xmm, Cseries, Cerf)
    axis_x3[k] = X
    axis_y3[k] = Cerf

print('\n20 horas')
print('Posição (mm)\tC_Series\t\tC_Erf')
csv4 += 'Posição (mm);C_Series;C_Erf\n'
for l in range(1, n + 1, 1):
    X = (l - 1) * dx
    Arg = X / (2.0 * sqrt(D * tempo4))
    Cerf = Cs - (Cs - C0) * erf(Arg)
    Xmm = X * 10 ** 3

    p = (L - X) / L
    soma = 0
    for r in range(1, 5000, 1):
        mi = (2 * r - 1) * (pi / 2)
        soma = soma + ((-1) ** (r + 1)) * ((cos(mi * p) * exp((-mi ** 2) * tau4)) / mi)
    Cseries = Cs - (Cs - C0) * 2 * soma

    print('%.4f\t\t\t%.4f\t\t\t%.4f' % (Xmm, Cseries, Cerf))
    csv4 += '%.4f;%.4f;%.4f\n' % (Xmm, Cseries, Cerf)
    axis_x4[l] = X
    axis_y4[l] = Cerf

with open('difusao1.csv', 'w') as file:
    file.write(csv1)

with open('difusao2.csv', 'w') as file:
    file.write(csv2)

with open('difusao3.csv', 'w') as file:
    file.write(csv3)

with open('difusao4.csv', 'w') as file:
    file.write(csv4)

plt.title('Difusao do C no Aço ASTM 1025')
plt.plot(axis_x1, axis_y1, 'b-')
plt.plot(axis_x2, axis_y2, 'r-')
plt.plot(axis_x3, axis_y3, 'g-')
plt.plot(axis_x4, axis_y4, 'y-')
plt.axis([0.0, 0.02, 0.0, 0.6])
plt.legend(('1 hora', '2 horas', '10 horas', '20 horas'))
plt.ylabel('Concentracao em Percentagem Peso C, [%p C]')
plt.xlabel('Distancia X [m]')
plt.show();
