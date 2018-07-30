#CRIANDO OS GRAVÍMETROS E SUAS POSIÇÕES
from random import randint
from random import random
import numpy as np
#ng = int(input('Digite o número de gravímetros'))  #VOU COMENTAR ESSE CARA POR ENQUANTO
min = randint(0,10) #Sorteando a posição inicial do Gravímetro.
max = randint(0,10) #Sorteando a posição final do Gravímetro.

if min > max:       #Caso o mínimo e o máximo estiverem trocados.
    min = max
    max = min
else:
    min = min
    max = max
if min == max:      #Caso o mínimo e o máximo forem iguais.
    min = randint(0,max)
    max = randint(min,10)
print(min, max)
x = np.linspace(min, max, num = 10)  #Vetor contendo as posições x do gravímetro.
coord = list(range(10))   #Estou supondo 10 Gravímetros pois me parece um bom número para conferir o funcionamento do programa.
for n in range(0,10):
    coord[n] = (x[n],0)   #A coordenada z vai ser sempre 0, pois os Gravímetros estão na superfície.
#print(coord)

#DEFININDO ALGUMAS CONSTANTES
G = 6.674*(10**(-11))
pi = 3.1415

#CRIANDO A ESFERA DE DENSIDADE P(RHO) E RAIO R NA POSIÇÃO (Xc,Zc)
R = random()
p = random()  #Não sei se esses valores estão bons, mas acredito que números entre 0 e 1 são valores razoáveis.
xc = randint(min,max)
zc = -10        #Vou deixar esses caras fixos por enquanto, depois posso randomizar as posições deles.
m = (4*pi*p*R**3)/3

#CÁLCULO DAS ANOMALIAS MEDIDAS EM CADA GRAVÍMETRO
r = list(range(10))
g = list(range(10))
for n in range(0,10):
    r[n] = ((xc - x[n])**2 + zc**2)**(1/2)
    if r[n] == 0:
        print('Alguma coordenada está errada')
    g[n] = abs(-(G*m)/r[n]**3)
print('O vetor das anomalias gravimétricas é {}'.format(g))
