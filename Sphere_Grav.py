xg = float(input('Digite a coordenada x do gravímetro:'))
yg = float(input('Digite a coordenada y do gravímetro:'))
zg = float(input('Digite a coordenada z do gravímetro:'))
xc = float(input( 'Digite a coordenada x do centro da esfera'))
yc = float(input( 'Digite a coordenada y do centro da esfera'))
zc = float(input('Digite a coordenada z do centro da esfera'))
r =  ((xg - xc)**2 + (yg - yc)**2 + (zg - zc)**2)**(1/2)
G = 6.674*(10**(-11))
if r**2 == 0:
    print('Erro de coordenadas')
R = float(input('Digite o Raio da esfera:'))
p = float(input('Digite a densidade da esfera'))
pi = 3.1415
m = (4*pi*p*R**3)/3
g = -(G*m)/r**3
print('g é igual a {}'.format(g))
