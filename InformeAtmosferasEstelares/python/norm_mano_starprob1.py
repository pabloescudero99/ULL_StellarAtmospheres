import numpy as np
import matplotlib.pyplot as plt
from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from astropy.modeling import models
from astropy import units as u
from scipy.interpolate import CubicSpline, interp1d

lines_spectral = np.array([3934, 3970, 4045, 4102, 4226, 4300, 4341, 4861, 4471, 4541, 5167, 5172, 5183])
lines_spectral_names = ['Ca II', r'H$_{\epsilon}$', 'FeI', r'H$_{\delta}$', 'Ca I/II','G b', r'H$_{\gamma}$', r'H$_{\beta}$', 'HeI', 'HeII', '', 'Mg', '']
#lineas de magnesio para 5200 A: 5167, 5172, 5183

#para que salgan las coordenadas flotantes
# def on_hover(event):
#     if event.inaxes:  # Verifica que el evento esté dentro de un eje
#         x, y = event.xdata, event.ydata
#         ax.set_title(f'Coordenadas: x={x:.2f}, y={y:.2f}')

# fig, ax = plt.subplots()

# # Establece el evento de "hover" sobre el gráfico
# fig.canvas.mpl_connect('motion_notify_event', on_hover)

#coge los datos de la estrella 1 o 2 de un archivo de texto
data=np.genfromtxt('C:/Users/alich/ULL/1master/Atmosferas_estelares/entregable1/starprob1.dat', dtype=float)
x=data[:,0]
y=data[:,1]

#el for y el if sirve para poner esto pero en el visual studio no coge las &.
#cogemos solo este intervalo para normalizar porque las estrellas de referencia van de 3900 a 5000.
# x=x[(x>3900)&(x<5000)]
# y=y[(x>3900)&(x<5000)]

selected_rows=[]
min_value=3900
max_value=5200

for row in data:
    if min_value <= row [0] <= max_value:
        selected_rows.append(row)
selected_rows=np.array(selected_rows)

x = selected_rows[:,0]
y = selected_rows[:,1]

# #plotear el espectro sin normalizar
# plt.plot(x,y) #scatter es puntitos y plot son lineas
# plt.xlabel('Longitud de onda (Amstrongs)')
# plt.ylabel('no se aun')
# plt.title('Espectro estrella 1')
# plt.grid(True)
# plt.show()

#para hacer el fit hay que hacer la lina a mano:
#cont_onda= (3901, 3913, 4003, 4038, 4088.7, 4317.17, 4504.9, 4795.3, 4894.9, 4943.5, 4963.5, 4974.08, 5033, 5046.02, 5061.03, 5093)
#el continuo de onda para 5200 A
cont_onda= (3901, 3913, 4003, 4038, 4088.7, 4317.17, 4504.9, 4795.3, 4894.9, 4943.5, 4963.5, 4974.08, 5033, 5046.02, 5061.03, 5093, 5118.14, 5199)
cont_onda = [wave for wave in cont_onda if min(x) <= wave <= max(x)]

# Ajuste del continuum con una spline cúbica utilizando las longitudes de onda reales
cs = CubicSpline(x, y)

# Evaluar la spline en los puntos de longitud de onda 'continuum_wavelengths'
cont_flux = cs(cont_onda)
interp_func = interp1d(cont_onda, cont_flux, kind='linear', fill_value='extrapolate')
cont_flux_interp = interp_func(x)
spectrum = Spectrum1D(spectral_axis=x*u.AA, flux=y*u.Jy)
normalized_spectrum=y/cont_flux_interp

# # plt.plot(cont_onda,cont_flux,marker='x',color='red')
# # #graficamos
# # plt.figure()
# #sin normalizar y con la linea del fit
# plt.plot(x, y, label='Star 1')
# plt.plot(x, cont_flux_interp, color = 'darkorange', label='Continuum fitting')
# plt.xlabel('$\lambda$ ($\t{\AA}$)')
# plt.legend(loc='upper right')
# plt.ylabel("Flux")

# plt.show()

#Comparamos la estrella 1 con las estrellas de tipo K
K0V=np.genfromtxt('tipo_k\HD3651_K0V.dat')
x_K0V=K0V[:,0]
y_K0V=K0V[:,1]
K0I=np.genfromtxt('tipo_k\HD12014_K0Ib.dat')
x_K0I=K0I[:,0]
y_K0I=K0I[:,1]
K5V=np.genfromtxt('tipo_k\HD201091_K5V.dat')
x_K5V=K5V[:,0]
y_K5V=K5V[:,1]
K4p5Ib=np.genfromtxt('tipo_k\HD219978_K4p5Ib.dat')
x_K4p5Ib=K4p5Ib[:,0]
y_K4p5Ib=K4p5Ib[:,1]


# #normalizado
plt.plot(x, normalized_spectrum, color="r", label='Star 1')
#plt.plot(x_K4p5Ib,y_K4p5Ib, color='blue', ls='dotted', label='HD219978_K4p5Ib')
for i in range(0, np.size(lines_spectral)): 
        plt.axvline(lines_spectral[i], 0, 1, color = 'dimgrey', alpha = 0.5)
        plt.text(lines_spectral[i] + 0.5, max(normalized_spectrum), lines_spectral_names[i], fontsize = 9)
plt.xlabel('$\lambda$ ($\t{\AA}$)')
plt.ylabel("Normalized flux")
plt.xlim(3900, 5200)
plt.title("Star 1")
#plt.legend(loc="upper right")
plt.grid(False)
plt.show()