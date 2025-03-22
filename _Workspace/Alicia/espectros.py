import numpy as np
import matplotlib.pyplot as plt
from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from astropy.modeling import models
from astropy import units as u
from scipy.interpolate import CubicSpline, interp1d

lines_spectral = np.array([3970, 4102, 4341, 4861, 4471, 4541])
lines_spectral_names = [r'H$_{\epsilon}$', r'H$_{\delta}$', r'H$_{\gamma}$', r'H$_{\beta}$', 'HeI', 'HeII']

# #para que salgan las coordenadas flotantes
# def on_hover(event):
#     if event.inaxes:  # Verifica que el evento esté dentro de un eje
#         x, y = event.xdata_1, event.ydata_1
#         ax.set_title(f'Coordenadas: x={x:.2f}, y={y:.2f}')

# fig, ax = plt.subplots()

# # Establece el evento de "hover" sobre el gráfico
# fig.canvas.mpl_connect('motion_notify_event', on_hover)


#coge los datos de la estrella 1 o 2 de un archivo de texto
#delimeter indica cual es el delimitado entre una columna y otra
#dtype es el tipo de datos que coge
data_1=np.genfromtxt('C:/Users/Pablo/OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/starprob1.dat', dtype=float)
data_2=np.genfromtxt('C:/Users/Pablo/OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/starprob2.dat', dtype=float)
x_1=data_1[:,0]
y_1=data_1[:,1]

x_2=data_2[:,0]
y_2=data_2[:,1]


# selected_rows_1=[]
# selected_rows_2=[]

# min_value=3900
# max_value=5100

# for row in data_1:
#     if min_value <= row [0] <= max_value:
#         selected_rows_1.append(row)

# selected_rows_1=np.array(selected_rows_1)

# for row in data_2:
#     if min_value <= row [0] <= max_value:
#         selected_rows_2.append(row)

# selected_rows_2=np.array(selected_rows_2)




# x_1 = selected_rows_1[:,0]
# y_1 = selected_rows_1[:,1]

# x_2 = selected_rows_2[:,0]
# y_2 = selected_rows_2[:,1]





#plotear el espectro sin normalizar
plt.plot(x_1,y_1) #scatter es puntitos y plot son lineas
plt.xlabel('$\lambda$ ($\t{\AA}$)')
plt.ylabel('Flux')
# plt.title('Spectrum Star 1')
plt.figure()

plt.plot(x_2,y_2) #scatter es puntitos y plot son lineas
plt.xlabel('$\lambda$ ($\t{\AA})$')
plt.ylabel('Flux')
# plt.title('Spectrum Star 2')
plt.show()