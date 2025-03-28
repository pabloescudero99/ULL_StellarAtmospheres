import numpy as np
import matplotlib.pyplot as plt
from specutils import Spectrum1D
from specutils.fitting import fit_generic_continuum
from astropy.modeling import models
from astropy import units as u
from scipy.interpolate import CubicSpline, interp1d

#para que salgan las coordenadas flotantes
# def on_hover(event):
#     if event.inaxes:  # Verifica que el evento esté dentro de un eje
#         x, y = event.xdata, event.ydata
#         ax.set_title(f'Coordenadas: x={x:.2f}, y={y:.2f}')

# fig, ax = plt.subplots()

# # Establece el evento de "hover" sobre el gráfico
# fig.canvas.mpl_connect('motion_notify_event', on_hover)


#coge los datos de la estrella 1 o 2 de un archivo de texto
#delimeter indica cual es el delimitado entre una columna y otra
#dtype es el tipo de datos que coge
data=np.genfromtxt('C:/Users/Pablo/OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/starprob2.dat', dtype=float)
x=data[:,0]
y=data[:,1]

# Datos del espectro completo
x_full=data[:,0]
y_full=data[:,1]
rows_full=[]
min_value_full=3540
max_value_full=7400

for row in data:
    if min_value_full <= row [0] <= max_value_full:
        rows_full.append(row)

rows_full=np.array(rows_full)
x_full=rows_full[:,0]
y_full=rows_full[:,1]



#el for y el if sirve para poner esto pero en el visual studio no coge las &.
#cogemos solo este intervalo para normalizar porque las estrellas de referencia van de 3900 a 5000.
# x=x[(x>3900)&(x<5000)]
# y=y[(x>3900)&(x<5000)]

selected_rows=[]

min_value=3900
max_value=5100
# max_value=5000

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

cont_onda= (3742.1,3781.7,3825.8,3859.1,3903.2,3950,3986.9,4014.8,4051.7,4081.4,4113.8,4149.8,4182.2,4215.5,4256.9,4297.4,4358.6,
            4402.7,4456.7,4535,4579.1,4658.3,4759.1,4822.1,4909.4,4928.3,4997.6,5021,5147.9,5184.8,5224.4,5286.5,5435.9,5516,5603.3,5656.4,
            5714.9,5804.9,5880.5,5948.9,5985.8,6076.7,6146.9,6249.5,6328.7)

cont_onda = [wave for wave in cont_onda if min(x) <= wave <= max(x)]

# Ajuste del continuum con una spline cúbica utilizando las longitudes de onda reales
cs = CubicSpline(x, y)

# Evaluar la spline en los puntos de longitud de onda 'continuum_wavelengths'
cont_flux = cs(cont_onda)

interp_func = interp1d(cont_onda, cont_flux, kind='cubic', fill_value='extrapolate')
cont_flux_interp = interp_func(x)

spectrum = Spectrum1D(spectral_axis=x*u.AA, flux=y*u.Jy)

normalized_spectrum=y/cont_flux_interp

lines_spectral = np.array([3970, 4102, 4341, 4861, 4471, 4541, 4026,4009])
lines_spectral2= np.array([4552, 4128, 4481, 4089, 3995, 3934, 4267])

lines_spectral_names = [r'H$_{\epsilon}$', r'H$_{\delta}$', r'H$_{\gamma}$', r'H$_{\beta}$', 'HeI', 'HeII', '', 'He I']
lines_spectral_names2 = ['Si III', 'Si II', 'Mg II', 'Si IV', 'N II', 'Ca II', 'C II']



# #grafciamos
# # plt.figure()
# # sin normalizar y con la linea del fit
# plt.plot(x, y, label='Star 2')
# plt.plot(x, cont_flux_interp, color = 'darkorange', label='Continuum fitting')
# plt.xlabel('$\lambda$ ($\t{\AA}$)')
# plt.legend(loc='upper right')
# plt.ylabel("Flux")
# plt.show()

# # plt.figure()
# b0v=np.genfromtxt('C:/Users/Pablo\OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/espectros_referencia\ALS_18929_B0V.dat')
# x_b0v=b0v[:,0]
# y_b0v=b0v[:,1]
b5v=np.genfromtxt('C:/Users/Pablo\OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/espectros_referencia\HD3369_B5V.dat')
x_b5v=b5v[:,0]
y_b5v=b5v[:,1]
# b5i=np.genfromtxt('C:/Users/Pablo\OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/espectros_referencia\HD13267_B5Ia.dat')
# x_b5i=b5i[:,0]
# y_b5i=b5i[:,1]
# b0i=np.genfromtxt('C:/Users/Pablo\OneDrive - Universidad de La Laguna/Master/Cuatri1/AtmosferasEstelares/Entregables/Entregable1/espectros_referencia\HD150898_B0Ib.dat')
# x_b0i=b0i[:,0]
# y_b0i=b0i[:,1]


plt.plot(x,normalized_spectrum, color="magenta", label='Star 2')
# plt.plot(x_b0v,y_b0v,ls='dotted',color ='b', label='ALS_18929_B0V')
# plt.plot(x_b0i,y_b0i,ls='dotted', color ='b', label ='HD150898_B0Ib')
# plt.plot(x_b5i,y_b5i,ls='dotted', color = 'b', label  ='HD13267_B5Ia')
# plt.plot(x_b5v,y_b5v,ls='dotted', color='b', label='HD3369_B5V')
for i in range(0, np.size(lines_spectral)): 
        plt.axvline(lines_spectral[i], 0, 1, color = 'orange', alpha = 0.3)
        plt.text(lines_spectral[i] + 0.5, max(y_b5v), lines_spectral_names[i], fontsize = 12)
for j in range(0, np.size(lines_spectral2)):
        plt.axvline(lines_spectral2[j], 0, 1, color = 'green', alpha = 0.3)
        plt.text(lines_spectral2[j] + 0.5, min(y_b5v), lines_spectral_names2[j], fontsize = 12) #normalized_spectrum

plt.xlabel('$\lambda$ ($\t{\AA}$)')
plt.ylabel('Normalized flux')
plt.xlim(3900, 5100)
plt.legend(loc='lower right')
# plt.title('zoom C II')

plt.show()