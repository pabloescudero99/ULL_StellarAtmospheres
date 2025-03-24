# Packages used in the project (necessary importing here)
import numpy as np
import matplotlib.pyplot as plt
from astropy.io.ascii import read
from scipy import constants as const

# Constants used in the project
k_B_Julios = const.Boltzmann # J/K
k_B = k_B_Julios * 1e7 # erg / K
k_B_eV = k_B_Julios *6.242e18 # eV/K
c_SI = const.speed_of_light
c = c_SI * 100 # cm/s
h_SI = const.h # J*S
h = h_SI *1e7 # erg*s
R_Rydb = 1.0968e5 #cm-1
e = const.physical_constants['elementary charge'][0] * 10 * c_SI # Fr (Franklin cgs) = dyn^(1/2) * cm
m_e = 9.1094e-28 # g

# # Firsts functions
def pob_e_libres(P_e, T):
    # Calculation of electronic population
    # Pe needed in erg/cm-¿?
    return P_e/(k_B*T)

# Eq Saha
def Saha(j,T,Ne):
    if j == -1: # H-/HI
            U = 1/2
            # Ionization energy H-
            x = 0.755 *1.602e-12 #erg
    if j == 0: # HI/p(HII)
            U = 2/1
            x = 13.6 *1.602e-12  #erg
    return 2.07e-16*Ne*U*T**(-3/2)*np.exp(x/(k_B*T))

# Eq Boltzmann
def Boltzmann(i,T):
    U = 2
    if i == 1:
         x = 0 *1.602e-12 #erg
    if i == 2:
         x = 10.2 *1.602e-12 #erg
    if i == 3:
         x = 12.09 * 1.602e-12 #erg
    g = 2*i**2
    return g/U*np.exp(-x/(k_B*T))

# Eq Rydberg
def Rydberg(n,m):
    # n,m could be np.inf
    term = R_Rydb * (1/n**2 - 1/m**2)
    return 1/term

###################################################################################

# # Other functions

# Electron scattering opacity
def k_e(Ne):
    sigma_e = 6.25e-25
    return sigma_e*Ne

# H- opacity for f-f processes
def k_ff_Hmen(T,l_cm,Pe,N_HI):
    l = l_cm *1e8 # Armstrongs
    theta = 5040/T
    f0 = -2.2763 - 1.6850*np.log10(l) + 0.76661*np.log10(l)**2 - 0.053346*np.log10(l)**3
    f1 = 15.2827 - 9.2846*np.log10(l) + 1.99381*np.log10(l)**2 - 0.142631*np.log10(l)**3
    f2 = -197.789 + 190.266*np.log10(l) - 67.9775*np.log10(l)**2 + 10.6913*np.log10(l)**3 - 0.625151*np.log10(l)**4
    sigma = 1e-26*10**(f0 + f1*np.log10(theta) + f2*np.log10(theta)**2)
    return Pe*sigma*N_HI

# H- opacity for b-f processes
def k_bf_Hmen(l_cm,Ni,nu,T, l0):
    kappas = []
    l = l_cm*1e8 # A
    a_0 = 1.99654
    a_1 = -1.18267e-5
    a_2 = 2.64243e-6
    a_3 = -4.40524e-10
    a_4 = 3.23992e-14
    a_5 = -1.39568e-18
    a_6 = 2.78701e-23
    for i in range(len(nu)):
        if nu[i] >= c/l0:
            sigma = (a_0 + a_1*l[i] + a_2*l[i]**2 + a_3*l[i]**3 + a_4*l[i]**4 + a_5*l[i]**5 + a_6*l[i]**6)*1e-18 #cm**2
            k = sigma*Ni*(1 - np.exp(-h*nu[i]/(k_B*T)))
        else:
            k = 0
        kappas.append(k)
    # print(kappas)
    return np.asarray(kappas)

# Opacity for f-f processes
def k_ff(l,T,nu,Ne,Nk):
    g_ff = 1 + (0.3456/(l*R_Rydb)**(1/3))*(l*k_B*T/(h*c) + 1/2)
    sigma = 3.7e8/(T**(1/2)*nu**3)*g_ff
    return sigma*Ne*Nk*(1 - np.exp(-h*nu/(k_B*T)))

# Opacity for b-f processes
def k_bf(n,l,nu,Ni,T, l1, l2, l3):
    kappas = []
    for i in range(len(nu)):
        g_bf = 1 - 0.3456/(l[i]*R_Rydb)**1/3 *(l[i]*R_Rydb/n**2 - 1/2)
        sigma = 2.815e29/(n**5*nu[i]**3)*g_bf
        if (n==1 and nu[i] >= c/l1):
            k = sigma*Ni*(1 - np.exp(-h*nu[i]/(k_B*T)))
        elif (n==2 and nu[i] >= c/l2):
            k = sigma*Ni*(1 - np.exp(-h*nu[i]/(k_B*T)))
        elif (n==3 and nu[i] >= c/l3):
            k = sigma*Ni*(1 - np.exp(-h*nu[i]/(k_B*T)))
        else:
            k = 0
        kappas.append(k)
    return np.array(kappas)

# Lines absorption
def k_linea(n,m,n_n,n_m):
  if n == 1 and m == 2:
      g_bb = 0.717 # Ref on overleaf: BakerMezel1938. Tab 2
      g_n, g_m = 2*n**2, 2*m**2 # Eq H
      f = 2**5 * g_bb * (1/n**2 - 1/m**2)**(-3) / (3**(3/2) * np.pi * n**5 * m**3)
      sigma = np.pi * e**2 * f / (m_e * c)
      k = sigma * (n_n - g_n * n_m / g_m)
  elif n == 1 and m == 3:
      g_bb = 0.765 # Ref on overleaf: BakerMezel1938. Tab 2
      g_n, g_m = 2*n**2, 2*m**2 # Eq H
      f = 2**5 * g_bb * (1/n**2 - 1/m**2)**(-3) / (3**(3/2) * np.pi * n**5 * m**3)
      sigma = np.pi * e**2 * f / (m_e * c)
      k = sigma * (n_n - g_n * n_m / g_m)
  elif n == 2 and m == 3:
      g_bb = 0.869 - 3/m**3 # Subject notes for Balmer series
      g_n, g_m = 2*n**2, 2*m**2 # Eq H
      f = 2**5 * g_bb * (1/n**2 - 1/m**2)**(-3) / (3**(3/2) * np.pi * n**5 * m**3)
      sigma = np.pi * e**2 * f / (m_e * c)
      k = sigma * (n_n - g_n * n_m / g_m)
  return k