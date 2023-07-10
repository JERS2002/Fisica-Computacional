import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Leer datos del archivo .dat
filename = "veltot1.dat"
data = np.loadtxt(filename)

# Calcular histograma con 30 bins
hist, bin_edges = np.histogram(data, bins=30)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Normalize frequencies
total_data_points = len(data)
normalized_hist = hist / total_data_points

# Leer datos del segundo archivo .dat
filename2 = "vel0tot1.dat"
data2 = np.loadtxt(filename2)

# Calcular histograma con 30 bins y frecuencias normalizadas
hist2, bin_edges2 = np.histogram(data2, bins=10)
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Normalize frequencies
total_data_points2 = len(data2)
normalized_hist2 = hist2 / total_data_points2

def maxwell_boltzmann(x, a, scale):
    return a * np.sqrt(2/np.pi) * (x**2) * np.exp(-x**2 / (2*scale**2))

# Ajustar histograma a la distribución de Maxwell-Boltzmann
params, _ = curve_fit(maxwell_boltzmann, bin_centers, normalized_hist)
a, scale = params

# Generar la curva ajustada (Maxwell-Boltzmann)
maxwell_boltzmann_curve = maxwell_boltzmann(bin_centers, a, scale)

# Definir función de ajuste para la gaussiana de equipartición
T=0.71
def equipartition_maxwell(x, T):
    return 0.15* np.sqrt(2/np.pi)*(1/T)**1.5 * x**2 *np.exp(-x**2 / (2 * T))

# Generar la curva de equipartición
equipartition_curve = equipartition_maxwell(bin_centers, T)

# Graficar histograma y curva ajustada
plt.bar(bin_centers, normalized_hist, width=np.diff(bin_edges)-0.03, align='center', alpha=0.5, label='t=20-50')
plt.bar(bin_centers2, normalized_hist2, width=np.diff(bin_edges2)-0.03, align='center', alpha=0.5, label='t=0')
plt.plot(bin_centers,maxwell_boltzmann_curve, color='red', label='Ajuste')
plt.plot(bin_centers, equipartition_curve, color='green', label='Equipartición')
plt.xlabel('$v$')
plt.ylabel('P($v$)')
plt.title('Distribución del módulo de la velocidad')
plt.legend()

# Guardar la gráfica como imagen PNG
output_filename = "histogramavtot1.png"
plt.savefig(output_filename)
plt.show()
# Mostrar mensaje de éxito
print("Parámetros del ajuste:")
print("Amplitud (a):", a)
print("Escala (scale):", scale)