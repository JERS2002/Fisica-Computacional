import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

# Leer datos del archivo .dat
filename = "velx4.dat"
data = np.loadtxt(filename)

# Calcular histograma con 30 bins
hist, bin_edges = np.histogram(data, bins=30)
bin_centers = (bin_edges[:-1] + bin_edges[1:]) / 2

# Normalize frequencies
total_data_points = len(data)
normalized_hist = hist / total_data_points

# Leer datos del segundo archivo .dat
filename2 = "vel0x4.dat"
data2 = np.loadtxt(filename2)

# Calcular histograma con 30 bins y frecuencias normalizadas
hist2, bin_edges2 = np.histogram(data2, bins=10)
bin_centers2 = (bin_edges2[:-1] + bin_edges2[1:]) / 2

# Normalize frequencies
total_data_points2 = len(data2)
normalized_hist2 = hist2 / total_data_points2

# Definir la función de ajuste (Gaussiana)
def gaussian(x, a, std):
    return a * np.exp(-(x)**2 / (2 * std**2))

# Ajustar histograma a la Gaussiana
params, _ = curve_fit(gaussian, bin_centers, normalized_hist)
a, std = params

# Generar la curva ajustada (Gaussiana)
gaussian_curve = gaussian(bin_centers, a, std)

# Definir función de ajuste para la gaussiana de equipartición
T=8
stdequi=np.sqrt(T)

def equipartition_gaussian(x, stdequi):
    return 0.55/ (stdequi * np.sqrt(2 * np.pi)) * np.exp(-x**2 / (2 * stdequi**2))

# Generar la curva de equipartición
equipartition_curve = equipartition_gaussian(bin_centers, stdequi)

# Graficar histograma y curva ajustada
plt.bar(bin_centers, normalized_hist, width=np.diff(bin_edges)-0.03, align='center', alpha=0.5, label='t=20-50')
plt.bar(bin_centers2, normalized_hist2, width=np.diff(bin_edges2)-0.03, align='center', alpha=0.5, label='t=0')
plt.plot(bin_centers, gaussian_curve, color='red', label='Ajuste Gaussiana')
plt.plot(bin_centers, equipartition_curve, color='green', label='Gaussiana equipartición')
plt.xlabel('$v_x$')
plt.ylabel('P($v_x$)')
plt.title('Distribución de la componente x de la velocidad')
plt.legend()

# Guardar la gráfica como imagen PNG
output_filename = "histogramavx4.png"
plt.savefig(output_filename)
plt.show()
# Mostrar mensaje de éxito
print("Gráfica del histograma y ajuste Gaussiano guardada como", output_filename)
print("Parámetros del ajuste:")
print("Amplitud (a):", a)
print("Desviación estándar (std):", std)