import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import linregress

# Leer datos del archivo
filename = "T-P.txt"
data = np.loadtxt(filename)

# Extraer columnas x e y
x = data[:, 0]
y = data[:, 1]

# Realizar regresión lineal
slope, intercept, r_value, p_value, std_err = linregress(x, y)

# Calcular valores de la recta de regresión
x_regression = np.linspace(np.min(x), np.max(x), 100)
y_regression = slope * x_regression + intercept

# Calcular errores en la pendiente e intercepto
sum_x_sq = np.sum((x - np.mean(x)) ** 2)
slope_error = std_err * np.sqrt(1 / len(x) + np.mean(x) ** 2 / sum_x_sq)
intercept_error = std_err * np.sqrt(np.sum(x ** 2) / (len(x) * sum_x_sq))


# Graficar puntos y recta de regresión
plt.scatter(x, y)
plt.plot(x_regression, y_regression, color='red', label='Regresión lineal')
plt.xlabel('Temperatura')
plt.ylabel('Presión')
plt.title('Ecuación de Estado')
plt.legend()

# Mostrar la gráfica por pantalla
plt.show()

# Guardar la gráfica como imagen PNG
output_filename = "ajusteTP.png"
plt.savefig(output_filename)

# Mostrar mensaje de éxito
print("Gráfica de la regresión lineal guardada como", output_filename)

# Imprimir los valores de los parámetros de la regresión lineal
print("Parámetros de la regresión lineal:")
print("Pendiente (slope):", slope)
print("Error en la pendiente (slope_error):", slope_error)
print("Intercepto (intercept):", intercept)
print("Error en el intercepto (intercept_error):", intercept_error)
print("Coeficiente de correlación (r-value):", r_value)
print("Valor-p (p-value):", p_value)
print("Error estándar (std_err):", std_err)
