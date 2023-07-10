import matplotlib.pyplot as plt
import numpy as np

# Leer datos del archivo
filename = "temperaturacuadrado.dat"
data = np.loadtxt(filename)

# Extraer columnas x e y
x = data[:, 0]
y = data[:, 1]

# Graficar la función
plt.plot(x, y)
plt.xlabel('Tiempo')
plt.ylabel('Temperatura')
plt.title('Evolución temporal de la temperatura')

# Mostrar la gráfica por pantalla
plt.show()

# Guardar la gráfica como imagen PNG
output_filename = "tempcuadrado.png"
plt.savefig(output_filename)

# Mostrar mensaje de éxito
print("Gráfica guardada como", output_filename)
