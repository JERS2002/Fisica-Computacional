import matplotlib.pyplot as plt
import numpy as np

# Ruta del archivo .dat
archivo = 'energiasv4.dat'

# Leer los datos del archivo .dat
datos = np.loadtxt(archivo)

# Obtener los valores para cada columna
x = datos[:, 0]
y1 = datos[:, 1]
y2 = datos[:, 2]
y3 = datos[:, 3]

# Crear la figura y los ejes
fig, ax = plt.subplots()

# Graficar las tres funciones
ax.plot(x, y1, color='red', linewidth=0.6, label='Enegía cinética')
ax.plot(x, y2, color='yellow', linewidth=0.6, label='Energía potencial')
ax.plot(x, y3, color='orange', linewidth=0.6, label='Energía total')

# Agregar etiquetas y título
ax.set_xlabel('Tiempo')
ax.set_ylabel('Energía')
ax.set_title('Evolución de las energías cinética, potencial y total para $v_0=4$')

# Agregar una leyenda
ax.legend()

# Guardar la gráfica en formato PNG
plt.savefig('Energiav4.png')

# Mostrar la gráfica en pantalla
plt.show()
