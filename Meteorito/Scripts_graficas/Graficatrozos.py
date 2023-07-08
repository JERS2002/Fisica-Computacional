import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
# Leer el archivo de entrada
archivo = 'trozos.dat'  # Ruta del archivo de entrada
datos = np.loadtxt(archivo)

# Extraer los datos de las columnas
x = datos[:, 0]
y1 = datos[:, 1]
y2 = datos[:, 2]

# Crear la figura y los ejes
fig, ax = plt.subplots()

# Graficar las funciones
ax.plot(x, y1, label='Trayectoria Trozo 1')
ax.plot(x, y2, label='Trayectoria Trozo 2')

# Agregar la línea horizontal punteada en y=0.3
ax.axhline(y=0.0166, linestyle='--', color='r', label="Radio terrestre")

# Etiquetas de los ejes y título
ax.set_xlabel('Tiempo (Días)')
ax.set_ylabel('Distancia a la Tierra ($d_{TL}$)')
ax.set_title('Energía de detonación E= 13,49 Gt')

ax.set_ylim(0.0, 15)
# Leyenda
ax.legend()

ax_zoom = inset_axes(ax, width="30%", height="30%", loc="center left")

ax_zoom.plot(x, y1, label='Función 1')
ax_zoom.plot(x, y2, label='Función 2')
ax_zoom.axhline(y=0.0166, linestyle='--', color='r', label="Radio terrestre")
# Establecer los límites del eje x y el eje y para hacer zoom en una zona específica
zoom_x_start = 2.5 # Valor de inicio del zoom en el eje x
zoom_x_end = 2.8# Valor de fin del zoom en el eje x
zoom_y_start = 0.0 # Valor de inicio del zoom en el eje y
zoom_y_end = 0.03 # Valor de fin del zoom en el eje y
ax_zoom.set_xlim(zoom_x_start, zoom_x_end)
ax_zoom.set_ylim(zoom_y_start, zoom_y_end)

ax_zoom.set_xticks([])
ax_zoom.set_yticks([])


# Guardar la gráfica en formato PNG
plt.savefig('TrozosE13,49.png')

# Mostrar la gráfica en pantalla
plt.show()
