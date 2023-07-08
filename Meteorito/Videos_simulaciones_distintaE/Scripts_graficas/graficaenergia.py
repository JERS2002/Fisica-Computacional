import matplotlib.pyplot as plt

def representar_funcion(x, y):
    plt.plot(x, y)
    plt.xlabel('Días')
    plt.ylabel('MJ/kg')
    plt.title('Gasto energético acumulado en función del tiempo')
    plt.savefig('graficaenergia.png')  # Guardar la gráfica como archivo PNG

def leer_datos(nombre_archivo):
    x = []
    y = []

    with open(nombre_archivo, 'r') as archivo:
        for linea in archivo:
            valores = linea.split()
            x.append(float(valores[0]))
            y.append(float(valores[1]))

    return x, y

# Nombre del archivo con los datos
nombre_archivo = 'energia.dat'

# Leer datos del archivo
x, y = leer_datos(nombre_archivo)

# Representar datos como una función y guardar la gráfica
representar_funcion(x, y)

