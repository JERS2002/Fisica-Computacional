# NOTA PRELIMINAR: Consultar la ayuda de gnuplot para aquello que podamos necesitar.
# Ejemplo: Si necesitamos ayuda sobre el estilo de línea teclearemos en nuestro terminal de gnuplot 'help linestyle'

# ==============================  Ejemplo de representación: Péndulo con ajuste ==============================================

# A partir de aquí se muestra un ejemplo (periodo de un pendulo simple al cuadrado frente a su longitud) donde 
# se representará un fichero de puntos experimentales con sus errores, 
# añadiendo además una curva de ajuste (lineal). Mostramos diferentes estilos, tamaños o colores para las líneas y los puntos.

# ------------- Fichero de datos (medidos) de entrada --------------------
# Fichero de entrada de datos con cuatro columnas: x (longitud), y (periodo al cuadrado), error x, error y 
# El fichero 'datos.dat' se encontrará en la misma carpeta que este script
fichero = "C:/Users/Ernesto/Desktop/Sistemasolar/energia.dat"

# --- Nombre del fichero de salida, correspondiente a la gráfica final en formato pdf ------
salida = "C:/Users/Ernesto/Desktop/Sistemasolar/Graficaenergia.png"

# ----------- Nombre de las etiquetas de los ejes en el gráfico final ----------------------
set xlabel "t (dias)" # para la longitud  
set ylabel "E=T+V" # para el cuadrado del periodO

## ----------------- Ejes logarítmicos -----------------------------------------------------
# Con los siguientes comandos (desactivados en el ejemplo) podemos imponer ejes logarítmicos
# set logscale x
# set logscale y

# ------------------------------------   Ejemplos de ajuste ---------------------------------------------------------------------- 
# Primero definimos la funcion de ajuste f(x), en el primer caso una recta, 
# para después indicar los parámetros de ajuste (en el primer caso a y b), indicando además que ajustaremos el fichero de entrada. 
# El ajuste genera un fichero en la carpeta en la que nos encontramos llamado 'fit.log', en el que se recogen todos los
# datos del ajuste.
#f(x)=a+b*x # En sentido estricto, deberíamos usar un ajuste lineal para el péndulo: relación T^2 con L
#fit f(x) fichero using 1:2:4 yerror via a,b # Atención: debemos normalizar los errores de a y b de gnuplot
# fit f(x) fichero via a, b
# Otro ejemplo (desactivado) con un polinomio de segundo grado y 3 parámetros de ajuste:
# f(x)=a+b*x+c*x**2
# fit f(x) fichero via a, b, c

# ----------------- Representación del fichero de puntos experimentales -----------------------------------------------------------
# Imponemos los rangos de los ejes, la leyenda ('Valor Medido'), además de pedir que se muestren los errores (w xyerror)
# Elegimos ademos el color de los puntos (rojo en el ejemplo)
# Elegimos el tamaño de los puntos mediante ps (pointsize): ps 1 en el ejemplo
# Elegimos el tipo de puntos mediante pt (pointtype): pt 7 en el ejemplo (círculos llenos)
plot [0:60000] [0:8e39]  fichero with linespoints    lc rgb "red" ps 0 pt 7,  \

 
 # --------------- Operaciones con columnas: algunos ejemplos (desactivados) -------------------------------------------------------
# 1- Como mostrar la primera columna y la segunda columna multiplicada por 2, sin errores:
# plot [0:1.75][0:14] fichero u 1:(2*$2) title "Valor Medido" lc rgb "red" ps 1 pt 7
# 2- Como mostrar la primera columna elevada a 2 y la raíz cuadrada de la segunda, sin errores:
# plot [0:4][0:4] fichero u ($1**2):(sqrt($2)) title "Valor Medido"  lc rgb "red" ps 1 pt 7
# 3- Como mostrar la primera columna y el logaritmo neperiano de la segunda, sin errores:
# plot [0:1.75][0:7] fichero u ($1):(log($2)) title "Valor Medido"  lc rgb "red" ps 1 pt 7


# ----------------- Representación de la función de ajuste ------------------------------------------------------------------------
# Sobre el plot anterior añadimos la función de ajuste f(x) usando replot
# Primero fijamos el color (negro en el ejemplo) 
#set style lines 1 lt rgb "black" 
#replot f(x) title "Ajuste Lineal" with lines linestyle 1 dt 1 lw 2 # línea continua sin puntos (dt 1) y anchura 3 (lw 3) 
# replot f(x) title "Ajuste Lineal" with lines linestyle 1 dt 2 lw 3 # línea a trazos sin puntos (dt 2) y anchura 3 (lw 3)
# replot f(x) title "Ajuste Lineal" with lines linestyle 1 dt 3 lw 2 # línea 'punteada' sin puntos (dt 3) y anchura 2 (lw 2)  
# replot f(x) title "Ajuste Lineal" with linespoints linestyle 1 lw 1 # línea continua con puntos (linespoints) y anchura 1 (lw 1)

## ----------------- Leyenda ---------------------------------------------------------------------------
# Posicionamos la leyenda arriba y a la izquierda para que no se solape con los datos en la imagen final
set key bottom
set key right
set key Left
set key nobox
set key noautotitles
## ----------------- Terminal de salida ----------------------------------------
# Imponemos el terminal de salida como pdf, con color y tipo de letra (y tamaño)
set terminal png 
set output salida
replot
