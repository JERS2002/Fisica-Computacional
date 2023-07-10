#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <algorithm>
#include <random>

using namespace std;


#define ITER_MAX 30000
#define L 10
#define mindist 0.9
#define N 16
#define h 0.002

//Definimos funciones para calcular aceleracion en tiempo t con ley de Gravitacion de Newton
//Aceleracion eje x
double modulo (int i, int j, double r[N][2])
{
double modulo;
modulo= sqrt(pow(r[i][0]-r[j][0],2) + pow(r[i][1]-r[j][1],2));
return modulo;
}

double distancia (int i, int j, double r[N][2]){
    
    double rij_x, rij_y;

    rij_x= min({abs(r[i][0]-r[j][0]), abs(r[i][0]-r[j][0]-L), abs(r[i][0]-r[j][0]+L) });
    rij_y = min({abs(r[i][1]-r[j][1]), abs(r[i][1]-r[j][1]-L), abs(r[i][1]-r[j][1]+L) });

    
    return sqrt(rij_x*rij_x+rij_y*rij_y);
}

//Calculamos la aceleracion teniendo en cuenta las condiciones de contorno periodicas
void aceleracion(double r[N][2], double a[N][2])
{
 double suma;
 suma=0.0;

 for (int k=0;k<2; k++){

 for (int i=0; i<N; i++){
 for (int j=0; j<N; j++) {
 if (j!=i){
    if (distancia(i, j, r)<3){
        if (abs(r[i][k]-r[j][k])<abs(r[i][k]-r[j][k]-L)&& abs(r[i][k]-r[j][k]) < abs(r[i][k]-r[j][k]+L)){
   suma= suma +((r[i][k]-r[j][k])*(2/(pow(distancia(i,j,r), 14))-1/(pow(distancia(i,j,r), 8))));
        }
   else {
    suma= suma -((r[i][k]-r[j][k])/(abs(r[i][k]-r[j][k]))*(min({abs(r[i][k]-r[j][k]), abs(r[i][k]-r[j][k]-L), abs(r[i][k]-r[j][k]+L) }))*(2/(pow(distancia(i,j,r), 14))-1/(pow(distancia(i,j,r), 8))));
   }
 }
 }
 }
 a[i][k]=24*suma; 
 suma = 0.0;
 }

}
return;
 }




//Funciones para calcular la energía potencial Vi de cada particula y la energia total V
void V_i(double r[N][2], double Vi[N])
{
double suma;
suma=0.0;
for (int i=0; i<N; i++){
 for (int j=0; j<N; j++) {
 if (j!=i){
    if (distancia(i, j, r)<3){
   suma= suma +(1/(pow(distancia(i,j,r),12))-1/(pow(distancia(i,j,r),6))); 
    }
 }
 }
 Vi[i]= 4*suma; 
 suma = 0.0;
 }
return;
}
double Vtot(double Vi[N]){
   double V;
   V=0.0;
    for (int i=0; i<N; i++){
    V= V+Vi[i];
    }
    return V;
}

//Funcion para calcular la energia cinetica Ti de cada particula y la energia total T
void T_i(double v[N][2], double Ti[N]){
for (int i=0; i<N; i++){
    Ti[i]=(v[i][0]*v[i][0]+v[i][1]*v[i][1])/2;
}
return;
}

double Ttot(double Ti[N]){
   double T;
   T=0.0;
    for (int i=0; i<N; i++){
    T= T+Ti[i];
    }
    return T;
}

//Funcion calcular r(h) por Desarrollo en Serie de Taylor
void r_h(double r[N][2], double v[N][2], double a[N][2]){
for (int j=0; j<2; j++ ){
for (int i=0; i<N; i++){
    r[i][j]=fmod((r[i][j]+h*v[i][j]+h*h*0.5*a[i][j]),L);
    if (r[i][j]<0) r[i][j]=r[i][j]+L;
    }
}
return;
}

//Funcion calcular wi= vi+ h/2 *a_i
void w_i(double v[N][2], double a[N][2], double w[N][2]){
 for (int j=0; j<2; j++ ){
for (int i=0; i<N; i++){
    w[i][j]=v[i][j]+0.5*h*a[i][j];
    }
}
return;   
}

// Funcion calcular v_i(h)= w_i(h)+ h/2 *a_i(h)
void v_i(double w[N][2], double a[N][2], double v[N][2]){
    for (int j=0; j<2; j++ ){
for (int i=0; i<N; i++){
    v[i][j]=w[i][j]+ 0.5*h*a[i][j];
    }
}
return;  
}



//Inicializamos los valores de la aceleracion a 0
void inicializar_aceleracion (double a[N][2]){
    for (int i = 0; i < N; i++) {
        for (int j = 0; j < 2; j++) {
a[i][j]=0.0;
}
}
}


//INICIALIZAMOS POSICIONES ALEATORIAS. Estableceremos una distancia minima entre particulas ya que a pequenias distancias las repulsiones se hacen muy grandes
void inicializar_posicion_aleatoria (double r[N][2]){
    random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> pos(0.0, 10.0);
    bool aceptable;
// Generamos posicion de una particula
r[0][0]=pos(gen);
r[0][1]=pos(gen);
    for (int i = 1; i < N; i++) {
        aceptable=false;
        while (aceptable==false)
        {
            aceptable=true;
            r[i][0]=pos(gen);
            r[i][1]=pos(gen);
            for (int j = i-1; j >= 0; j--){
                if (distancia(i,j,r)<mindist) aceptable=false;  //Se compara con la distancia minima aceptable, si no se genera otra posicion para la particula
            }
        }
        
        
}
}

//Inicializar posicion prueba
void inicializar_posicionprueba (double r[N][2]){
    r[0][0]=9.8;
    r[0][1]=3;
    r[1][0]=0.3;
    r[1][1]=3.5;
}

//INICIALIZAMOS VELOCIDADES ALEATORIAS
void inicializar_velocidad (double v[N][2]){
        random_device rd;
    mt19937 gen(rd());
    uniform_real_distribution<double> vel(0.0, 2.0);
    double angulo;
    double modulo;
    modulo=5.0;
    for (int i = 0; i < N; i++) {   
     angulo=M_PI*vel(gen);
    v[i][0]=modulo*cos(angulo);
    v[i][1]=modulo*sin(angulo);

}
}

//Algoritmo VERLET
void Verlet (double r[N][2], double a[N][2], double v[N][2], double w[N][2]){
//Calculamos r_i(t+h) y w_i (t+h)
r_h (r,v,a);
w_i (v,a,w);

// Calculamos a_i(t+h)
aceleracion (r,a);


// Calculamos v_i (t+h)
v_i (w,a,v);
}

//Calcular presion
double presion(double r[N][2], double v[N][2], double a[N][2])
{double fuerza, nuevapos[2];
fuerza=0;

for (int j=0; j<2; j++ ){
for (int i=0; i<N; i++){
    nuevapos[j]=r[i][j]+h*v[i][j]+h*h*0.5*a[i][j];
    if (nuevapos[0]>L || nuevapos[0]<0){
        fuerza=fuerza+abs(2*v[i][0]);
    }
    if (nuevapos[1]>L || nuevapos[1]<0){
        fuerza=fuerza+abs(2*v[i][1]);
        
    }
    }
}
return fuerza;


}





//FUNCION PRINCIPAL
int main ()
{
    double r[N][2], a[N][2], v[N][2], w[N][2], Vi[N], Ti[N], t, TEMP_media, PRESION_media;
    //PASO TEMPORAL
    TEMP_media=0;

//GENERAMOS POSICIONES Y VELOCIDADES INICIALES ALEATORIAS
inicializar_posicion_aleatoria(r);
inicializar_velocidad(v);
//Evaluamos aceleracion en t=0
aceleracion (r,a);
// Escribimos la posicion inicial en el archivo donde guardaremos las posiciones
    ofstream file("posiciones.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "posiciones.dat" << "para leer" << endl;
    }
    for (int i = 0; i < N; i++) {

     file << r[i][0] << ", " << r[i][1] << endl;
     }
     file << endl; 

//FICHERO AUXILIAR
ofstream fich("datos.dat");
    if (!fich) {
        cerr << "No se pudo abrir " << "datos.dat" << "para leer" << endl;
    }
    for (int i = 0; i < N; i++) {

     fich << "0 "<< r[i][0] << ", " << r[i][1] <<  ", "  << v[i][0] << ", " << v[i][1] <<  ", " << a[i][0] << ", " << a[i][1] << "distancia " <<  distancia(0, 1, r)<< endl;
     }
     fich << endl; 



//ITERACIONES PARA SACAR LAS POSICIONES

for(int i=1; i<=ITER_MAX; i++){
  t=i*h;
    //REALIZAMOS ALGORITMO VERLET
    Verlet(r, a, v, w);


//ESCRIBIMOS POSICIONES
for (int j = 0; j < N; j++) {

file << r[j][0] << ", " << r[j][1] << endl;
}
file << endl; 

  for (int j= 0; j < N; j++) {

     fich<< t << " "<<  r[j][0] << ", " << r[j][1] <<  ", "  << v[j][0] << ", " << v[j][1] <<  ", " << a[j][0] << ", " << a[j][1] <<endl;
     }
     fich << endl; 

T_i (v, Ti);


//PROMEDIAMOS VELOCIDAD
if (t>=20.0 && t<=50.0){
TEMP_media=TEMP_media+Ttot(Ti)/N;
}
//PROMEDIAMOS PRESION 
if (t>=20.0 && t<=50.0){
PRESION_media=PRESION_media+presion(r,v,a);
}


//t=t+h e iteramos
}

//Sacamos temperatura media entre t=20 y t=50
TEMP_media=TEMP_media/15001;
PRESION_media=PRESION_media/(4*L*h*15001);
cout <<"T: "<< TEMP_media << " P: "<< PRESION_media;
return 0;

}