#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

using namespace std;

#define G 6.67e-11
#define Ms 1.99e30
#define c 1.496e11




//Definimos funciones para calcular aceleracion en tiempo t con ley de Gravitacion de Newton
//Aceleracion eje x
double modulo (int i, int j, double r[9][2])
{
double modulo;
modulo= sqrt(pow(r[i][0]-r[j][0],2) + pow(r[i][1]-r[j][1],2));
return modulo;
}


void aceleracion_x (double r[9][2], double a[9][2], double m[9])
{
 double suma;
 suma=0.0;
 for (int i=0; i<9; i++){
 for (int j=0; j<9; j++) {
 if (j!=i){
   suma= suma +(m[j]*(r[i][0]-r[j][0]))/(pow(modulo(i,j,r), 3));
 }
 }
 a[i][0]=-suma; 
 suma = 0.0;
 }
return;
 }

//Aceleracion eje y
void aceleracion_y (double r[9][2], double a[9][2], double m[9])
{ 
double suma;
 suma=0.0;
 for (int i=0; i<9; i++){
 for (int j=0; j<9; j++) {
 if (j!=i){
   suma= suma +(m[j]*(r[i][1]-r[j][1]))/(pow(modulo(i,j,r), 3));
  
 }
 }
 a[i][1]=-suma; 
 suma = 0.0;
 }
return;
}
//Funciones para calcular la energía potencial Vi de cada planeta y la energia total V
void V_i(double m[9], double r[9][2], double Vi[9])
{
double suma;
suma=0.0;
for (int i=0; i<9; i++){
 for (int j=0; j<9; j++) {
 if (j!=i){
   suma= suma +m[j]/(modulo(i,j,r)); 
 }
 }
 Vi[i]= -G *Ms*Ms*m[i]*suma/c; 
 suma = 0.0;
 }
return;
}
double Vtot(double Vi[9]){
   double V;
   V=0.0;
    for (int i=0; i<9; i++){
    V= V+Vi[i];
    }
    return V;
}

//Funcion para calcular la energia cinetica Ti de cada planeta y la energia total T
void T_i(double m[9], double v[9][2], double Ti[9]){
for (int i=0; i<9; i++){
    Ti[i]=Ms*Ms*G*m[i]*(v[i][0]*v[i][0]+v[i][1]*v[i][1])/(2*c);
}
return;
}

double Ttot(double Ti[9]){
   double T;
   T=0.0;
    for (int i=0; i<9; i++){
    T= T+Ti[i];
    }
    return T;
}

//Funcion calcular r(h) por Desarrollo en Serie de Taylor
void r_h(double r[9][2], double v[9][2], double a[9][2], double h){
for (int j=0; j<2; j++ ){
for (int i=0; i<9; i++){
    r[i][j]=r[i][j]+h*v[i][j]+h*h*a[i][j]/2;
    }
}
return;
}

//Funcion calcular wi= vi+ h/2 *a_i
void w_i(double v[9][2], double a[9][2], double w[9][2], double h){
 for (int j=0; j<2; j++ ){
for (int i=0; i<9; i++){
    w[i][j]=v[i][j]+h*a[i][j]/2;
    }
}
return;   
}

// Funcion calcular v_i(h)= w_i(h)+ h/2 *a_i(h)
void v_i(double w[9][2], double a[9][2], double v[9][2], double h){
    for (int j=0; j<2; j++ ){
for (int i=0; i<9; i++){
    v[i][j]=w[i][j]+ h*a[i][j]/2;
    }
}
return;  
}



//Leer matriz
void leer_matriz(double matriz[9][2], const string& nombrearch) {
    ifstream infile(nombrearch);
    if (!infile) {
        cerr << "No se pudo abrir '" << nombrearch << "para leer" << endl;
        return;
    }
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 2; j++) {
            infile >> matriz[i][j];
        }
    }
    infile.close();
}



 //Leer Masas
void leer_masas(double *arr, const char *nombrearch) {
    ifstream file(nombrearch);
    if (file.is_open()) {
        for (int i = 0; i < 9; i++) {
            file >> arr[i];
        }
        file.close();
    } else {
        cerr << "No se pudo abrir " << nombrearch << std::endl;
    }
}

//Inicializamos los valores de la aceleracion a 0
void inicializar_aceleracion (double a[9][2]){
    for (int i = 0; i < 9; i++) {
        for (int j = 0; j < 2; j++) {
a[i][j]=0.0;
}
}
}

void inicializar_booleano(bool periodo[9])
{
for (int i = 0; i < 9; i++) {
periodo[i]=true;
}
}
 //Main
int main ()
{
    double r[9][2], a[9][2], m[9], v[9][2], w[9][2], Vi[9], Ti[9], h;
    int iteracion, iteracion_max;
    bool periodo[9];
    iteracion=0;
    iteracion_max=10000;
    h=0.1;

    leer_masas (m, "masasestabilidad.txt");
    leer_matriz (r, "Posiciones.txt");
    leer_matriz (v, "Velocidades.txt");
    inicializar_aceleracion (a);
    inicializar_booleano(periodo);
    
// Escribimos la posicion inicial en el archivo donde guardaremos las posiciones
    ofstream file("estabilidad_data.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "planets_data.dat" << "para leer" << endl;
    }
    for (int i = 0; i < 9; i++) {

     file << r[i][0] << ", " << r[i][1] << endl;
     }
     file << endl; 

    //Abrimos archivo para escribir periodos
 ofstream fich("periodos_estabilidaddata.dat");
    if (!fich) {
        cerr << "No se pudo abrir " << "periodos_data.dat" << "para leer" << endl;
    }

//Abrimos fichero para escribir energia
ofstream fichenerg("energiaestabilidad.dat");
    if (!fichenerg) {
        cerr << "No se pudo abrir " << "periodos_data.dat" << "para leer" << endl;
    }
//Escribimos energia en t=0
V_i(m, r, Vi);
T_i (m, r, Ti);
fichenerg <<"0 " << Ttot(Ti)+Vtot(Vi) << endl;

//Evaluamos aceleracion en t=0
aceleracion_x (r,a,m);
aceleracion_y (r,a,m);

for(int i=0; i<iteracion_max; i++){
    iteracion= iteracion + 1;

//Calculamos r_i(t+h) y w_i (t+h)
r_h (r,v,a,h);
w_i (v,a,w,h);



for (int k=1; k<9; k++){
    if ((r[k][1]<0) && periodo[k]) 
    {
        fich<< iteracion*h*2.0*58.1 << endl;
        fich << endl;
        periodo[k]=false;
    }   
}


if((i%10)==0){
for (int j = 0; j < 9; j++) {

file << r[j][0] << ", " << r[j][1] << endl;
}
file << endl; 
}
// Calculamos a_i(t+h)
aceleracion_x (r,a,m);
aceleracion_y (r,a,m);

// Calculamos v_i (t+h)
v_i (w,a,v,h);

//Escribimos energia
V_i(m, r, Vi);
T_i (m, r, Ti);
fichenerg <<iteracion*h* 58.1 << " " << Ttot(Ti)+Vtot(Vi) << endl;


//t=t+h e iteramos
}
return 0;
}