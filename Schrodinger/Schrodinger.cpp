#include <iostream>
#include <conio.h>
#include <cmath>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <time.h>
#include "complex.h"




using namespace std;

#define N 100
#define TIEMPO 1000
#define lambda 15

void potencialinicial(double V[N+1], double k0tilde);
void funcioninicial(fcomplex  phi[N+1], double k0tilde, double h);

int main (){
    
int ncic;
double h, k0, k0tilde, stilde, V[N+1], x, norma;
fcomplex phi[N+1], alpha[N], Beta[N], chi[N+1], gamma[N+1], A0[N], uno, coefb;
//Damos parametros iniciales y generamos parametros

h=0.01;

ncic=1;
k0=2.0*M_PI*ncic/(N*h);
k0tilde=k0*h;
stilde=1/(40*k0tilde*k0tilde);
alpha[N-1]= Complex(0.0, 0.0);
Beta[N-1]=Complex(0.0, 0.0);
chi[0]=Complex(0.0, 0.0);
chi[N]=Complex(0.0, 0.0);
uno=Complex(1.0,0.0);
coefb=Complex(0.0, 4.0/stilde);
norma=0.0;

//Inicializamos el potencial y la funcion de onda
potencialinicial(V, k0tilde);
funcioninicial(phi, k0tilde, h);

//Imprimimos estado inicial
ofstream file("schrodinger_data.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "schrodinger_data.dat" << "para leer" << endl;
    }

for (int j = 0; j <= N; j++){
x=j*h;
//file << x << ", "<< Cabs(phi[j])*Cabs(phi[j]) << ", " <<  phi[j].r*phi[j].r << ", " << phi[j].i*phi[j].i << endl;
file << x << ", " << V[j] << ", " << Cabs(phi[j])*Cabs(phi[j]) << ", "   << phi[j].i*phi[j].i << ", "  << phi[j].r*phi[j].r << endl;
norma=norma + Cabs(phi[j])*Cabs(phi[j]);
}

file << endl;
cout<< "La norma es " << norma << endl;

//Calculamos los alpha
for (int j=N-1; j>0; j--){
A0[j]= Complex(-2.0-V[j], 2.0/stilde);
gamma[j]=Cdiv(uno, Cadd(A0[j], alpha[j]));
alpha[j-1]=RCmul(-1.0,  gamma[j]);
}





for (int n=0; n<TIEMPO; n++){
  norma=0;
//Calculamos los Beta
for (int j=N-1; j>0; j--){
Beta[j-1]=Cmul(gamma[j], Csub(Cmul(coefb, phi[j] ), Beta[j]));
}

//Calculamos chi
for (int j=0; j<N; j++){
  chi[j+1]=Cadd(Cmul(alpha[j],chi[j]),Beta[j]);
}

//Calculamos phi
for (int j=1; j<N; j++){
  phi[j]=Csub(chi[j], phi[j]);
}

//Imprimimos estado
for (int j = 0; j <= N; j++){
x=j*h;
//file << x << ", "<< Cabs(phi[j])*Cabs(phi[j]) << ", " <<  phi[j].r*phi[j].r << ", " << phi[j].i*phi[j].i << endl;
file << x << ", " << V[j] << ", " << Cabs(phi[j])*Cabs(phi[j]) << ", "   << phi[j].i*phi[j].i << ", "  << phi[j].r*phi[j].r << endl;
norma=norma + Cabs(phi[j])*Cabs(phi[j]);
}

file << endl;
cout<< "La norma es " << norma << endl;

}

}
 

void potencialinicial(double V[N+1], double k0tilde){
for(int j=0; j<N+1; j++){
    if (j<2*N/5){
    V[j]=0    ;
    }
    else if (j>3*N/5){
      V[j]=0  ;
    }
    else{
    V[j]=lambda*k0tilde*k0tilde;
    }
}
return;
}
void funcioninicial(fcomplex phi[N+1], double k0tilde, double h){
for(int j=0; j<N+1; j++){
    if (j==0){
    phi[j]=Complex(0.0,0.0)  ;
    }
    else if (j==N){
      phi[j]=Complex(0.0,0.0)  ;
    }
    else{
    phi[j]= Cgauss(k0tilde*j, exp((-8.0*(4.0*j-N)*(4.0*j-N))/(1.0*N*N))/3.3283);
    }
}
}