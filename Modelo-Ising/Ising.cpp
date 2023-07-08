#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <time.h>

using namespace std;
#define N 100
#define Pasos 800
#define T 3
//Funcion para inicializar configuracion ordenada (Todo 1)
void inicord(int s[N][N])
{

for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
s[i][j]=1;
    }
}
    return;
}



//Funcion para inicializar configuracion desordenada (aleatoria)
void inicdesord(int s[N][N])
{
for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
        s[i][j]=-1 + 2*(rand()%2);

    }
}    
return;    
}

double energia(int s[N][N], int n, int m)
{
double E;
int nmas, nmenos, mmas, mmenos;
nmas=n+1;
nmenos=n-1;
mmas=m+1;
mmenos=m-1;

if (nmas==N){
    nmas=0;
}
if (nmenos==-1){
    nmas=N-1;
}
if (mmas==N){
    mmas=0;
}
if (mmenos==-1){
    nmas=N-1;
}

E=2.0*s[n][m]*(s[nmas][m]+s[nmenos][m]+s[n][mmas]+s[n][mmenos]);

return E;
}


int main (){
double E, p, chi;
int s[N][N], n, m;
srand(time(NULL));


//Iniciamos el sistema ordenado (Todo 1)
// inicord(s);

//Iniciamos el sistema desordenado 
 inicdesord(s);

//Imprimimos estado inicial
ofstream file("ising_data.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "ising_data.dat" << "para leer" << endl;
    }
for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++){
if (j==N-1){
file << s[i][j] << endl;
}
else{
file << s[i][j] << ", ";
}

}
}
file << endl;




for (int k=1; k<=Pasos; k++)
{
for (int t=1; t<=N*N; t++){
m=rand()%N;
n=rand()%N;


E=energia(s, n, m);


if (exp(-E*1.0/T)<1){
    p=exp(-E*1.0/T);
}
else{
    p=1;
}


chi=rand()*1.0/RAND_MAX;
if (chi<p){
    s[n][m]=-s[n][m];
}




}
//Imprimimos estado 

for (int i = 0; i < N; i++) {
for (int j = 0; j < N; j++){
if (j==N-1){
file << s[i][j] << endl;
}
else{
file << s[i][j] << ", ";
}

}
}
file << endl;

}



return 0;

}