#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>

using namespace std;

#define PASOS 5000
#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define dTL 3.844e8
#define omega 2.6627e-6
#define RT 6.37816e6
#define RL 1.7374e6
#define m 20000

double f1(double pr);
double f2(double r, double pphi);
double f3 (double r, double phi, double pr, double pphi, double delta, double mu, double h, double i);
double f4 (double r, double phi, double pr, double pphi, double delta, double mu, double h, double i);


int main ()
{

double r, phi, pr, pphi, v0, theta0, delta, h, mu, k1[4], k2[4], k3[4], k4[4], Hprime;
r= RT/dTL;
phi=90* M_PI/180;

v0=11200;
theta0= 50*M_PI/180;



pr=v0*cos(theta0-phi)/(dTL);
pphi=-r*v0*sin(theta0-phi)/(dTL);




delta= G*MT/(dTL*dTL*dTL);

mu=0.07349/5.9736;
h= 60.0;

Hprime=pr*pr*m*dTL*dTL/2+ pphi*pphi*m*dTL*dTL/(2*r*r) - G*m*MT/(r*dTL)- G*m*ML/(dTL*sqrt(r*r+1-2*r*cos(phi)))-omega*m*dTL*dTL*pphi;
//Escribimos H' en fichero
ofstream fich("Hprime.dat");
    if (!fich) {
        cerr << "No se pudo abrir " << "cohete.dat" << "para leer" << endl;
    }

fich <<"0  "<< Hprime << endl;
// Escribimos la posicion inicial en el archivo donde guardaremos las posiciones
    ofstream file("cohete.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "cohete.dat" << "para leer" << endl;
    }



file << r*cos(phi) << ", " <<r*sin(phi)  << endl;
file << 1.0 << ", " << 0.0 << endl ;
file << 0.0  << ", " << 0.0 << endl ;

file << endl; 

 for (int i=1; i<PASOS; i++)
 {
    //Evaluamos k(1)_n
    k1[0]=h*f1(pr);
    k1[1]=h*f2(r, pphi);
    k1[2]=h*f3(r, phi, pr, pphi, delta, mu, h, 1.0*i);
    k1[3]=h*f4(r, phi, pr, pphi, delta, mu, h, 1.0*i);

    //Evaluamos k(2)_n
   k2[0]=h*f1(pr+k1[2]/2);
    k2[1]=h*f2(r+k1[0]/2, pphi+k1[3]/2);
    k2[2]=h*f3(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2,pphi+k1[3]/2, delta, mu, h, 1.0*i+1/2);
    k2[3]=h*f4(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2,pphi+k1[3]/2, delta, mu, h, 1.0*i+1/2);

   // Evaluamos k(3)_n
k3[0]=h*f1(pr+k2[2]/2);
    k3[1]=h*f2(r+k2[0]/2, pphi+k2[3]/2);
    k3[2]=h*f3(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, mu, h, 1.0*i+1/2);
    k3[3]=h*f4(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, mu, h, 1.0*i+1/2);


   // Evaluamos k(4)_n
k4[0]=h*f1(pr+k1[2]/2);
    k4[1]=h*f2(r+k3[0], pphi+k3[3]);
    k4[2]=h*f3(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, mu, h, 1.0*i+1);
    k4[3]=h*f4(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, mu, h, 1.0*i+1);

    //Evaluamos y_n (t+h)
    for (int n=0; n<4; n++){
        r=r+(k1[0]+k2[0]+k3[0]+k4[0])/6;
        phi=phi+(k1[1]+k2[1]+k3[1]+k4[1])/6;
        pr=pr+(k1[2]+k2[2]+k3[2]+k4[2])/6;
        pphi=pphi+(k1[3]+k2[3]+k3[3]+k4[3])/6;
    }

// Imprimimos Hprime
Hprime=pr*pr*m*dTL*dTL/2+ pphi*pphi*m*dTL*dTL/(2*r*r) - G*m*MT/(r*dTL)- G*m*ML/(dTL*sqrt(r*r+1-2*r*cos(phi-omega*i*h)))-omega*m*dTL*dTL*pphi;
fich << i << "  " << Hprime <<endl;
//Imprimimos estado
if((i%5)==0){
file << r*cos(phi) << ", " <<r*sin(phi)  << endl;
file << cos(omega*i*h) << ", " << sin(omega*i*h) << endl ;
file << 0.0  << ", " << 0.0 << endl ;

file << endl; 
 }
  }




}

double f1(double pr){
    return pr;
}

double f2(double r, double pphi){
    return pphi/(r*r);
}
double f3 (double r, double phi, double pr, double pphi, double delta, double mu, double h, double i){
    return pphi*pphi/(r*r*r) - delta*(1/(r*r)+mu*(r-cos(phi-omega*h*i))/(pow (1+r*r-2*r*cos(phi-omega*h*i), 3/2)));
}
double f4 (double r, double phi, double pr, double pphi, double delta, double mu, double h, double i){
  return -delta*mu*r*sin(phi-omega*h*i)/(pow (1+r*r-2*r*cos(phi-omega*h*i), 3/2));
}