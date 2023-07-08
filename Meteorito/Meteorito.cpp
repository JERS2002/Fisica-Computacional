#include <iostream>
#include <conio.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <fstream>
#include <cmath>


using namespace std;

#define PASOS 20000
#define G 6.67e-11
#define MT 5.9736e24
#define ML 0.07349e24
#define Mmet 8.7127e15
#define dTL 3.844e8
#define omega 2.6627e-6
#define RT 6.37816e6
#define RL 1.7374e6
#define m 20000
#define radiometeorito 10000
#define T 13.49e6 //ENERGIA DE LA DETONACION EN KILOTONES

double navef1(double pr);
double navef2(double r, double pphi);
double navef3 (double r, double phi, double pr, double pphi, double delta, double ML_red, double h, double i);
double navef4 (double r, double phi, double pr, double pphi, double delta, double ML_red, double h, double i);
void rgnave(double h, double& r, double& phi, double& pr, double& pphi, double ML_red , double Mmet_red , double delta, double phi_L0, double rmet,  int i );
void rgmet(double h, double& prmet, double& rmet, double delta);
void ppolar_cart (double& px, double& py, double r, double phi, double pr, double pphi, double& v);
void pcart_polar (double px, double py, double r, double phi, double& pr, double& pphi);
void impulso(double modulo, double angulo, double px, double py, double r, double phi, double& pr, double& pphi, bool dirmov, double v);
double trozometf1 (double prtrozo);
double metf2 (double rmet, double delta);
double trozometf3 (double  rtrozo, double phitrozo, double prtrozo, double pphitrozo, double delta, double ML_red, double h, double phi_L0, int i);
double trozometf4 (double  rtrozo, double phitrozo, double prtrozo, double pphitrozo, double delta, double ML_red, double h, double phi_L0, int i);
void rgtrozos(double h, double& r, double& phi, double& pr, double& pphi, double ML_red, double delta, double phi_L0,  int i );

int main ()
{

double r, phi, pr, pphi, px, py, v0, theta0, delta, h, ML_red, Mmet_red, rmet, prmet, phimet, pphimet, phi_L0, Hprime, nave_met, E, deltav, v, rmin, rmetmin;
nave_met=10;

//bool dirmov, angulodado; //Variables utilizadas para elegir si el impulso se da en la direccion de movimiento o en una direccion dada
//dirmov=true;
//angulodado=false;

//VARIABLES AUXILIARES

delta= G*MT/(dTL*dTL*dTL);
ML_red=ML/MT;
Mmet_red=Mmet/MT;
h= 60.0; //PASO TEMPORAL

//CONDICIONES INICIALES EN EL LANZAMIENTO
r= RT/dTL;
phi=90* M_PI/180;
phi_L0=6 * M_PI/180;

v0=11200;
theta0= 55.75*M_PI/180; //Theta0 da la direccion de la velocidad inicial, ver diapositivas cohete
//Cálculo energía
E=0.5*v0*v0;
deltav=431.2968;

pr=v0*cos(theta0-phi)/(dTL);
pphi=r*v0*sin(theta0-phi)/(dTL);
//CONDICIONES INICIALES METEORITO
rmet=10.0;
prmet= -5000/dTL;

/*
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
*/
// Escribimos la posicion inicial en el archivo donde guardaremos las posiciones
    ofstream file("meteorito.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "meteorito.dat" << "para leer" << endl;
    }

file << r*cos(phi) << ", " <<r*sin(phi)  << endl; //POSICION NAVE
file << 0.0 << ", " << rmet << endl;   //POSICION METEORITO
file << cos(phi_L0)<< ", " << sin(phi_L0) << endl ; //POSICION LUNA
file << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

file << endl; 

   ofstream fichenerg("energia.dat");
    if (!file) {
        cerr << "No se pudo abrir " << "meteorito.dat" << "para leer" << endl;
    }
    ofstream fichtrozos("trozos.dat");
    if (!fichtrozos) {
        cerr << "No se pudo abrir " << "trozos.dat" << "para leer" << endl;
    }



// Fichero auxiliar con datos que usaremos para ajustar las trayectorias
    ofstream fich("datos.dat");
    if (!fich) {
        cerr << "No se pudo abrir " << "datos.dat" << "para leer" << endl;
    }
fich << "TIEMPO (s): " << 0 << endl; 
fich << r*cos(phi) << ", " <<r*sin(phi) << endl; //POSICION NAVE
fich << 0.0 << ", " << rmet << endl;   //POSICION METEORITO
fich << cos(phi_L0)<< ", " << sin(phi_L0) << endl ; //POSICION LUNA
fich << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

fich << endl; 

fichenerg<< "0" << " " << "0 "  <<endl;

int i=1;
//APROXIMACION AL METEORITO PASANDO POR LA LUNA
 while (nave_met> radiometeorito/dTL && i<PASOS)
 {
   // EVOLUCION NAVE Y METEORITO
    rgnave(h, r, phi, pr, pphi, ML_red, Mmet_red, delta, phi_L0, rmet,  i);
    rgmet(h, prmet, rmet, delta);
    nave_met=sqrt(rmet*rmet+r*r-2*r*rmet*sin(phi));
   
//IMPULSOS
if (i>=2699 && i<= 2700){
    cout << i*h << " ";
impulso(-1.122e-6, 0, px, py, r, phi, pr, pphi, true, v);
 cout << endl; 
 ppolar_cart(px, py, r, phi, pr, pphi, v);
    E=E+0.5*deltav*(2*v+deltav);
}

if (i>=10389 && i<= 10401){
    cout << i*h << " ";
impulso(1.122e-6, -90, px, py, r, phi, pr, pphi, false, v);
cout << endl; 
ppolar_cart(px, py, r, phi, pr, pphi, v);
    E=E+0.5*deltav*(2*v+deltav);
}

if (i==10462){
    cout << i*h << " ";
    impulso(1.122e-6, 0, px, py, r, phi, pr, pphi, false, v);
    cout << endl; 
    ppolar_cart(px, py, r, phi, pr, pphi, v);
    E=E+0.5*deltav*(2*v+deltav);
}
if (i==10472){
    cout << i*h << " ";
    impulso(1.122e-6, -26.26, px, py, r, phi, pr, pphi, false, v);
    cout << endl; 
    ppolar_cart(px, py, r, phi, pr, pphi, v);
    E=E+0.5*deltav*(2*v+deltav);
}




 //SACAMOS DATOS A FICHERO
    if((i%50)==0){
     file << r*cos(phi) << ", " << r*sin(phi)  << endl; //POSICION NAVE
     file << 0.0 << ", " << rmet << endl;   //POSICION METEORITO
     file << cos(phi_L0 + omega*h*i)<< ", " << sin(phi_L0 + omega*h*i) << endl ; //POSICION LUNA
     file << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

     file << endl; 
    }
ppolar_cart(px,py,r, phi, pr, pphi,v);
if(i>=10000 && i<=11000){
//FICHERO DE DATOS AUXILIAR
     fich << "TIEMPO (s): " << i*h << endl; 
     fich << "iteracion " << i << endl ;
     fich << r*cos(phi) << ", " <<r*sin(phi) << " " << px<< " "<< py  << endl; //POSICION NAVE
     fich << 0.0 << ", " << rmet << " " << prmet << endl;   //POSICION METEORITO
     fich << cos(phi_L0 + omega*h*i)<< ", " << sin(phi_L0 + omega*h*i) << endl ; //POSICION LUNA
     fich << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

     fich << endl; 
     }   

//Fichero de datos energia
fichenerg<< i*h/86400 << " " << E/(pow(10, 6))  <<endl; //Fichero contiempo en días y energía en MJ

     i++;
 }



//CUANDO nave_met sea menor o igual que el radio del meteorito, ocurre la DETONACION
//Cuando ocurra la detonacion, igualamos las coordenadas de la nave a las del meteorito, ya que para evitar introducir mas variables, utilizaremos las anteriores para cada trozo del meteorito
cout<< "TIEMPO aterrizaje (s):  " << i*h << endl;

r=rmet;
pr=prmet/2;
prmet=pr;
phi=90* M_PI/180;//Trozo 1
phimet=90* M_PI/180; //Trozo 2
pphi=sqrt(T*4.184e12/(2*Mmet))*rmet/(dTL); //Momento tras la explosion para trozo 1
pphimet=-pphi; // Idem trozo 2

rmetmin=100;
rmin=100;


//Una vez impuestas las condiciones iniciales tras la detonacion, estudiamos las trayectorias de los trozos
int j;
j=i;
for (i=j+1; i<=j+PASOS; i++){

    rgtrozos(h, r, phi, pr, pphi, ML_red, delta, phi_L0, i);
    rgtrozos (h, rmet, phimet, prmet, pphimet, ML_red, delta, phi_L0, i);
if (rmet<rmetmin){
    rmetmin=rmet;
}
if (r<rmin){
    rmin=r;
}

 //SACAMOS DATOS A FICHERO 
    if((i%50)==0){
     file << r*cos(phi) << ", " << r*sin(phi)  << endl; //POSICION NAVE
     file << rmet*cos(phimet) << ", " << rmet*sin(phimet) << endl;   //POSICION METEORITO
     file << cos(phi_L0 + omega*h*i)<< ", " << sin(phi_L0 + omega*h*i) << endl ; //POSICION LUNA
     file << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

     file << endl; 

//FICHERO DE DATOS AUXILIAR
     fich << "TIEMPO (s): " << i*h << endl; 
     fich << "iteracion " << i << endl ;
     fich << r*cos(phi) << ", " <<r*sin(phi) << " " << pr << " "<< pphi  << endl; //POSICION NAVE
     fich << rmet*cos(phimet)<< ", " <<  rmet*sin(phimet) << endl;   //POSICION METEORITO
     fich << cos(phi_L0 + omega*h*i)<< ", " << sin(phi_L0 + omega*h*i) << endl ; //POSICION LUNA
     fich << 0.0  << ", " << 0.0 << endl ; // POSICION TIERRA

     fich << endl; 
     }   
 //Fichero graficas posiciones trozos
fichtrozos<< (i-j)*h/86400 << " " << rmet << " " << r << endl;
}

cout<< "Distancia minima trozo 1:" <<rmetmin;
cout<< "Distancia minima trozo 2:" <<rmin;
}







// FUNCIONES EMPLEADAS (ecuaciones del movimiento )


double navef1(double pr){
    return pr;
}

double navef2(double r, double pphi){
    return pphi/(r*r);
}
double navef3 (double r, double phi, double pr, double pphi, double delta, double ML_red, double h, double phi_L0, double Mmet_red, double rmet, int i){
    return pphi*pphi/(r*r*r) - delta*(1/(r*r)+ML_red*(r-cos(phi- phi_L0 -omega*h*i))/(pow (1+r*r-2*r*cos(phi-phi_L0-omega*h*i), 1.5)) +Mmet_red*(r-rmet*sin(phi))/(pow (rmet*rmet+r*r-2*r*rmet*sin(phi), 1.5)));
}
double navef4 (double r, double phi, double pr, double pphi, double delta, double ML_red, double h, double phi_L0, double Mmet_red,  double rmet,  int i){
  return delta*(-ML_red*r*sin(phi-phi_L0-omega*h*i)/(pow (1+r*r-2*r*cos(phi- phi_L0 -omega*h*i), 1.5)) +Mmet_red*r*cos(phi)/(pow (rmet*rmet+r*r-2*r*rmet*sin(phi), 1.5)));
}

double metf1 (double prmet){
    return prmet;
}
double metf2 (double rmet, double delta){
return -delta/(rmet*rmet);
}

double trozometf1 (double prtrozo){
    return prtrozo;
}
double trozometf2 (double rtrozo, double pphitrozo){
    return pphitrozo/(rtrozo*rtrozo);
}
double trozometf3 (double  rtrozo, double phitrozo, double prtrozo, double pphitrozo, double delta, double ML_red, double h, double phi_L0, int i){
return pphitrozo*pphitrozo/(rtrozo*rtrozo*rtrozo) - delta*(1/(rtrozo*rtrozo)+ML_red*(rtrozo-cos(phitrozo- phi_L0 -omega*h*i))/(pow (1+rtrozo*rtrozo-2*rtrozo*cos(phitrozo-phi_L0-omega*h*i), 1.5)));
}
double trozometf4 (double  rtrozo, double phitrozo, double prtrozo, double pphitrozo, double delta, double ML_red, double h, double phi_L0, int i) {
 return delta*(-ML_red*rtrozo*sin(phitrozo-phi_L0-omega*h*i)/(pow (1+rtrozo*rtrozo-2*rtrozo*cos(phitrozo- phi_L0 -omega*h*i), 1.5)));
}

//RUNGE_KUTTA PARA LA NAVE
void rgnave(double h, double& r, double& phi, double& pr, double& pphi, double ML_red , double Mmet_red , double delta, double phi_L0, double rmet,  int i ){
    double k1[4], k2[4], k3[4], k4[4];
 //Evaluamos k(1)_n
    k1[0]=h*navef1(pr);
    k1[1]=h*navef2(r, pphi);
    k1[2]=h*navef3(r, phi, pr, pphi, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i);
    k1[3]=h*navef4(r, phi, pr, pphi, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i);

    //Evaluamos k(2)_n
   k2[0]=h*navef1(pr+k1[2]/2);
    k2[1]=h*navef2(r+k1[0]/2, pphi+k1[3]/2);
    k2[2]=h*navef3(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2, pphi+k1[3]/2, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1/2);
    k2[3]=h*navef4(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2, pphi+k1[3]/2, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1/2);

   // Evaluamos k(3)_n
k3[0]=h*navef1(pr+k2[2]/2);
    k3[1]=h*navef2(r+k2[0]/2, pphi+k2[3]/2);
    k3[2]=h*navef3(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1/2);
    k3[3]=h*navef4(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1/2);


   // Evaluamos k(4)_n
k4[0]=h*navef1(pr+k1[2]/2);
    k4[1]=h*navef2(r+k3[0], pphi+k3[3]);
    k4[2]=h*navef3(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1);
    k4[3]=h*navef4(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, ML_red, h, phi_L0, Mmet_red, rmet, 1.0*i+1);

 
        r=r+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
        phi=phi+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        pr=pr+(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        pphi=pphi+(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;

 }

 //RUNGE KUTTA METEORITO

 void rgmet(double h, double& prmet, double& rmet, double delta){
 double k1[2], k2[2], k3[2], k4[2];
  //Evaluamos k(1)_n
    k1[0]=h*metf1(prmet);
    k1[1]=h*metf2(rmet, delta);
 //Evaluamos k(2)_n
   k2[0]=h*metf1(prmet+k1[1]/2);
    k2[1]=h*metf2(rmet+k1[0]/2, delta);
    //Evaluamos k(3)_n
   k3[0]=h*metf1(prmet+k2[1]/2);
    k3[1]=h*metf2(rmet+k2[0]/2, delta);
    //Evaluamos k(4)_n
    k4[0]=h*metf1(prmet+k3[1]);
    k4[1]=h*metf2(rmet+k3[0], delta);

    rmet=rmet+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
    prmet=prmet+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
 }

// RUNGE KUTTA PARA LOS TROZOS DE METEORITO
void rgtrozos(double h, double& r, double& phi, double& pr, double& pphi, double ML_red, double delta, double phi_L0,  int i ){
    double k1[4], k2[4], k3[4], k4[4];
 //Evaluamos k(1)_n
    k1[0]=h*trozometf1(pr);
    k1[1]=h*trozometf2(r, pphi);
    k1[2]=h*trozometf3(r, phi, pr, pphi, delta, ML_red, h, phi_L0, 1.0*i);
    k1[3]=h*trozometf4(r, phi, pr, pphi, delta, ML_red, h, phi_L0, 1.0*i);

    //Evaluamos k(2)_n
   k2[0]=h*trozometf1(pr+k1[2]/2);
    k2[1]=h*trozometf2(r+k1[0]/2, pphi+k1[3]/2);
    k2[2]=h*trozometf3(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2, pphi+k1[3]/2, delta, ML_red, h, phi_L0,  1.0*i+1/2);
    k2[3]=h*trozometf4(r+k1[0]/2, phi+k1[1]/2, pr+k1[2]/2, pphi+k1[3]/2, delta, ML_red, h, phi_L0, 1.0*i+1/2);

   // Evaluamos k(3)_n
k3[0]=h*trozometf1(pr+k2[2]/2);
    k3[1]=h*trozometf2(r+k2[0]/2, pphi+k2[3]/2);
    k3[2]=h*trozometf3(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, ML_red, h, phi_L0, 1.0*i+1/2);
    k3[3]=h*trozometf4(r+k2[0]/2, phi+k2[1]/2, pr+k2[2]/2,pphi+k2[3]/2, delta, ML_red, h, phi_L0,  1.0*i+1/2);


   // Evaluamos k(4)_n
k4[0]=h*trozometf1(pr+k1[2]/2);
    k4[1]=h*trozometf2(r+k3[0], pphi+k3[3]);
    k4[2]=h*trozometf3(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, ML_red, h, phi_L0, 1.0*i+1);
    k4[3]=h*trozometf4(r+k3[0], phi+k3[1], pr+k3[2],pphi+k3[3], delta, ML_red, h, phi_L0, 1.0*i+1);

 
        r=r+(k1[0]+2*k2[0]+2*k3[0]+k4[0])/6;
        phi=phi+(k1[1]+2*k2[1]+2*k3[1]+k4[1])/6;
        pr=pr+(k1[2]+2*k2[2]+2*k3[2]+k4[2])/6;
        pphi=pphi+(k1[3]+2*k2[3]+2*k3[3]+k4[3])/6;

 }
//Transformaciones de momento de cartesianas a polares y viceversa
void ppolar_cart (double& px, double& py, double r, double phi, double pr, double pphi, double& v){
px= pr*cos(phi)-pphi*sin(phi)/r;
py= pr*sin(phi)+pphi*cos(phi)/r;
v=sqrt(px*px+py*py)*dTL;
}
void pcart_polar (double px, double py, double r, double phi, double& pr, double& pphi){
pr=px*cos(phi)+py*sin(phi);
pphi=(py*cos(phi)-px*sin(phi))*r;
}

//FUNCION PARA DAR IMPULSO

void impulso(double modulo, double angulo, double px, double py, double r, double phi, double& pr, double& pphi, bool dirmov, double v){
     ppolar_cart(px, py, r, phi, pr, pphi, v);
    
     if (dirmov==true){  //SI DIRMOV==TRUE SE OBVIARA EL ANGULO DADO Y EL IMPULSO SE REALIZARA EN LA DIRECCION DEL MOVIMIENTO
        if (px>0){
            if (py>0){
            angulo =180*atan(abs(py/px))/M_PI;
            }
            else{
                angulo= -180*atan(abs(py/px))/M_PI;
            }
        }               
        else {
            if (py>0){
            angulo =180-180*atan(abs(py/px))/M_PI;
            }
            else{
                angulo= -180 +180*atan(abs(py/px))/M_PI;
            }
        }
     }
    cout<< angulo;
    px=px+modulo*cos(angulo*M_PI/180);
    py=py+modulo*sin(angulo *M_PI/180); 
    pcart_polar(px, py, r, phi, pr, pphi);
}