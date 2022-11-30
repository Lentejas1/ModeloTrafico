#include <cmath>
#include <iostream>

using namespace ::std;
#define N 2000 /* N es el número de pasos*/

double m = 1500, L = 4, C = 20000, v = 120 / 3.6;
double d = 26, t_0 = C / (m * d);
double dt = 1. / N;
double X[N][5]; /* Fila tiempo, columna coche*/
double V[N][5]; /* Fila tiempo, columna coche*/


double rk4(){
    for (int i = 0; i < N; i++) { /* Para cada t*/
        for (int j = 1; j < 5; j++){ /* Para cada K y L */
            if (j == 1){
                for (int k = 1; k < 5; k++){
                    K1 = V[i][k];
                    L1 = f();
                }
                /*Que haga lo de K1 y L1*/
            }
            else if (j == 2){
                /*Que haga lo de K2 y L2*/
            }
            else if (j == 3){
                /*Que haga lo de K3 y L3*/
            }
            else{
                /*Que haga lo de  K4 y L4*/
            }

        }
    }
}

int main(){
    double t_n, tau_n;
    double k[4][N], l[4][N]; /*Filas son coches y columans pasos de tiempo*/

return 0;
}