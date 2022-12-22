#include <cmath>
#include <iostream>

using namespace ::std;
#define N 1000 /* N es el número de pasos*/

double m = 1500, l = 4, C = 20000, v = 120 / 3.6;
double d = 26, t_0 = C / (m * d);
double dt = 1. / N;
double t_c=2;


void respuesta_instantanea(){
    double X[N][5]; /* Fila tiempo, columna coche*/
    double V[N][5]; /* Fila tiempo, columna coche*/
    int i_norm = int(t_c / dt);
    for (int i=0; i<N-1; i++){
        double K[5][4], L[5][4]; /*Filas son coches y columnas K sub algo*/
        double v_norm = v * m / C;
        V[0][0] = v_norm;
        X[0][0] = v_norm * i * dt;
        if (i <= i_norm){
            X[i+1][0] = v_norm * (i+1) * dt, V[i+1][0] = v_norm;
        }
        else {
            X[i+1][0] = X[i][0] + 0.2 * v_norm * dt, V[i+1][0] = 0.2 * v_norm;
        }
        X[0][0]=0;
        for (int j = 1; j < 5; j++){
            X[0][j] = X[0][j-1] - d / d - l / d;
            V[0][j] = v_norm;
        }



        K[0][0] = 0, K[0][1] = 0 + dt / 2 * K[0][0]; // Hay que restar un índice (
        K[0][2] = 0 + dt / 2 * K[0][1], K[0][3] = dt * K[0][2];

        L[0][0] = v_norm, L[0][1] = v_norm + dt / 2 * L[0][0];
        L[0][2] = v_norm + dt / 2 * K[0][1], K[0][3] = v_norm + dt * K[0][2];

        for (int j=1; j<5; j++){
            K[j][0] = -(V[i][j] - V[i][j-1])/(fabs(X[i][j] - X[i][j-1]));
            L[j][0] = V[i][j];

            K[j][1] = -((V[i][j] + dt*0.5*K[j][0])- (V[i][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i][j-1] + dt*0.5*L[j-1][0])));
            L[j][1] = V[i][j] + dt*0.5*L[j][0];

            K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i][j-1] + dt*0.5*L[j-1][1])));
            L[j][2] = V[i][j] + dt*0.5*K[j][1];

            K[j][3] = -((V[i][j] + dt*K[j][2])- (V[i][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i][j-1] + dt*L[j-1][2])));
            L[j][3] = V[i][j] + dt*K[j][2];
        }
        for (int k=1; k < 5; k++){
            X[i+1][k] = X[i][k] + dt / 6 * (L[k][0] * 2*L[k][1] + 2*L[k][2]+L[k][3]);
            V[i+1][k] = V[i][k] + dt / 6 * (K[k][0] * 2*K[k][1] + 2*K[k][2]+K[k][3]);
        }
    }

    FILE* data;
    data = fopen("Instantaneo.csv", "w");
    fprintf(data, "t (s), x_0 (m), x_1 (m), x_2 (m), x_3 (m), x_4 (m)\n");

    for (int i = 0; i < N - 1; i++){
            fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i*t_0, X[i][0]*d,X[i][1]*d,X[i][2]*d,X[i][3]*d,X[i][4]*d);
    }
    fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf", dt*N*t_0, X[N-1][0]*d,X[N-1][1]*d,X[N-1][2]*d,X[N-1][3]*d,X[N-1][4]*d);

    FILE* data2;
    data2 = fopen("Velocidades.csv", "w");
    fprintf(data2, "t (s), v_0 (m/s), v_1 (m/s), v_2 (m/s), v_3 (m/s), v_4 (m/s)\n");
    double v_norm = v * m / C;
    for (int i = 0; i < N - 1; i++){
        fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i*t_0, V[i][0]*v_norm,V[i][1]*v_norm,V[i][2]*v_norm,V[i][3]*v_norm,V[i][4]*v_norm);
    }
    fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf", dt*N*t_0, V[N-1][0]*v_norm,V[N-1][1]*v_norm,V[N-1][2]*v_norm,V[N-1][3]*v_norm,V[N-1][4]*v_norm);

}
int main(){
    cout << N * t_0;

    respuesta_instantanea();
return 0;
}