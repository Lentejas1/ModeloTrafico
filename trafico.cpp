#include <cmath>
#include <iostream>
#include <bits/stdc++.h>

using namespace ::std;
#define N 5000 /* N es el número de pasos*/

double m = 1500, l = 4, C = 20000, v = 120 / 3.6;
double d = 26, t_0 = C / (m * d);
double dt = 1. / 500;
double t_c= 2 * t_0;
int already_shown = 0;

void respuesta_instantanea(int reaccion){
    double X[N][5];       /* Fila tiempo, columna coche*/
    double V[N][5];       /* Fila tiempo, columna coche*/
    double K[5][4], L[5][4];          /*Filas son coches y columnas K sub algo*/
    int i_col = int(t_c / dt);
    double v_norm = v * m / C;

    V[0][0] = v_norm;            /*Condiciones iniciales del coche 0*/
    X[0][0] = 0;

    K[0][0] = 0;
    K[0][1] = 0 + dt / 2 * K[0][0];    /*K's del coche 0*/
    K[0][2] = 0 + dt / 2 * K[0][1];
    K[0][3] = dt * K[0][2];



    for (int j = 1; j < 5; j++){
        X[0][j] = X[0][j-1] - (l+d) / d;             /*Posiciones iniciales de los otros 4*/
        V[0][j] = v_norm;                            /*Velocidades iniciales de los otros 4*/
    }

    for (int i=0; i<N-1; i++){
        //Antes del canvio de velocidad
        if (i < i_col){
            X[i+1][0] = X[i][0] + v_norm * dt;
            V[i+1][0] = v_norm;

            L[0][0] = v_norm;
            L[0][1] = v_norm + dt / 2 * K[0][0];
            L[0][2] = v_norm + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_norm + dt * K[0][2];


            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i][j - 1]) / (fabs(X[i][j] - X[i][j - 1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i][j - 1] + dt * 0.5 * K[j - 1][0])) / (fabs((X[i][j] + dt * 0.5 * L[j][0]) - (X[i][j - 1] + dt * 0.5 * L[j - 1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i][j - 1] + dt * 0.5 * K[j - 1][1])) / (fabs((X[i][j] + dt * 0.5 * L[j][1]) - (X[i][j - 1] + dt * 0.5 * L[j - 1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i][j - 1] + dt * K[j - 1][2])) / (fabs((X[i][j] + dt * L[j][2]) - (X[i][j - 1] + dt * L[j - 1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
            //Despues del cambio de velocidad
        else {

            double v_new = 0.2*v_norm;
            X[i+1][0] = X[i][0] + v_new * dt;
            V[i+1][0] = v_new;

            L[0][0] = v_new;
            L[0][1] = v_new + dt / 2 * K[0][0];
            L[0][2] = v_new + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_new + dt * K[0][2];


            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i-reaccion][j-1])/(fabs(X[i][j] - X[i-reaccion][j-1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i-reaccion][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i-reaccion][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i-reaccion][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i-reaccion][j-1] + dt*L[j-1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
    }

    FILE* data;
    data = fopen("Instantaneo1.csv", "w");
    fprintf(data, "t (s), x_0 (m), x_1 (m), x_2 (m), x_3 (m), x_4 (m)\n");

    for (int i = 0; i < N - 1; i++){
        fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, X[i][0]*d,X[i][1]*d,X[i][2]*d,X[i][3]*d,X[i][4]*d);
    }
    fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, X[N-1][0]*d,X[N-1][1]*d,X[N-1][2]*d,X[N-1][3]*d,X[N-1][4]*d);

    FILE* data2;
    data2 = fopen("Velocidades1.csv", "w");
    fprintf(data2, "t (s), v_0 (m/s), v_1 (m/s), v_2 (m/s), v_3 (m/s), v_4 (m/s)\n");
    for (int i = 0; i < N - 1; i++){
        fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, V[i][0]* C/m,V[i][1]* C/m,V[i][2]* C/m,V[i][3]* C/m,V[i][4]* C/m);

        if (V[i][1]*C/m < 33.2){
            if (already_shown == 0){
                cout << "El coche 1 empieza a frenar en t="<<i*dt / t_0 << "s\n";
                already_shown = 1;
            }
        }
    }
    fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, V[N-1][0]* C/m,V[N-1][1]* C/m,V[N-1][2]* C/m,V[N-1][3]* C/m,V[N-1][4]* C/m);

}

void movimiento_exponencial(int reaccion){
    double X[N][5];       /* Fila tiempo, columna coche*/
    double V[N][5];       /* Fila tiempo, columna coche*/
    double K[5][4], L[5][4];          /*Filas son coches y columnas K sub algo*/
    int i_col = int(t_c / dt);
    double v_norm = v * m / C;

    V[0][0] = v_norm;            /*Condiciones iniciales del coche 0*/
    X[0][0] = 0;

    K[0][0] = 0;
    K[0][1] = 0 + dt / 2 * K[0][0];    /*K's del coche 0*/
    K[0][2] = 0 + dt / 2 * K[0][1];
    K[0][3] = dt * K[0][2];



    for (int j = 1; j < 5; j++){
        X[0][j] = X[0][j-1] - (l+d) / d;             /*Posiciones iniciales de los otros 4*/
        V[0][j] = v_norm;                            /*Velocidades iniciales de los otros 4*/
    }

    for (int i=0; i<N-1; i++){
        //Antes del canvio de velocidad
        if (i <= i_col){
            X[i+1][0] = X[i][0] + v_norm * dt;
            V[i+1][0] = v_norm;

            L[0][0] = v_norm;
            L[0][1] = v_norm + dt / 2 * K[0][0];
            L[0][2] = v_norm + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_norm + dt * K[0][2];

            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i][j-1])/(fabs(X[i][j] - X[i][j-1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i][j-1] + dt*0.5*L[j-1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i][j-1] + dt*0.5*L[j-1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i][j-1] + dt*L[j-1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
            //Despues del cambio de velocidad
        else {
            double v_new = v_norm*( 1 - (dt*i - t_c)*pow(M_E,1 - (dt*i - t_c)) );
            X[i+1][0] = X[i][0] + v_new * dt;
            V[i+1][0] = v_new;

            K[0][0] = v_norm * (i*dt - t_c - 1) * pow(M_E, 1-(dt*i-t_c));
            K[0][1] = v_norm * (i*dt - t_c - 1) * pow(M_E, 1-(dt*i-t_c)) + dt/2 * K[0][0];
            K[0][2] = v_norm * (i*dt - t_c - 1) * pow(M_E, 1-(dt*i-t_c)) + dt/2 * K[0][1];
            K[0][3] = v_norm * (i*dt - t_c - 1) * pow(M_E, 1-(dt*i-t_c)) * dt * K[0][2];

            L[0][0] = v_new;
            L[0][1] = v_new + dt / 2 * K[0][0];
            L[0][2] = v_new + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_new + dt * K[0][2];

            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i-reaccion][j-1])/(fabs(X[i][j] - X[i-reaccion][j-1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i-reaccion][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i-reaccion][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i-reaccion][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i-reaccion][j-1] + dt*L[j-1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
    }

    FILE* data;
    data = fopen("Instantaneo2.csv", "w");
    fprintf(data, "t (s), x_0 (m), x_1 (m), x_2 (m), x_3 (m), x_4 (m)\n");

    for (int i = 0; i < N - 1; i++){
        fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, X[i][0]*d,X[i][1]*d,X[i][2]*d,X[i][3]*d,X[i][4]*d);
    }
    fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, X[N-1][0]*d,X[N-1][1]*d,X[N-1][2]*d,X[N-1][3]*d,X[N-1][4]*d);

    FILE* data2;
    data2 = fopen("Velocidades2.csv", "w");
    fprintf(data2, "t (s), v_0 (m/s), v_1 (m/s), v_2 (m/s), v_3 (m/s), v_4 (m/s)\n");
    for (int i = 0; i < N - 1; i++){
        fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, V[i][0]* C/m,V[i][1]* C/m,V[i][2]* C/m,V[i][3]* C/m,V[i][4]* C/m);
    }
    fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, V[N-1][0]* C/m,V[N-1][1]* C/m,V[N-1][2]* C/m,V[N-1][3]* C/m,V[N-1][4]* C/m);

}

void movimiento_sinusoidal(int reaccion, double omega){
    double X[N][5];       /* Fila tiempo, columna coche*/
    double V[N][5];       /* Fila tiempo, columna coche*/
    double K[5][4], L[5][4];          /*Filas son coches y columnas K sub algo*/
    int i_col = int(t_c / dt);
    double v_norm = v * m / C;

    V[0][0] = v_norm;            /*Condiciones iniciales del coche 0*/
    X[0][0] = 0;

    K[0][0] = 0;
    K[0][1] = 0 + dt / 2 * K[0][0];    /*K's del coche 0*/
    K[0][2] = 0 + dt / 2 * K[0][1];
    K[0][3] = dt * K[0][2];



    for (int j = 1; j < 5; j++){
        X[0][j] = X[0][j-1] - (l+d) / d;             /*Posiciones iniciales de los otros 4*/
        V[0][j] = v_norm;                            /*Velocidades iniciales de los otros 4*/
    }

    for (int i=0; i<N-1; i++){
        //Antes del canvio de velocidad
        if (i <= i_col){
            X[i+1][0] = X[i][0] + v_norm * dt;
            V[i+1][0] = v_norm;

            L[0][0] = v_norm;
            L[0][1] = v_norm + dt / 2 * K[0][0];
            L[0][2] = v_norm + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_norm + dt * K[0][2];

            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i][j-1])/(fabs(X[i][j] - X[i][j-1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i][j-1] + dt*0.5*L[j-1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i][j-1] + dt*0.5*L[j-1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i][j-1] + dt*L[j-1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
            //Despues del cambio de velocidad
        else {
            double v_new = v_norm * (1-0.8*pow(sin(omega*(i*dt-t_c)),2));
            double a_new = v_norm*(-1.6*sin(omega*(i*dt-t_c))*cos(omega*(i*dt-t_c))*omega);
            X[i+1][0] = X[i][0] + v_new * dt;
            V[i+1][0] = v_new;

            K[0][0] = a_new;
            K[0][1] = a_new + dt/2 * K[0][0];
            K[0][2] = a_new + dt/2 * K[0][1];
            K[0][3] = a_new + dt * K[0][2];
            L[0][0] = v_new;
            L[0][1] = v_new + dt / 2 * K[0][0];
            L[0][2] = v_new + dt / 2 * K[0][1];            /*L's del coche 0*/
            L[0][3] = v_new + dt * K[0][2];

            for (int j=1; j<5; j++){
                K[j][0] = -(V[i][j] - V[i-reaccion][j-1])/(fabs(X[i][j] - X[i-reaccion][j-1]));
                L[j][0] = V[i][j];

                K[j][1] = -((V[i][j] + dt*0.5*K[j][0]) - (V[i-reaccion][j-1] + dt*0.5*K[j-1][0]))/(fabs((X[i][j] + dt*0.5*L[j][0]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][0])));
                L[j][1] = V[i][j] + dt*0.5*K[j][0];

                K[j][2] = -((V[i][j] + dt*0.5*K[j][1])- (V[i-reaccion][j-1] + dt*0.5*K[j-1][1]))/(fabs((X[i][j] + dt*0.5*L[j][1]) - (X[i-reaccion][j-1] + dt*0.5*L[j-1][1])));
                L[j][2] = V[i][j] + dt*0.5*K[j][1];

                K[j][3] = -((V[i][j] + dt*K[j][2]) - (V[i-reaccion][j-1] + dt*K[j-1][2]))/(fabs((X[i][j] + dt*L[j][2]) - (X[i-reaccion][j-1] + dt*L[j-1][2])));
                L[j][3] = V[i][j] + dt*K[j][2];
            }

            for (int k=1; k < 5; k++){
                X[i+1][k] = X[i][k] + (dt/6)*(L[k][0] + 2*L[k][1] + 2*L[k][2] + L[k][3]);
                V[i+1][k] = V[i][k] + (dt/6)*(K[k][0] + 2*K[k][1] + 2*K[k][2] + K[k][3]);
            }
        }
    }

    FILE* data;
    data = fopen("Instantaneo3.csv", "w");
    fprintf(data, "t (s), x_0 (m), x_1 (m), x_2 (m), x_3 (m), x_4 (m)\n");

    for (int i = 0; i < N - 1; i++){
        fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, X[i][0]*d,X[i][1]*d,X[i][2]*d,X[i][3]*d,X[i][4]*d);
    }
    fprintf(data, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, X[N-1][0]*d,X[N-1][1]*d,X[N-1][2]*d,X[N-1][3]*d,X[N-1][4]*d);

    FILE* data2;
    data2 = fopen("Velocidades3.csv", "w");
    fprintf(data2, "t (s), v_0 (m/s), v_1 (m/s), v_2 (m/s), v_3 (m/s), v_4 (m/s)\n");
    for (int i = 0; i < N - 1; i++){
        fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf\n", dt*i/t_0, V[i][0]* C/m,V[i][1]* C/m,V[i][2]* C/m,V[i][3]* C/m,V[i][4]* C/m);
    }
    fprintf(data2, "%lf, %lf, %lf, %lf, %lf, %lf", dt*(N-1)/t_0, V[N-1][0]* C/m,V[N-1][1]* C/m,V[N-1][2]* C/m,V[N-1][3]* C/m,V[N-1][4]* C/m);

}

int main(){
    double reaction_time = 0;
    int reaction_steps = floor(reaction_time * t_0 / dt);
    double omega = 1. / 4 * M_PI;
    cout << "El tiempo de reaccion es de " << reaction_time << " s \n";
    respuesta_instantanea(reaction_steps);
    movimiento_sinusoidal(reaction_steps, omega);

    return 0;
}