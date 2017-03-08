#include <stdio.h>
#include <assert.h>
#include "CircuitSolver.h"

int epsilon = 0.0000001;

void addRow (int row1, int row2, int N, double Matrix[N][N+1]){
    for (int i = 0; i <= N; i++){                                //have to account for constants in this, N is limit of coefficient matrix//
        Matrix[row1][i] += Matrix[row2][i];
    }
}

void scalarMultiplyRow(int row, double multiplier, int N, double Matrix[N][N+1]){
    for (int i = 0; i <= N; i++){                                           //have to account for constants in this, N is limit of coefficient matrix//
        Matrix[row][i] *= multiplier;
    }
}

void switchRows (int row1, int row2, int N, double Matrix[N][N+1]){
    for (int i = 0; i <= N; i++){
        double temp = Matrix[row1][i];
        Matrix[row1][i] = Matrix[row2][i];
        Matrix[row2][i] = temp;
    }
}

void checkForZeros (int column, int N, double Matrix[N][N+1]){ //ensures there are no zeros on the main diagonal//
    int epsilon = 0.0000001;
    for (int i = column; i < N; i++){
        if (Matrix[i][column] > epsilon){
            switchRows(column, i, N, Matrix);
            break;
        }
    }
    assert(Matrix[column][column] > epsilon || Matrix[column][column] * -1 > epsilon);
}

void solveMatrix (double* coefficientMatrix, double* constants, double* solutions, int N){
    double augmentedMatrix[N][N+1];
    for (int i = 0; i < N; i++){
        for (int j = 0; j < N; j++){
            augmentedMatrix[i][j] = coefficientMatrix[i*N + j];
        }
    }
    for (int i = 0; i < N; i++){
        augmentedMatrix[i][N] = constants[i];
    }
    printMatrix(N, augmentedMatrix);
    for (int i = 0; i < N; i++){     //loops through coefficients, columns//
        checkForZeros(i, N, augmentedMatrix);
        for (int j = i+1; j < N; j++){
            if (augmentedMatrix[j][i]){
                scalarMultiplyRow(i, (-1*augmentedMatrix[j][i]/augmentedMatrix[i][i]), N, augmentedMatrix);
                addRow(j, i, N, augmentedMatrix);
            }
        }  
    }
    printMatrix(N, augmentedMatrix);
    //Matrix should now be upper diagonal with no zeros on main diagonal, checkForZeros assert would have triggered//
    for (int i = N - 1; i > -1; i--){
        for (int j = i-1; j > -1; j--){
            if (augmentedMatrix[j][i]){
                scalarMultiplyRow(i, (-1*augmentedMatrix[j][i]/augmentedMatrix[i][i]), N, augmentedMatrix);
                addRow(j, i, N, augmentedMatrix);
            }
        }
    }
    //Should now be diagonal//
    for (int i = 0; i < N; i++){
        assert(augmentedMatrix[i][i] > epsilon);
        scalarMultiplyRow(i, 1/augmentedMatrix[i][i], N, augmentedMatrix);
    }
    printMatrix(N, augmentedMatrix);
    //last column should be answers//
}

void printMatrix (int N, double Matrix[N][N+1]){
    for (int i = 0; i < N; i++){
        printf("[ ");
        for (int j = 0; j < N; j++){
            if (j == N-1){
                printf("%g", Matrix[i][j]);
            }
            else{
                printf("%g ", Matrix[i][j]);
            }
        }
        printf(" | %g]\n", Matrix[i][N]);
    }
}