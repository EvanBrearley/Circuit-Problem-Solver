#ifndef CIRCUITSOLVER_H
#define CIRCUITSOLVER_H

typedef struct {
  double current, voltage;
  int direction, polarity, node1, node2;
} generalElement;

typedef struct {
  int node1, node2;
  double resistance;
} resistors;

typedef struct {
  int node1, node2, polarity;
  double voltage;
} voltsrc;

typedef struct {
  int node1, node2, direction;
  double current;
} currentsrc;

int epsilon;

void KCL();
void KVL();
void nodeAnalysis();
void meshAnalysis();
void swap (int *a, int *b);
void floatSwap (double *a, double *b);
void AnalyseMesh(int meshnum, int nodenum, int *meshes, int voltageNum, voltsrc voltagesource[], int resistorNum, resistors resistor[], double meshcurrent[], int mesh, double* constant);
void analyseNode(double nodeVoltage[], resistors resistor[], voltsrc voltagesource[], currentsrc currentsource[], int resistorNum, int currentNum, int voltageNum, int i);
void saveEqn(double *constant, int meshnum, int mesh, double meshcurrent[], double *resultantMatrix, double *constantVector);
void saveNodeEqn(double *constant, int nodeNum, int row, double variable[], double *resultantMatrix, double *constantVector);
int checkOtherMeshes(int meshnum, int nodenum, int *meshes, int node1, int node2, int mesh, int direction, int* sameDirection);
void addRow (int row1, int row2, int N, double Matrix[N][N+1]);
void scalarMultiplyRow(int row, double multiplier, int N, double Matrix[N][N+1]);
void switchRows (int row1, int row2, int N, double Matrix[N][N+1]);
void checkForZeros (int column, int N, double Matrix[N][N+1]);
void solveMatrix (double* coefficientMatrix, double* constants, double* solutions, int N);
void printMatrix (int N, double Matrix[N][N+1]);

#endif