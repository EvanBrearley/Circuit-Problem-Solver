#include <stdio.h>
#include <assert.h>
#include "CircuitSolver.h"

int main (void) {
    int functionType;
    printf("Enter 1 for KCL, 2 for KVL,3 - Nodal Analysis, 4 - Mesh Analysis\n");
    scanf("%d", &functionType);
    if (functionType == 1){
        KCL();
    }
    if (functionType == 2){
        KVL();
    }
    if (functionType == 3){
        nodeAnalysis();
    }
    if (functionType == 4){
        meshAnalysis();
    }
    return 0;
}

void KCL () {
    int nodeNum;
    double constant = 0;
    int knownNum = 0;
    int varnum = 0;

    printf("How many nodes analysed?");
    scanf("%d", &nodeNum);
    
    int node[nodeNum];
    
    printf("Which nodes would you like analysed?");
    for (int i = 0; i < nodeNum; i++) {
        scanf("%d", &node[i]);
    }

    double CoefficientMatrix[nodeNum][nodeNum+1]; //extra column for constants//
    double y[nodeNum]; //holds constants//
    double x[nodeNum]; //for the solution//

    //setup unknowns
    printf("How many variables\n");
    scanf("%d", &varnum);
    generalElement variable[varnum];
    printf("Their directions?\n");
    for (int i = 0; i < varnum; i++) {
        scanf("%d", &variable[i].direction);
        variable[i].current = 0;
    }
    //setup knowns//
    printf("How many known currents\n");
    scanf("%d", &knownNum);
    generalElement known[knownNum];
    printf("Type in their connected nodes, directions and what nodes they enter\n");
    for (int i = 0; i < knownNum; i++){
        scanf("%d %d %d %lf", &known[i].node1, &known[i].node2, &known[i].direction, &known[i].current);
    }

    for (int j = 0; j < nodeNum; j++) {
        constant = 0;
        for (int i = 0; i < varnum; i++) {
            if (variable[i].node1 == node[j] || variable[i].node2 == node[j]) {
                variable[i].current = 1;
            }
        }
        for (int i = 0; i < varnum; i++) {
            if (variable[i].direction != node[j]) {
                 variable[i].current *= -1;
            }
        }

        for (int i = 0; i < knownNum; i++) {
            if ((known[i].node1 == node[j]) || (known[i].node2 == node[j])) {
                if (known[i].direction == node[j]) {
                    constant -= known[i].current;
                }
                else {
                    constant += known[i].current;
                }
            }
        }

        for (int i = 0; i < varnum; i++) {
            CoefficientMatrix[j][i] = variable[i].current;
            variable[i].current = 0;
        }
        y[j] = constant;
    }
    solveMatrix((double *)CoefficientMatrix, (double *)y, (double *)x, nodeNum);
}

void KVL() {
    int meshNum; 
    printf("How many meshes are being analysed?");
    scanf("%d", &meshNum);
    double CoefficientMatrix[meshNum][meshNum]; //assuming meshNum = varnum, #equations = #variables//
    double y[meshNum]; //holds constants//
    double x[meshNum]; //for the solution//
    double constant;
    int  knownNum, varnum, nodenum;

    printf("How many variables\n");
    scanf("%d", &varnum);
    generalElement variable[varnum];
    printf("Enter the two connected nodes, and the polarities\n");
    for (int i = 0; i < varnum; i++) {
        scanf("%d %d %d", &variable[i].node1, &variable[i].node2, &variable[i].polarity);
        variable[i].voltage = 0;
    }

    printf("How many known voltages?\n");
    scanf("%d", &knownNum);
    generalElement known[knownNum]; //number of known voltages, rename//
    printf("Enter the two connected nodes, their polarities, and their voltages.\n");
    for (int i = 0; i < knownNum; i++) {
        scanf("%d %d %d %lf", &known[i].node1, &known[i].node2, &known[i].polarity, &known[i].voltage);
    }

    for (int mesh = 0; mesh < meshNum; mesh++) {
        printf("For Mesh %d\n", mesh + 1);
        constant = 0;

        printf("# Nodes in this mesh?");
        scanf("%d", nodenum);
        int node[nodenum + 1];

        printf("Enter the nodes that comprise the desired loop, in order, enter the first node at the end as well.");
        for (int i = 0; i < nodenum; i++){
            scanf("%d", &node[i]);
        }
        for (int i = 0; i < varnum; i++) {
            for (int j = 0; j < nodenum; j++) {
                if (variable[i].node1 == node[j] || variable[i].node2 == node[j]) {
                    if (variable[i].node1 == node[j + 1] || variable[i].node2 == node[j + 1]) {
                        variable[i].voltage = 1;
                    }
                }
            }
        }
        for (int i = 0; i < varnum; i++) {
            for (int j = 0; j < nodenum; j++) {
                if (variable[i].node2 == node[j]) {
                    if (variable[i].node1 == node[j + 1]) {
                        swap(&variable[i].node1, &variable[i].node2);
                    }
                }
            }
        }
        for (int i = 0; i < nodenum; i++) {
            for (int j = 0; j < varnum; j++) {
                if (variable[j].node1 == node[i]) {
                    if (node[i] == variable[j].polarity) {
                        variable[j].voltage *= -1;
                    }
                }
            }
        }
        for (int i = 0; i <  knownNum; i++) {
            for (int j = 0; j < nodenum; j++) {
                if (known[i].node2 == node[j]) {
                    if (known[i].node1 == node[j + 1]) {
                        swap(&known[i].node1, &known[i].node2);
                    }
                }
            }
        }
        for (int i = 0; i < nodenum; i++) {
            for (int j = 0; j <  knownNum; j++) {
                if (known[j].node1 == node[i] && known[j].node2 == node[i + 1]) {
                    if (node[i] == known[j].polarity) {
                        constant += known[j].voltage;
                    }
                    else {
                        constant -= known[j].voltage;
                    }
                }
            }
        }
        for (int i = 0; i < varnum; i++) {
            CoefficientMatrix[mesh][i] = variable[i].voltage;
            variable[i].voltage = 0;
        }
        y[mesh] = constant;
        constant = 0;
    }
    solveMatrix((double *)CoefficientMatrix, (double *)y, (double *)x, meshNum);
}

void nodeAnalysis() {
    int resistorNum, voltageNum, currentNum, nodenum;
    printf("How many nodes in the circuit?");
    scanf("%d", &nodenum);
    double CoefficientMatrix[nodenum - 1][nodenum - 1]; //assuming meshNum = varnum, #equations = #variables//
    double y[nodenum - 1]; //holds constants//
    double x[nodenum - 1]; //for the solution//

    double nodeVoltage[nodenum];

    for (int i = 0; i < nodenum; i++) {
        nodeVoltage[i] = 0; // nodeVoltage[0] used as the constant//
    }

    printf("How many resistors?");
    scanf("%d", &resistorNum);
    resistors resistor[resistorNum];
    for (int i = 0; i < resistorNum; i++) {
        printf("For %dth resistor, input the nodes it is connected to, and the resistance, node 0 is the ground, they must be sequential from zero.\n",i);
        scanf("%d %d %lf", &resistor[i].node1, &resistor[i].node2, &resistor[i].resistance);
    }

    printf("How many voltage sources?");
    scanf("%d", &voltageNum);
    voltsrc voltagesource[voltageNum];
    for (int i = 0; i < voltageNum; i++) {
        printf("For %dth voltage source, input the nodes it is connected to, the node with the higher potential,and the voltage.\n", i);
        scanf("%d %d %d %lf", &voltagesource[i].node1, &voltagesource[i].node2, &voltagesource[i].polarity, &voltagesource[i].voltage);
    }
    
    for (int i = 0; i < voltageNum; i++) {
        if (!voltagesource[i].node1) {
            swap(&voltagesource[i].node1, &voltagesource[i].node2);
        }
    }

    printf("How many current sources?");
    scanf("%d", &currentNum);
    currentsrc currentsource[currentNum];
    for (int i = 0; i < currentNum; i++) {
        printf("For %dth current source, input the nodes it is connected to, the direction of the current and its magnitude.\n", i);
        scanf("%d %d %d %lf", &currentsource[i].node1, &currentsource[i].node2, &currentsource[i].direction, &currentsource[i].current);
    }

    for (int i = 1; i < nodenum; i++) {
        int isVoltSource = 0;
        int analysed = 0;
        for (int j = 0; j < voltageNum; j++) {
            if (voltagesource[j].node1 == i) {
                isVoltSource++;
                if (voltagesource[j].node2) {
                    analyseNode(nodeVoltage, resistor, voltagesource, currentsource, resistorNum, currentNum, voltageNum, voltagesource[j].node2);
                    if (!analysed) {
                         analyseNode(nodeVoltage, resistor, voltagesource, currentsource, resistorNum, currentNum, voltageNum, i);
                         analysed++;
                    }
                }
            }
            if (voltagesource[j].node2 == i) {
                isVoltSource++;
            }
        }   

        if (analysed) {
            saveNodeEqn(&nodeVoltage[0], nodenum, i - 1, nodeVoltage, (double *)CoefficientMatrix, (double *)y);
            analysed = 0;
        }

        if (!isVoltSource) {
          analyseNode(nodeVoltage, resistor, voltagesource, currentsource, resistorNum, currentNum, voltageNum, i);
          saveNodeEqn(&nodeVoltage[0], nodenum, i - 1, nodeVoltage, (double *)CoefficientMatrix, (double *)y);
        }

        for (int j = 0; j < nodenum; j++) {
          nodeVoltage[j] = 0;
        }

        for (int j = 0; j < voltageNum; j++) {
            if (voltagesource[j].node1 == i) {
                if (voltagesource[j].node1 == voltagesource[j].polarity) {
                    nodeVoltage[voltagesource[j].node1] = 1;
                    nodeVoltage[voltagesource[j].node2] = -1;
                }
                else {
                    nodeVoltage[voltagesource[j].node2] = 1;
                    nodeVoltage[voltagesource[j].node1] = -1;
                }
                nodeVoltage[0] = voltagesource[j].voltage; //overwrites potential bug from nodeVoltage 0 getting a value//
                saveNodeEqn(&nodeVoltage[0], nodenum, voltagesource[j].node2 - 1, nodeVoltage, (double *)CoefficientMatrix, (double *)y);
            }
        }

        for (int j = 0; j < nodenum; j++) {
            nodeVoltage[j] = 0;
        }
    }
    solveMatrix((double *)CoefficientMatrix, (double *)y, (double *)x, nodenum - 1);
}

void analyseNode(double nodeVoltage[], resistors resistor[], voltsrc voltagesource[], currentsrc currentsource[], int resistorNum, int currentNum, int voltageNum, int i) {
    for (int j = 0; j < resistorNum; j++) {
        if (resistor[j].node2 == i) {
            swap(&resistor[j].node1, &resistor[j].node2);
        }
        if (resistor[j].node1 == i) {
            nodeVoltage[i] += (1 / resistor[j].resistance);
            if (resistor[j].node2) {
                nodeVoltage[resistor[j].node2] -= (1 / resistor[j].resistance);
            }
        }
    }
    for (int j = 0; j < currentNum; j++) {
        if (currentsource[j].node2 == i) {
            swap(&currentsource[j].node1, &currentsource[j].node2);
        }
        if (currentsource[j].node1 == i) {
            if (currentsource[j].direction != i) {
                currentsource[j].current *= -1;
            }
            nodeVoltage[0] += currentsource[j].current;
        }
    }
}

void meshAnalysis() {
    int resistorNum, voltageNum, currentNum, nodenum, meshnum;
    double CoefficientMatrix[meshnum][meshnum];
    double y[meshnum]; //holds constants//
    double x[meshnum]; //for the solution//

    
    printf("How many resistors?");
    scanf("%d", &resistorNum);
    resistors resistor[resistorNum];
    for (int i = 0; i < resistorNum; i++) {
        printf("For %dth resistor, input the nodes it is connected to, and the resistance, node 0 is the ground, they must be sequential from zero.\n",i);
        scanf("%d %d %lf", &resistor[i].node1, &resistor[i].node2, &resistor[i].resistance);
    }
    
    printf("How many voltage sources?");
    scanf("%d", &voltageNum);
    voltsrc voltagesource[voltageNum];
    for (int i = 0; i < voltageNum; i++) {
        printf("For %dth voltage source, input the nodes it is connected to, the node with the higher potential,and the voltage.\n", i);
        scanf("%d %d %d %lf", &voltagesource[i].node1, &voltagesource[i].node2, &voltagesource[i].polarity, &voltagesource[i].voltage);
    }
    
    printf("How many current sources?");
    scanf("%d", &currentNum);
    currentsrc currentsource[currentNum];
    for (int i = 0; i < currentNum; i++) {
        printf("For %dth current source, input the nodes it is connected to, the direction of the current and its magnitude.\n", i);
        scanf("%d %d %d %lf", &currentsource[i].node1, &currentsource[i].node2, &currentsource[i].direction, &currentsource[i].current);
    }

    printf("How many meshes are there, how many nodes are there in total and how many nodes are in the largest one?\n");
    scanf("%d %d", &meshnum, &nodenum);
    printf("Enter the nodes (0 can't be used) that comprise each mesh, with the first node also repeated at the end, enter 0 to advance to the next mesh\n");
    scanf("%d %d", &meshnum, &nodenum);
    int meshes[meshnum][nodenum + 1]; //extra slot to hold first node repeated at end//
    double meshcurrent[meshnum];
    for (int i = 0; i < meshnum; i++) {
        printf("for %dth mesh\n", i);
        for (int j = 0; j < nodenum+1; j++){
            scanf("%d", &meshes[i][j]);
            if (!meshes[i][j]) {
                break;
            }
        }
    }
    
    double constant;
    int eqnConsistent;
    for (int i = 0; i < meshnum; i++) {
        for (int i = 0; i < meshnum; i++) {
            meshcurrent[i] = 0;
        }
        int sameDirection = 0;
        int noCurrentSources = 1;
        constant = 0;
        eqnConsistent = 0;
        for (int j = 0; j < nodenum; j++) {
            if (!meshes[i][j + 1]) {
                break;
            }
            for (int k = 0; k < currentNum; k++) {
                if (currentsource[k].node2 == meshes[i][j] && currentsource[k].node1 == meshes[i][j + 1]) {
                    swap(&currentsource[k].node1, &currentsource[k].node2);
                }
                if (currentsource[k].node1 == meshes[i][j] && currentsource[k].node2 == meshes[i][j + 1]) {
                    if (currentsource[k].direction == meshes[i][j]) {
                        sameDirection = 1;
                    }
                    else {
                        sameDirection = 0;
                    }
                    noCurrentSources--;
                    int otherMesh = checkOtherMeshes(meshnum, nodenum, (int *)meshes, currentsource[k].node1, currentsource[k].node2, i, currentsource[k].direction, &sameDirection);
                    if (otherMesh > -1) {
                        AnalyseMesh(meshnum, nodenum, (int *)meshes, voltageNum, voltagesource, resistorNum, resistor, meshcurrent, i, &constant);
                        AnalyseMesh(meshnum, nodenum, (int *)meshes, voltageNum, voltagesource, resistorNum, resistor, meshcurrent, otherMesh, &constant);
                        saveEqn(&constant, meshnum, i, meshcurrent, (double *)CoefficientMatrix, (double *)y);
                        meshcurrent[i] = 1;
                        meshcurrent[otherMesh] = 1;
                        if (!sameDirection) {
                            meshcurrent[otherMesh] *= -1;
                        }
                        constant = currentsource[k].current;
                        saveEqn(&constant, meshnum, otherMesh, meshcurrent, (double *)CoefficientMatrix, (double *)y);
                    }
                    else {
                        meshcurrent[i] = 1;
                        constant = currentsource[k].current;
                        if (currentsource[k].direction == meshes[i][j]) {
                            constant *= -1;
                        }
                        saveEqn(&constant, meshnum, i, meshcurrent, (double *)CoefficientMatrix, (double *)y);
                    }
                }
              }
        }
        if (noCurrentSources) {
            AnalyseMesh(meshnum, nodenum, (int *)meshes, voltageNum, voltagesource, resistorNum, resistor, meshcurrent, i, &constant);
            saveEqn(&constant, meshnum, i, meshcurrent, (double *)CoefficientMatrix, (double *)y);
        }
    }
    solveMatrix((double *)CoefficientMatrix, (double *)y, (double *)x, meshnum);
}

void AnalyseMesh(int meshnum, int nodenum, int *meshes, int voltageNum, voltsrc voltagesource[], int resistorNum, resistors resistor[], double meshcurrent[], int mesh, double *constant) {
    int i = mesh;
    for (int j = 0; j < nodenum; j++) {
        if (!meshes[i * nodenum + i + j + 1] || !meshes[i * nodenum + i + j]) {
            break;
        }
        for (int k = 0; k < voltageNum; k++) {
            if (voltagesource[k].node2 == meshes[i * nodenum + i + j]) { //prevent an accidental double trigger, by sorting first, then checking and breaking on first instance//
                if (voltagesource[k].node1 == meshes[i * nodenum + i + j + 1]) {
                    swap(&voltagesource[k].node1, &voltagesource[k].node2);
                }
            }
            if (voltagesource[k].node1 == meshes[i * nodenum + i + j]) {
                if (voltagesource[k].node2 == meshes[i * nodenum + i + j + 1]) {
                    if (voltagesource[k].node1 == voltagesource[k].polarity) {
                        *constant -= voltagesource[k].voltage;
                    }
                    else {
                        *constant += voltagesource[k].voltage;
                    }
                    break;
                }
            }
        }
        for (int k = 0; k < resistorNum; k++) {
            if (resistor[k].node2 == meshes[i * nodenum + i + j]) { //prevent an accidental double trigger, by sorting first, then checking and breaking on first instance//
                if (resistor[k].node1 == meshes[i * nodenum + i + j + 1]) {
                  swap(&resistor[k].node1, &resistor[k].node2);
                }
            }
            if (resistor[k].node1 == meshes[i * nodenum + i + j]) {
                if (resistor[k].node2 == meshes[i * nodenum + i + j + 1]) {
                    meshcurrent[i] += resistor[k].resistance;
                    int placeholderSameDirection = 0;
                    int otherMesh = checkOtherMeshes(meshnum, nodenum, (int *)meshes, resistor[k].node1, resistor[k].node2, i, 0, &placeholderSameDirection);
                    if (otherMesh > -1) {
                        meshcurrent[otherMesh] -= resistor[k].resistance;
                    }
                }
            }
        }
    }
}

int checkOtherMeshes(int meshnum, int nodenum, int *meshes, int node1, int node2, int mesh, int direction, int* sameDirection) {
    int otherDirection = 1;
    for (int i = 0; i < mesh; i++) { //check meshes before current i//
        for (int j = 0; j < nodenum; j++) {
            if (!meshes[i * nodenum + i + j + 1] || !meshes[i * nodenum + i + j]) {
                break;
            }
            if (node1 == meshes[i * nodenum + i + j] || node2 == meshes[i * nodenum + i + j]) {
                if (node2 == meshes[i * nodenum + i + j + 1] || node1 == meshes[i * nodenum + i + j + 1]) {
                    if (direction == meshes[i * nodenum + i + j]) {
                        otherDirection = 1;
                    }
                    else {
                        otherDirection = 0;
                    }
                    if (otherDirection != *sameDirection) {
                        *sameDirection = 0;
                    }
                    return i;
                }
            }
        }
    }
    for (int i = mesh + 1; i < meshnum; i++) { //check meshes before current 1//
        for (int j = 0; j < nodenum; j++) {
            if (!meshes[i * nodenum + i + j + 1] || !meshes[i * nodenum + i + j]) {
                break;
            }
            if (node1 == meshes[i * nodenum + i + j] || node2 == meshes[i * nodenum + i + j]) {
                if (node1 == meshes[i * nodenum + i + j + 1] || node2 == meshes[i * nodenum + i + j + 1]) {
                    if (direction == meshes[i * nodenum + i + j]) {
                        otherDirection = 1;
                    }
                    else {
                        otherDirection = 0;
                    }
                    if (otherDirection != *sameDirection) {
                        *sameDirection = 0;
                    }
                    return i;
                }
            }
        }
    }
    return -1;
}

void saveEqn(double *constant, int variableNum, int row, double variable[], double *resultantMatrix, double *constantVector) {
    for (int i = 0; i < variableNum; i++) {
        resultantMatrix[row * variableNum + i] = variable[i];
        variable[i] = 0;
    }
    constantVector[row] = *constant;
    *constant = 0;
}

void saveNodeEqn(double *constant, int nodeNum, int row, double variable[], double *resultantMatrix, double *constantVector) {
    for (int i = 1; i < nodeNum; i++) {
        resultantMatrix[row * nodeNum - row + i - 1] = variable[i];
        variable[i] = 0;
    }
    constantVector[row] = variable[0];
    variable[0] = 0;
}

void swap (int *a, int *b) {
    int temp = *a;
    *a = *b;
    *b = temp;
}