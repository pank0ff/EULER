#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//������� ������� �������� ��������� ������� ������ � LDLT ������������
class GaussianSolver {
private:
    double** initMatrix;
    double** curMatrix;
    double* solutions;
    int dim;
    int dimExp;

public:
    //������: ������ ����� - ������� �������, ����� ����� ������� �������� �������
    explicit GaussianSolver(ifstream& file);

    GaussianSolver(double** matrix, int dimension);

    void solve();               //����� ������

    void solveViaLDLT();        //����� LDLT ������������ (������ ��� ������������ ������)

    bool isSymmetric();

    double* getSolutions();

    void printInitial();

    void printSolutions();

    double getDiscrepancy();    //���������� ����� ������� �������

    double getRelativeError();  //���������� ������������� �����������

    ~GaussianSolver();

};
