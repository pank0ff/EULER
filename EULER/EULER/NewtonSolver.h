#pragma once
#include <iomanip>
#include "Expression.h"
#include "GaussianSolver.h"

using namespace std;


//������� ������ ���������� ��������� ������� �������
class NewtonSolver {
private:
    int dim;
    Expression** functions;
    double* solutions;
    int steps;
    int kMax;
    double eps;

public:
    NewtonSolver(int dim, Expression** functions, int kMax, double* initApprox, double eps);

    void solve();

    double* getSolutions();

    int getStepsAmount();

    ~NewtonSolver();
}; 
