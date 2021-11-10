#pragma once
#include <iostream>
#include <iomanip>
#include "Expression.h"
#include "NewtonSolver.h"

using namespace std;

//������� ������ ���������������� ��������� ������� ������� ������
class ImplicitEulerSolver {
private:

    //����� - ������� ��� Expression
    //��������� ������������ NewtonSolver ��� ������� ���������, ��� ������� �� �������� �� ����� ����������
    //��� ������������ ������ ��������� ��� ������ � �������� ������� (������� ��������� ������������� �� Expression)
    class DifferentialExpression : public Expression {
        Expression* function;       //���� ���������
        int dim;                    //���-�� ���������
        int i;                      //����� ���������
        double* solutions;          //������ ���������� Yk, ����������� �� ���� ����������, �� �� ���������� � NewtonSolver
        double tau;
        double t;
        double* argsWithT;          //����������� ������ ���������� (��������� ������� - �������� t)

        double adjust(double* args) {
            return args[i] - tau * function->F(args) - solutions[i];    //��������� ������� ����� ������ � NewtonSolver
        }

    public:

        DifferentialExpression(int dim, Expression* function, int i, double* params, double tau, double x) {
            this->dim = dim;
            this->function = function;
            this->i = i;
            this->solutions = params;
            this->tau = tau;
            this->t = x;
            argsWithT = new double[dim + 1];
        }

        //��������������� �����, ������� ���������� � NewtonSolver, ��� ��� ���������� ��������� ��� ���������
        double F(double* args) override {
            //��������� � ��� �������� t
            for (int j = 0; j < dim; ++j) {
                argsWithT[j] = args[j];
            }
            argsWithT[dim] = t;
            return adjust(argsWithT);
        }

        //���������� � ������������ � F ����������� ������ �������
        double D(int curArg, int argsAmount, double* args) override {
            return (curArg == i ? 1 : 0) - tau * function->D(curArg, argsAmount, args);
        }

        void setTau(double tau) {
            this->tau = tau;
        }

        void setX(double x) {
            this->t = x;
        }

        virtual ~DifferentialExpression() {
            delete[] argsWithT;
        }
    };

    int dim;
    Expression** functions;
    double xBeg;
    double xEnd;
    double eps;
    double tauMax;
    double tauMin;
    double* solutions;
    double* solutionsPrev;
    double* solutionsNext;

    void print(double x);

public:
    ImplicitEulerSolver(int dim, Expression** functions, double xBeg, double xEnd, double eps, double tauMax, double tauMin, double* initApprox);

    void solve();

    ~ImplicitEulerSolver();
};