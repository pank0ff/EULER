#pragma once
#include <iostream>
#include <iomanip>
#include "Expression.h"
#include "NewtonSolver.h"

using namespace std;

//Решение систем дифференциальных уравнений неявным методом Эйлера
class ImplicitEulerSolver {
private:

    //Класс - обертка для Expression
    //Позволяет использовать NewtonSolver для решения уравнений, вид которых не известен на этапе компиляции
    //Это обеспечивает единый интерфейс для явного и неявного методов (Задание уравнений наследованием от Expression)
    class DifferentialExpression : public Expression {
        Expression* function;       //Само уравнение
        int dim;                    //Кол-во уравнений
        int i;                      //Номер уравнения
        double* solutions;          //Вектор параметров Yk, вычисляются по ходу выполнения, но не изменяются в NewtonSolver
        double tau;
        double t;
        double* argsWithT;          //Расширенный вектор аргументов (последний элемент - параметр t)

        double adjust(double* args) {
            return args[i] - tau * function->F(args) - solutions[i];    //Уравнение которое нужно решить в NewtonSolver
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

        //Переопределнный метод, который вызывается в NewtonSolver, где ему передаются аргументы без параметра
        double F(double* args) override {
            //Добавляем к ним параметр t
            for (int j = 0; j < dim; ++j) {
                argsWithT[j] = args[j];
            }
            argsWithT[dim] = t;
            return adjust(argsWithT);
        }

        //Измененная в соответствии с F производная данной функции
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