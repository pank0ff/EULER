#include <iostream>
#include <cmath>
#include "Expression.h"
#include "ExplicitEulerSolver.h"
#include "ImplicitEulerSolver.h"

#define DIM 2
#define DIM2 3

#define L1 -100.0
#define L2 -200.0
#define L3 -300.0

using namespace std;

//Чтобы задать свою функцию необходимо отнаследоваться от Expression (последний элемент массива args - параметр t)

class F1 : public Expression {
    double F(double* args) override {
        return -args[0] * args[1] + (abs(args[2]) < 10e-6 ? 1 : sin(args[2]) / args[2]);
    }
};
class F2 : public Expression {
    double F(double* args) override {
        return -args[1] * args[1] + 3.125 * args[2] / (1 + args[2] * args[2]);
    }
};

class F3 : public Expression {
public:
    double F(double* args) override {
        return ((2 * L1 + 4 * L2) * args[0] + 2 * (L1 - L2) * args[1] + 2 * (L1 - L2) * args[2] + (4 * L1 + 2 * L2)) / 6;
    }
};
class F4 : public Expression {
public:
    double F(double* args) override {
        return (2 * (L1 - L2) * args[0] + (2 * L1 + L2 + 3 * L3) * args[1] + (2 * L1 + L2 - 3 * L3) * args[2] + (4 * L1 - L2 - 9 * L3)) / 6;
    }
};
class F5 : public Expression {
public:
    double F(double* args) override {
        return (2 * (L1 - L2) * args[0] + (2 * L1 + L2 - 3 * L3) * args[1] + (2 * L1 + L2 + 3 * L3) * args[2] + (4 * L1 - L2 + 9 * L3)) / 6;
    }
};


int main() {
    //Test 1

    Expression** funcs = new Expression * [DIM];
    funcs[0] = new F1;
    funcs[1] = new F2;

    double* initApprox = new double[DIM];
    initApprox[0] = 0;
    initApprox[1] = -0.412;

    ExplicitEulerSolver solver1(DIM, funcs, 0, 1, 1e-3, 0.001, initApprox);
    solver1.solve();

    //Test 2

    Expression** f = new Expression * [DIM2];
    f[0] = new F3;
    f[1] = new F4;
    f[2] = new F5;

    double* init = new double[DIM2];
    init[0] = 10;
    init[1] = 22;
    init[2] = 9;

    ImplicitEulerSolver solver2(DIM2, f, 0, 1, 1e-5, 1e-3, 1e-5, init);
    solver2.solve();

    return 0;
}