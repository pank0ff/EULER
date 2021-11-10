#include "Expression.h"

double Expression::D(int curArg, int argsAmount, double* args) {
    double* augArgs = new double[argsAmount];
    for (int l = 0; l < argsAmount; ++l) {
        if (l != curArg) {
            augArgs[l] = args[l];
        }
        else {
            augArgs[l] = args[l] + DERIVATIVE_DX;
        }
    }
    double result = (F(augArgs) - F(args)) / DERIVATIVE_DX;
    delete[] augArgs;
    return result;
}

Expression::~Expression() { }
