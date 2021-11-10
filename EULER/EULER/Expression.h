#pragma once
#define DERIVATIVE_DX 10e-9

class Expression {
public:

    virtual double F(double* args) = 0;

    virtual double D(int curArg, int argsAmount, double* args);

    virtual ~Expression();
};
