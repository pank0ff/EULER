#include "NewtonSolver.h"

NewtonSolver::NewtonSolver(int dim, Expression** functions, int kMax, double* initApprox, double eps) {
    this->dim = dim;
    this->functions = functions;
    this->kMax = kMax;
    solutions = new double[dim];
    for (int i = 0; i < dim; ++i) {
        solutions[i] = initApprox[i];
    }
    this->eps = eps;
    steps = 1;
}

void NewtonSolver::solve() {
    double** JF = new double* [dim];                 //ћатрица JF - матрица якоби
    for (int i = 0; i < dim; ++i) {
        JF[i] = new double[dim + 1];
    }

    double sig1;    // ритерии завершени€ вычислений
    double sig2;

    double* prevSolutions = new double[dim];
    while (true) {
        for (int i = 0; i < dim; ++i) {
            prevSolutions[i] = solutions[i];
        }

        for (int i = 0; i < dim; ++i) {             //«аполнени€ матрицы якоби
            for (int j = 0; j < dim + 1; ++j) {
                if (j < dim) {
                    JF[i][j] = functions[i]->D(j, dim, solutions);
                }
                else {
                    JF[i][j] = -functions[i]->F(solutions);
                }
            }
        }

        GaussianSolver solver(JF, dim);   //–ешение системы линейных уравнений с матрицей якоби
        solver.solve();
        double* dxk = solver.getSolutions();

        for (int i = 0; i < dim; ++i) {
            solutions[i] = solutions[i] + dxk[i];
        }

        sig1 = functions[0]->F(solutions);               //ѕервый критерий завершени€
        for (int i = 0; i < dim; ++i) {
            if (functions[i]->F(solutions) > sig1) {
                sig1 = functions[i]->F(solutions);
            }
        }

        //¬торой критерий завершени€
        double tmp = abs(solutions[0] - prevSolutions[0]);

        sig2 = abs(abs(solutions[0]) < 1 ? tmp : tmp / solutions[0]);
        for (int i = 1; i < dim; ++i) {
            tmp = abs(solutions[i] - prevSolutions[i]);
            if (abs(solutions[i]) < 1 and tmp > sig2) {
                sig2 = tmp;
            }
            else if (abs(solutions[i]) >= 1 and abs(tmp / solutions[i]) > sig2) {
                sig2 = abs(tmp / solutions[i]);
            }
        }

        steps++;

        if (steps >= kMax) {
            cerr << "WARNING! MAX STEPS REACHED!" << endl;
        }

        if ((sig1 <= eps && sig2 <= eps) || steps >= kMax) {
            break;
        }
    }


    delete[] prevSolutions;
    for (int i = 0; i < dim; ++i) {
        delete[] JF[i];
    }
    delete[] JF;
}

double* NewtonSolver::getSolutions() {
    return solutions;
}

int NewtonSolver::getStepsAmount() {
    return steps;
}

NewtonSolver::~NewtonSolver() {
    delete[] solutions;
}
