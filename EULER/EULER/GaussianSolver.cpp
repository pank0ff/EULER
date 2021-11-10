#include "GaussianSolver.h"

using namespace std;

GaussianSolver::GaussianSolver(ifstream& file) {
    if (!file.is_open()) {
        cerr << "FILE NOT FOUND" << endl;
        throw exception();
    }
    dim = 0;
    file >> dim;
    dimExp = dim + 1;
    initMatrix = new double* [dim];
    curMatrix = new double* [dim];
    for (int i = 0; i < dim; ++i) {
        initMatrix[i] = new double[dimExp];
        curMatrix[i] = new double[dimExp];
        for (int j = 0; j < dimExp; ++j) {
            file >> initMatrix[i][j];
            curMatrix[i][j] = initMatrix[i][j];
        }
    }
    solutions = nullptr;
}

GaussianSolver::GaussianSolver(double** matrix, int dimension) {
    dim = dimension;
    dimExp = dim + 1;
    initMatrix = new double* [dim];
    curMatrix = new double* [dim];
    for (int i = 0; i < dim; ++i) {
        initMatrix[i] = new double[dimExp];
        curMatrix[i] = new double[dimExp];
        for (int j = 0; j < dimExp; ++j) {
            initMatrix[i][j] = matrix[i][j];
            curMatrix[i][j] = initMatrix[i][j];
        }
    }
    solutions = nullptr;
}

void GaussianSolver::solve() {
    for (int k = 0; k < dim; ++k) {
        int iMax = k;                       //Поиск строки с максимальным жлементом и перестановка
        double vMax = curMatrix[k][k];
        for (int i = k; i < dim; ++i) {
            if (abs(curMatrix[i][k]) > abs(vMax)) {
                iMax = i;
                vMax = curMatrix[i][k];
            }
        }
        if (vMax == 0) {
            cerr << "SINGULAR MATRIX" << endl;
            throw exception();
        }
        if (iMax != k) {
            double* buf = curMatrix[k];
            curMatrix[k] = curMatrix[iMax];
            curMatrix[iMax] = buf;
        }

        for (int j = k; j < dimExp; ++j) {      //Прямой ход Гаусса
            curMatrix[k][j] /= vMax;
        }
        double factor;
        for (int i = k + 1; i < dim; ++i) {
            factor = curMatrix[i][k];
            for (int j = k; j < dimExp; ++j) {
                curMatrix[i][j] = curMatrix[i][j] - curMatrix[k][j] * factor;
            }
        }
    }

    //Обратных ход Гаусса
    solutions = new double[dim];
    for (int k = dim - 1; k >= 0; --k) {
        solutions[k] = curMatrix[k][dimExp - 1];
        for (int j = k + 1; j < dim; ++j) {
            solutions[k] -= curMatrix[k][j] * solutions[j];
        }
    }
}

void GaussianSolver::solveViaLDLT() {
    if (!isSymmetric()) {
        cerr << "CANT FACTORIZE NON-SYMMETRIC MATRIX" << endl;
        throw exception();
    }

    double** ATilda = new double* [dim];        //Матрица вспомогатльных величин
    for (int i = 0; i < dim; ++i) {
        ATilda[i] = new double[dim];
    }
    double** L = new double* [dim];             //Нижняя треульная матрица L
    for (int i = 0; i < dim; ++i) {
        L[i] = new double[dim];
        for (int j = 0; j < dim; ++j) {
            if (i == j) {
                L[i][j] = 1;
            }
            else if (i < j) {
                L[i][j] = 0;
            }
        }
    }
    double* D = new double[dim];                //"Сжатая" матрица диагональная матрица D

    for (int j = 0; j < dim; ++j) {             //Формирование матриц L и D
        double d = initMatrix[j][j];
        for (int k = 0; k < j; ++k) {
            d -= ATilda[j][k] * L[j][k];
        }
        D[j] = d;
        for (int i = j + 1; i < dim; ++i) {
            double aTilda = initMatrix[i][j];
            for (int k = 0; k < j; ++k) {
                aTilda -= ATilda[i][k] * L[j][k];
            }
            ATilda[i][j] = aTilda;

            L[i][j] = ATilda[i][j] / D[j];
        }
    }

    solutions = new double[dim];                //Решение трех простых систем линейных уравнений
    for (int i = 0; i < dim; ++i) {
        solutions[i] = initMatrix[i][dimExp - 1];
        for (int k = 0; k < i; ++k) {
            solutions[i] -= L[i][k] * solutions[k];
        }
    }

    for (int i = 0; i < dim; ++i) {
        solutions[i] /= D[i];
    }

    for (int i = dim - 1; i >= 0; --i) {
        for (int k = i + 1; k < dim; ++k) {
            solutions[i] -= L[k][i] * solutions[k];
        }
    }

    for (int i = 0; i < dim; ++i) {
        delete[] ATilda[i];
        delete[] L[i];
    }
    delete[] ATilda;
    delete[] L;
    delete[] D;
}

bool GaussianSolver::isSymmetric() {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dim; ++j) {
            if (initMatrix[i][j] != initMatrix[j][i]) {
                return false;
            }
        }
    }
    return true;
}

double* GaussianSolver::getSolutions() {
    return solutions;
}

void GaussianSolver::printInitial() {
    for (int i = 0; i < dim; ++i) {
        for (int j = 0; j < dimExp; ++j) {
            if (j == dimExp - 1) {
                cout << setw(10) << "|";
            }
            cout << setw(10) << initMatrix[i][j];
        }
        cout << endl;
    }
    cout << endl;
}

void GaussianSolver::printSolutions() {
    for (int i = 0; i < dim; ++i) {
        cout << "X(" << i << ") = " << solutions[i] << endl;
    }
    cout << endl;
}

double GaussianSolver::getDiscrepancy() {
    double maxDelta = 0;
    for (int i = 0; i < dim; ++i) {
        double result = initMatrix[i][dimExp - 1];
        for (int j = 0; j < dim; ++j) {
            result -= initMatrix[i][j] * solutions[j];
        }
        result = abs(result);
        if (result >= maxDelta) {
            maxDelta = result;
        }
    }
    return maxDelta;
}

double GaussianSolver::getRelativeError() {
    double** auxMatrix = new double* [dim];             //Формирование вспомогательной системы
    for (int i = 0; i < dim; ++i) {
        auxMatrix[i] = new double[dimExp];
        for (int j = 0; j < dim; ++j) {
            auxMatrix[i][j] = initMatrix[i][j];
        }
        double sum = 0;
        for (int j = 0; j < dim; ++j) {
            sum += initMatrix[i][j] * solutions[j];
        }
        auxMatrix[i][dimExp - 1] = sum;
    }

    GaussianSolver auxSolver(auxMatrix, dim);           //Решение вспомогательной системы
    auxSolver.solve();
    double* auxSolutions = auxSolver.getSolutions();

    double maxSolDelta = abs(solutions[0] - auxSolutions[0]);   //Вычисление погрешности
    double maxSol = abs(solutions[0]);
    for (int i = 1; i < dim; ++i) {
        if (abs(solutions[i]) > maxSol) {
            maxSol = abs(solutions[i]);
        }
        if (abs(solutions[i] - auxSolutions[i]) > maxSolDelta) {
            maxSolDelta = abs(solutions[i] - auxSolutions[i]);
        }
    }

    for (int i = 0; i < dim; ++i) {
        delete[] auxMatrix[i];
    }
    delete[] auxMatrix;

    delete[] auxSolutions;

    return maxSolDelta / maxSol;
}

GaussianSolver::~GaussianSolver() {
    for (int i = 0; i < dim; ++i) {
        delete[] curMatrix[i];
        delete[] initMatrix[i];
    }
    delete[] curMatrix;
    delete[] initMatrix;
    delete[] solutions;
}
