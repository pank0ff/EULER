#pragma once
#include <fstream>
#include <iostream>
#include <iomanip>
#include <cmath>

using namespace std;

//Решение системы линейных уравнений методом Гаусса и LDLT факторизации
class GaussianSolver {
private:
    double** initMatrix;
    double** curMatrix;
    double* solutions;
    int dim;
    int dimExp;

public:
    //Формат: Первое число - порядок матрицы, далее через проебел элементы матрицы
    explicit GaussianSolver(ifstream& file);

    GaussianSolver(double** matrix, int dimension);

    void solve();               //Метод Гаусса

    void solveViaLDLT();        //Метод LDLT факторизации (только для симметричных матриц)

    bool isSymmetric();

    double* getSolutions();

    void printInitial();

    void printSolutions();

    double getDiscrepancy();    //Вычисление нормы вектора невязки

    double getRelativeError();  //Вычисление относительной погрешности

    ~GaussianSolver();

};
