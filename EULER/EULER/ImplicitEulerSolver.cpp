#include "ImplicitEulerSolver.h"
#include "vector"
#include "pbPlots.hpp"
#include "supportLib.hpp"
ImplicitEulerSolver::ImplicitEulerSolver(int dim, Expression** functions, double xBeg, double xEnd, double eps,
    double tauMax, double tauMin, double* initApprox) {
    this->dim = dim;
    this->functions = functions;
    this->xBeg = xBeg;
    this->xEnd = xEnd;
    this->eps = eps;
    this->tauMax = tauMax;
    this->tauMin = tauMin;
    solutionsPrev = new double[dim];
    solutions = new double[dim + 1];    //Дополнительная ячейка для хранения параметра
    solutionsNext = new double[dim];
    for (int i = 0; i < dim; ++i) {
        solutionsPrev[i] = initApprox[i];
        solutions[i] = initApprox[i];
        solutionsNext[i] = initApprox[i];
    }
}

void ImplicitEulerSolver::solve() {
    double t = xBeg;
    double tNext;

    double tauPrev = tauMin;
    double tau = tauMin;
    double tauNext;

    solutions[dim] = t;

    double* epsK = new double[dim];
    vector<double> xxx(100000);
    vector<double> yy1(100000);
    vector<double> yy2(100000);
    vector<double> yy3(100000);
    int it = 0;
    Expression** adjFunctions = new Expression * [dim];   //Формирование нового вектора функций на основе исходного
    for (int i = 0; i < dim; ++i) {
        adjFunctions[i] = new DifferentialExpression(dim, functions[i], i, solutions, tau, t);
    }

    while (t < xEnd) {
        tNext = t + tau;

        for (int i = 0; i < dim; ++i) {     //Уставновка необходимых параметров для нового вектора функций
            dynamic_cast<DifferentialExpression*>(adjFunctions[i])->setTau(tau);
            dynamic_cast<DifferentialExpression*>(adjFunctions[i])->setX(t);
        }

        NewtonSolver solver(dim, adjFunctions, 100, solutions, eps);    //Решение системы нелинейных уравнений относительно solutionsNext
        solver.solve();

        for (int i = 0; i < dim; ++i) {
            solutionsNext[i] = solver.getSolutions()[i];
        }

        for (int i = 0; i < dim; ++i) {
            epsK[i] = -(tau / (tau + tauPrev)) * (solutionsNext[i] - solutions[i] - tau * (solutions[i] - solutionsPrev[i]) / tauPrev);
        }

        bool needRecount = false;   //При необходимости пересчета уменьшаем шаг
        for (int i = 0; i < dim; i++) {
            if (abs(epsK[i]) > eps and !needRecount) {
                tau /= 2;
                tNext = t;
                for (int j = 0; j < dim; j++) {
                    solutionsNext[j] = solutions[j];
                }
                needRecount = true;
            }
        }

        if (needRecount) {
            continue;
        };

        double* tauTmp = new double[dim];
        for (int i = 0; i < dim; i++) {     //Правило трех зон
            if (abs(epsK[i]) > eps) {
                tauTmp[i] = tau / 2;
            }
            else if (eps / 4 < abs(epsK[i]) && abs(epsK[i]) <= eps) {
                tauTmp[i] = tau;
            }
            else if (abs(epsK[i]) <= eps / 4) {
                tauTmp[i] = 2 * tau;
            }
            //tauTmp[i] = sqrt(eps/abs(epsK[i])) * tau; //По формуле
        }


        double curTauMin = tauTmp[0];       //Поиск минимального шага
        for (int i = 1; i < dim; i++) {
            if (curTauMin > tauTmp[i]) {
                curTauMin = tauTmp[i];
            }
        }
        delete[] tauTmp;
        tauNext = curTauMin;

        if (tauNext > tauMax) {
            tauNext = tauMax;
        }

        if ((t + tau > xEnd) and t < xEnd) {
            t = xEnd - tau;
        }

        print(t);
        xxx[it] = t;
        for (int i = 0; i < dim; i++) {
            solutionsPrev[i] = solutions[i];
            solutions[i] = solutionsNext[i];
        }

        tauPrev = tau;
        tau = tauNext;
        t = tNext;
        yy1[it] = solutions[0];
        yy2[it] = solutions[1];
        yy3[it] = solutions[2];
        it++;
    }

    delete[] epsK;
    for (int i = 0; i < dim; ++i) {
        delete adjFunctions[i];
    }


    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

    ScatterPlotSeries* series = GetDefaultScatterPlotSeriesSettings();
    series->xs = &xxx;
    series->ys = &yy1;
    //    series->xs = &xxx;
    //    series->ys = &yy3;
    series->linearInterpolation = true;
    series->lineType = toVector(L"dashed");
    series->lineThickness = 2;
    series->color = GetGray(0.3);

    ScatterPlotSettings* settings = GetDefaultScatterPlotSettings();
    settings->width = 600;
    settings->height = 400;
    settings->autoBoundaries = true;
    settings->autoPadding = true;
    settings->title = toVector(L"x^2 - 2");
    settings->xLabel = toVector(L"X axis");
    settings->yLabel = toVector(L"Y axis");
    settings->scatterPlotSeries->push_back(series);

    DrawScatterPlotFromSettings(imageReference, settings);

    vector<double>* pngdata = ConvertToPNG(imageReference->image);
    WriteToFile(pngdata, "example3.png");
    DeleteImage(imageReference->image);

    RGBABitmapImageReference* imageReference1 = CreateRGBABitmapImageReference();


    ScatterPlotSeries* series1 = GetDefaultScatterPlotSeriesSettings();
    series1->xs = &xxx;

    series1->ys = &yy2;
    series1->linearInterpolation = true;
    series1->lineType = toVector(L"dashed");
    series1->lineThickness = 2;
    series1->color = GetGray(0.3);

    ScatterPlotSettings* settings1 = GetDefaultScatterPlotSettings();
    settings1->width = 600;
    settings1->height = 400;
    settings1->autoBoundaries = true;
    settings1->autoPadding = true;
    settings1->title = toVector(L"x^2 - 2");
    settings1->xLabel = toVector(L"X axis");
    settings1->yLabel = toVector(L"Y axis");
    settings1->scatterPlotSeries->push_back(series1);

    DrawScatterPlotFromSettings(imageReference1, settings1);

    vector<double>* pngdata1 = ConvertToPNG(imageReference1->image);
    WriteToFile(pngdata1, "example31.png");
    DeleteImage(imageReference1->image);


    RGBABitmapImageReference* imageReference2 = CreateRGBABitmapImageReference();


    ScatterPlotSeries* series2 = GetDefaultScatterPlotSeriesSettings();
    series2->xs = &xxx;

    series2->ys = &yy3;
    series2->linearInterpolation = true;
    series2->lineType = toVector(L"dashed");
    series2->lineThickness = 2;
    series2->color = GetGray(0.3);

    ScatterPlotSettings* settings2 = GetDefaultScatterPlotSettings();
    settings2->width = 600;
    settings2->height = 400;
    settings2->autoBoundaries = true;
    settings2->autoPadding = true;
    settings2->title = toVector(L"x^2 - 2");
    settings2->xLabel = toVector(L"X axis");
    settings2->yLabel = toVector(L"Y axis");
    settings2->scatterPlotSeries->push_back(series2);

    DrawScatterPlotFromSettings(imageReference2, settings2);

    vector<double>* pngdata2 = ConvertToPNG(imageReference2->image);
    WriteToFile(pngdata2, "example32.png");
    DeleteImage(imageReference2->image);

    delete[] adjFunctions;
}

void ImplicitEulerSolver::print(double x) {
    cout << "X = " << x << ":   ";
    for (int i = 0; i < dim; ++i) {
        cout << "Y" << i << " = " << solutions[i] << "; ";
    }
    cout << endl;
}

ImplicitEulerSolver::~ImplicitEulerSolver() {
    delete[] solutionsPrev;
    delete[] solutions;
    delete[] solutionsNext;
}
