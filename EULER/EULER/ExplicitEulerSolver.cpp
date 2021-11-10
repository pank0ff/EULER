#include "ExplicitEulerSolver.h"
#include "vector"
#include "pbPlots.hpp"
#include "supportLib.hpp"
ExplicitEulerSolver::ExplicitEulerSolver(int dim, Expression** functions, double xBeg, double xEnd, double eps,
    double tauMax, double* initApprox) {
    this->dim = dim;
    this->functions = functions;
    this->xBeg = xBeg;
    this->xEnd = xEnd;
    this->eps = eps;
    this->tauMax = tauMax;
    solutions = new double[dim + 1];    //Дополнительная ячейка для параметра t
    for (int i = 0; i < dim; ++i) {
        solutions[i] = initApprox[i];
    }
}

void ExplicitEulerSolver::print(double x) {
    cout << "X = " << x << ":   ";
    for (int i = 0; i < dim; ++i) {
        cout << "Y" << i << " = " << solutions[i] << "; ";
    }
    cout << endl;
}

vector<double> xxx(100000);
vector<double> yy1(100000);
vector<double> yy2(100000);
vector<double> yy3(100000);
int it = 0;

void ExplicitEulerSolver::solve() {
    double* fSols = new double[dim];
    double xCur = xBeg;
    while (xCur <= xEnd) {
        print(xCur);
        xxx[it] = xCur;
        solutions[dim] = xCur;                      //Нахождение F(Yk, tk)
        for (int i = 0; i < dim; ++i) {
            fSols[i] = functions[i]->F(solutions);
        }

        double tau = eps / (abs(fSols[0]) + eps / tauMax);      //Вычисление шага
        for (int i = 0; i < dim; ++i) {
            if (eps / (abs(fSols[i]) + eps / tauMax) < tau) {
                tau = eps / (abs(fSols[i]) + eps / tauMax);
            }
        }

        for (int i = 0; i < dim; ++i) {
            solutions[i] = solutions[i] + tau * fSols[i];
        }
        yy1[it] = solutions[0];
        yy2[it] = solutions[1];
        //        yy3[it] = solutions[3];


        xCur += tau;
        it++;
    }
    RGBABitmapImageReference* imageReference = CreateRGBABitmapImageReference();

    ScatterPlotSeries* series = GetDefaultScatterPlotSeriesSettings();
    series->xs = &xxx;
    series->ys = &yy1;
    //    series->ys = &yy2;
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
    WriteToFile(pngdata, "example2.png");
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
    WriteToFile(pngdata1, "example21.png");
    DeleteImage(imageReference1->image);

    delete[] fSols;
}

ExplicitEulerSolver::~ExplicitEulerSolver() {
    delete[] solutions;
}
