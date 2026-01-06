#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

// Helper Functions
double transform_x(double t, double a, double b);
double scale_weight(double w, double a, double b);
// Gaussian Algorithms
double solve_gauss_2(function<double(double)> f, double a, double b);
double solve_gauss_3(function<double(double)> f, double a, double b);
double solve_gauss_4(function<double(double)> f, double a, double b);
// Comparison Algorithms
double solve_trapezoidal(function<double(double)> f, double a, double b, int n);
double solve_simpson(function<double(double)> f, double a, double b, int n);
// Math & Data
double target_function(double x);
double exact_solution(double a, double b);
void log_results(const string& filename, const vector<string>& methods, const vector<double>& errors);


// START OF IMPLEMENTATION

double transform_x(double t, double a, double b) {
    return ((b - a) * t + (b + a)) / 2.0;
}
double scale_weight(double sum, double a, double b) {
    return ((b - a) / 2.0) * sum;
}
