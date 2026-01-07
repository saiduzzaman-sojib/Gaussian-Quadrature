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

double solve_gauss_2(function<double(double)> f, double a, double b) {
    
    const double t = 0.577350269; 
    const double w = 1.0;

    double sum = w * f(transform_x(-t, a, b)) + 
                 w * f(transform_x(t, a, b));

    return scale_weight(sum, a, b);
}
// 3-Point Gauss-Legendre ---

double solve_gauss_3(function<double(double)> f, double a, double b) {
    
    const double t1 = 0.0; 
    const double t2 = 0.774596669; 
    
    
    const double w1 = 0.888888889; 
    const double w2 = 0.555555556; 

    double sum = w1 * f(transform_x(t1, a, b)) +
                 w2 * f(transform_x(-t2, a, b)) +
                 w2 * f(transform_x(t2, a, b));

    return scale_weight(sum, a, b);
}

// --- MEMBER 5: 4-Point Gauss-Legendre ---

double solve_gauss_4(function<double(double)> f, double a, double b) {
    
    const double t1 = 0.339981044; 
    const double t2 = 0.861136312; 

    const double w1 = 0.652145155;
    const double w2 = 0.347854845;

    double sum = w1 * f(transform_x(-t1, a, b)) + w1 * f(transform_x(t1, a, b)) +
                 w2 * f(transform_x(-t2, a, b)) + w2 * f(transform_x(t2, a, b));

    return scale_weight(sum, a, b);
}
// ---Trapezoidal Rule ---
// Used for comparison with Gaussian methods
double solve_trapezoidal(function<double(double)> f, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        sum += 2.0 * f(a + i * h);
    }
    return (h / 2.0) * sum;
}
// --- Simpson's 1/3 Rule ---
// Higher accuracy comparison

double solve_simpson(function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) n++; 
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        if (i % 2 == 0)
            sum += 2.0 * f(a + i * h);
        else
            sum += 4.0 * f(a + i * h);
    }

    return (h / 3.0) * sum;
}
// ---Simpson's 1/3 Rule ---
// Higher accuracy comparison

double solve_simpson(function<double(double)> f, double a, double b, int n) {
    if (n % 2 != 0) n++; 
    
    double h = (b - a) / n;
    double sum = f(a) + f(b);

    for (int i = 1; i < n; i++) {
        if (i % 2 == 0)
            sum += 2.0 * f(a + i * h);
        else
            sum += 4.0 * f(a + i * h);
    }

    return (h / 3.0) * sum;
}
