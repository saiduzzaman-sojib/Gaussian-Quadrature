#include <iostream>
#include <vector>
#include <cmath>
#include <functional>
#include <fstream>
#include <string>
#include <iomanip>

using namespace std;

double transform_x(double t, double a, double b);
double scale_weight(double w, double a, double b);
double solve_gauss_2(function<double(double)> f, double a, double b);
double solve_gauss_3(function<double(double)> f, double a, double b);
double solve_gauss_4(function<double(double)> f, double a, double b);
double solve_trapezoidal(function<double(double)> f, double a, double b, int n);
double solve_simpson(function<double(double)> f, double a, double b, int n);
double target_function(double x);
double exact_solution(double a, double b);
void log_results(const string& filename, const vector<string>& methods, const vector<double>& errors);

// Start Of Implementation
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


// 3-Point Gauss-Legendre
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

// 4-Point Gauss-Legendre 
double solve_gauss_4(function<double(double)> f, double a, double b) {
    
    const double t1 = 0.339981044; 
    const double t2 = 0.861136312; 
    const double w1 = 0.652145155;
    const double w2 = 0.347854845;
    double sum = w1 * f(transform_x(-t1, a, b)) + w1 * f(transform_x(t1, a, b)) +
                 w2 * f(transform_x(-t2, a, b)) + w2 * f(transform_x(t2, a, b));

    return scale_weight(sum, a, b);
}

//Trapezoidal Rule 
double solve_trapezoidal(function<double(double)> f, double a, double b, int n) {
    double h = (b - a) / n;
    double sum = f(a) + f(b);
    for (int i = 1; i < n; i++) {
        sum += 2.0 * f(a + i * h);
    }
    return (h / 2.0) * sum;
}

//Simpson's 1/3 Rule
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


// Test Function
double target_function(double x) {
    return x * exp(x);
}
double exact_solution(double a, double b) {
    double F_b = exp(b) * (b - 1.0);
    double F_a = exp(a) * (a - 1.0);
    return F_b - F_a;
}

//Data Logging 
void log_results(const string& filename, const vector<string>& methods, const vector<double>& errors) {
    ofstream file(filename);
    
    if (file.is_open()) 
    {
        file << "Method,Error\n"; 
        for (size_t i = 0; i < methods.size(); ++i) {
            file << methods[i] << "," << errors[i] << "\n";
        }
        file.close();
        cout << "[Data Saved]: Analysis exported to " << filename << endl;
    } 
    else 
    {
        cerr << "Error: Could not open file." << endl;
    }
}

//Main Driver Function
int main() 
{
    cout << "=== Gaussian Quadrature Project (Group B3) ===" << endl;

    double a = 0.0;
    double b = 1.5;
    
    double exact = exact_solution(a, b);
    cout << "Exact Integral Value: " << exact << "\n" << endl;

    double g2 = solve_gauss_2(target_function, a, b);
    double g3 = solve_gauss_3(target_function, a, b);
    double g4 = solve_gauss_4(target_function, a, b);
    double trap = solve_trapezoidal(target_function, a, b, 4);
    double simp = solve_simpson(target_function, a, b, 4);

    vector<double> results = {g2, g3, g4, trap, simp};

    vector<double> errors;
    errors.push_back(abs(g2 - exact));
    errors.push_back(abs(g3 - exact));
    errors.push_back(abs(g4 - exact));
    errors.push_back(abs(trap - exact));
    errors.push_back(abs(simp - exact));

    vector<string> names = {"Gauss-2", "Gauss-3", "Gauss-4", "Trapezoidal", "Simpson"};
    
   
    cout << left << setw(15) << "Method" << " | " << setw(10) << "Result" << " | " << "Error" << endl;
    cout << "----------------------------------------------" << endl;
    
    for(size_t i=0; i<names.size(); i++) 
    {
        cout << left << setw(15) << names[i] << " | " <<  
                setw(10) << setprecision(6) << results[i] << " | " << 
                scientific << errors[i] << defaultfloat << endl;
    }

    log_results("results.csv", names, errors);

    return 0;
}
