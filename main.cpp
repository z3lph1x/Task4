#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

double g(double x, double y, double u) {
    return 1 / (x * x) + exp(1) / (log(x)) * y * y - exp(u) * y;
}

vector<double> Runge_Kutta4(double a, double b, double y0, double u0, int N) {
    double h = (b - a) / N;
    vector<double> x;
    x.resize(N);
    x[0] = a;
    for (int i = 1; i < N; i++)
        x[i] = x[i - 1] + h;
    vector<double> y;
    y.resize(N);
    vector<double> u;
    u.resize(N);
    y[0] = y0;
    u[0] = u0;
    double k1, k2, k3, k4, l1, l2, l3, l4;
    for (int i = 1; i < N; i++) {
        k1 = h * u[i - 1];
        l1 = h * g(x[i - 1], y[i - 1], u[i - 1]);
        k2 = h * (u[i - 1] + l1 / 2);
        l2 = h * g(x[i - 1] + h / 2, y[i - 1] + k1 / 2, u[i - 1] + l1 / 2);
        k3 = h * (u[i - 1] + l2 / 2);
        l3 = h * g(x[i - 1] + h / 2, y[i - 1] + k2 / 2, u[i - 1] + l2 / 2);
        k4 = h * (u[i - 1] + l3);
        l4 = h * g(x[i - 1] + h / 2, y[i - 1] + k3 / 2, u[i - 1] + l3 / 2);
        y[i] = y[i - 1] + 1. / 6 * (k1 + 2 * (k2 + k3) + k4);
        u[i] = u[i - 1] + 1. / 6 * (l1 + 2 * (l2 + l3) + l4);
    }
    return y;
}

int main() {
    int N = 10;
    double e = 0.001;
    double a = exp(1);
    double b = exp(2);
    double h = (b - a) / N;
    double alpha1 = 5.;
    double ans1 = Runge_Kutta4(a, b, exp(1), alpha1 + h, N)[N - 1];//alpha1 + h
    double ans2 = Runge_Kutta4(a, b, exp(1), alpha1, N)[N - 1];//alpha
    double df = (ans1 - ans2) / h;
    double alpha2 = alpha1 - (ans2 - exp(2)) / df;
    alpha1 = alpha2;
    while (abs(ans2- 2*exp(1))> e){
        double ans1 = Runge_Kutta4(a, b, exp(1), alpha1 + h, N)[N - 1];//alpha1 + h
        double ans2 = Runge_Kutta4(a, b, exp(1), alpha1, N)[N - 1];//alpha
        df = (ans1 - ans2) / h;
        alpha2 = alpha1 - (ans2 - 2 * exp(2)) / df;
    }
    cout << alpha2 << "<- alpha; " << ans2 << " <- definition of function at the last iteration";
}
