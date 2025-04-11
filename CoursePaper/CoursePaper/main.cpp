#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <random>
#include <chrono>
#include "matplotlibcpp.h" // ���������� ��� ��������

namespace plt = matplotlibcpp;

class Integration {
private:
    double a, b;  // ������� ��������������
    double epsilon;  // ��������

public:
    Integration(double a, double b, double epsilon) : a(a), b(b), epsilon(epsilon) {}

    // ����� ��������
    double trapezoidal(const std::function<double(double)>& func, int n = 1000) const {
        double h = (b - a) / n;
        double result = 0.5 * (func(a) + func(b));
        for (int i = 1; i < n; ++i) {
            result += func(a + i * h);
        }
        return result * h;
    }

    // ����� ������� (��������)
    double simpson(const std::function<double(double)>& func, int n = 1000) const {
        if (n % 2 != 0) n++;  // n ������ ���� ������
        double h = (b - a) / n;
        double result = func(a) + func(b);
        for (int i = 1; i < n; i += 2) {
            result += 4 * func(a + i * h);
        }
        for (int i = 2; i < n - 1; i += 2) {
            result += 2 * func(a + i * h);
        }
        return result * h / 3.0;
    }

    // ����� �����-�����
    double monteCarlo(const std::function<double(double)>& func, int numPoints = 10000) const {
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(a, b);

        double sum = 0.0;
        for (int i = 0; i < numPoints; ++i) {
            double x = dist(gen);
            sum += func(x);
        }
        return (b - a) * sum / numPoints;
    }

    // ����� ������ (2 �����)
    double gauss(const std::function<double(double)>& func) const {
        static const double x1 = -1 / std::sqrt(3);
        static const double x2 = 1 / std::sqrt(3);
        double mid = (a + b) / 2.0;
        double halfLength = (b - a) / 2.0;
        return halfLength * (func(mid + halfLength * x1) + func(mid + halfLength * x2));
    }

    // ������������
    void visualize(const std::function<double(double)>& func, const std::string& methodName,
        const std::vector<double>& xPoints, const std::vector<double>& yPoints) const {
        plt::figure();
        plt::plot(xPoints, yPoints, "b-", { {"label", "Function"} });
        plt::fill_between(xPoints, yPoints, 0, { {"alpha", "0.3"} });
        plt::title("Integration Visualization: " + methodName);
        plt::xlabel("x");
        plt::ylabel("f(x)");
        plt::legend();
        plt::show();
    }
};

int main() {
    // ������ ������� �������������� � ��������
    Integration integrator(0, M_PI, 1e-6);

    // ������� ��� ��������������
    auto func = [](double x) { return std::sin(x); };

    // ������� ����� ��� ������������
    std::vector<double> xPoints, yPoints;
    for (double x = 0; x <= M_PI; x += 0.01) {
        xPoints.push_back(x);
        yPoints.push_back(func(x));
    }

    // ����� ��� ��������� ������� � ��������
    auto measureAndVisualize = [&](const std::string& methodName,
        const std::function<double()>& method) {
            auto start = std::chrono::high_resolution_clock::now();
            double result = method();
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            // ������ ����������
            std::cout << methodName << ": result = " << result
                << ", time = " << elapsed.count() << " seconds\n";

            // ������������
            integrator.visualize(func, methodName, xPoints, yPoints);
        };

    // ������ �������
    measureAndVisualize("Trapezoidal", [&]() { return integrator.trapezoidal(func); });
    measureAndVisualize("Simpson", [&]() { return integrator.simpson(func); });
    measureAndVisualize("Gauss", [&]() { return integrator.gauss(func); });
    measureAndVisualize("Monte Carlo", [&]() { return integrator.monteCarlo(func); });

    return 0;
}