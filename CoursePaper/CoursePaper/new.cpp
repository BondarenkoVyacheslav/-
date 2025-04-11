#include <iostream>
#include <cmath>
#include <functional>
#include <vector>
#include <random>
#include <chrono>
#include "matplotlibcpp.h" // Библиотека для графиков

namespace plt = matplotlibcpp;

// Класс для хранения пределов интегрирования
class IntegrationLimits {
public:
    double a, b;

    IntegrationLimits(double a, double b) : a(a), b(b) {}

    void setLimits(double newA, double newB) {
        a = newA;
        b = newB;
    }
};

// Базовый класс для всех методов интегрирования
class Integration {
protected:
    IntegrationLimits limits;
    double epsilon;

public:
    Integration(IntegrationLimits limits, double epsilon) : limits(limits), epsilon(epsilon) {}

    virtual double integrate(const std::function<double(double)>& func) const = 0; // Чисто виртуальный метод
};

// Класс трапеций
class TrapezoidalIntegration : public Integration {
public:
    TrapezoidalIntegration(IntegrationLimits limits, double epsilon) : Integration(limits, epsilon) {}

    double integrate(const std::function<double(double)>& func) const override {
        int n = 1000;
        double h = (limits.b - limits.a) / n;
        double result = 0.5 * (func(limits.a) + func(limits.b));
        for (int i = 1; i < n; ++i) {
            result += func(limits.a + i * h);
        }
        return result * h;
    }
};

// Класс метода Симпсона
class SimpsonIntegration : public Integration {
public:
    SimpsonIntegration(IntegrationLimits limits, double epsilon) : Integration(limits, epsilon) {}

    double integrate(const std::function<double(double)>& func) const override {
        int n = 1000;
        if (n % 2 != 0) n++;
        double h = (limits.b - limits.a) / n;
        double result = func(limits.a) + func(limits.b);
        for (int i = 1; i < n; i += 2) {
            result += 4 * func(limits.a + i * h);
        }
        for (int i = 2; i < n - 1; i += 2) {
            result += 2 * func(limits.a + i * h);
        }
        return result * h / 3.0;
    }
};

// Класс Монте-Карло
class MonteCarloIntegration : public Integration {
public:
    MonteCarloIntegration(IntegrationLimits limits, double epsilon) : Integration(limits, epsilon) {}

    double integrate(const std::function<double(double)>& func) const override {
        int numPoints = 10000;
        std::random_device rd;
        std::mt19937 gen(rd());
        std::uniform_real_distribution<> dist(limits.a, limits.b);

        double sum = 0.0;
        for (int i = 0; i < numPoints; ++i) {
            double x = dist(gen);
            sum += func(x);
        }
        return (limits.b - limits.a) * sum / numPoints;
    }
};

// Класс метода Гаусса
class GaussIntegration : public Integration {
public:
    GaussIntegration(IntegrationLimits limits, double epsilon) : Integration(limits, epsilon) {}

    double integrate(const std::function<double(double)>& func) const override {
        static const double x1 = -1 / std::sqrt(3);
        static const double x2 = 1 / std::sqrt(3);
        double mid = (limits.a + limits.b) / 2.0;
        double halfLength = (limits.b - limits.a) / 2.0;
        return halfLength * (func(mid + halfLength * x1) + func(mid + halfLength * x2));
    }
};

// Класс визуализации
class Visualization {
public:
    static void visualize(const std::function<double(double)>& func, const std::string& methodName,
        const std::vector<double>& xPoints, const std::vector<double>& yPoints) {
        plt::figure();
        plt::plot(xPoints, yPoints, "b-");

        for (size_t i = 0; i < xPoints.size(); ++i) {
            std::vector<double> xArea = { xPoints[i], xPoints[i] };
            std::vector<double> yArea = { 0, yPoints[i] };
            plt::plot(xArea, yArea, "b-");
        }

        plt::title("Integration Visualization: " + methodName);
        plt::xlabel("x");
        plt::ylabel("f(x)");
        plt::show();
    }
};

int main() {
    const int M_PI = 10;
    // Устанавливаем пределы интегрирования
    IntegrationLimits limits(0, M_PI);
    double epsilon = 1e-6;

    // Создаем методы интегрирования
    TrapezoidalIntegration trapezoidalIntegrator(limits, epsilon);
    SimpsonIntegration simpsonIntegrator(limits, epsilon);
    MonteCarloIntegration monteCarloIntegrator(limits, epsilon);
    GaussIntegration gaussIntegrator(limits, epsilon);

    // Функции для интегрирования
    auto func1 = [](double x) { return std::sin(x); };
    auto func2 = [](double x) { return x * x; };

    // Создаем точки для визуализации
    auto createPoints = [&](const std::function<double(double)>& func) {
        std::vector<double> xPoints, yPoints;
        for (double x = limits.a; x <= limits.b; x += 0.01) {
            xPoints.push_back(x);
            yPoints.push_back(func(x));
        }
        return std::make_pair(xPoints, yPoints);
        };

    auto measureAndVisualize = [&](Integration& integrator, const std::string& methodName,
        const std::function<double(double)>& func) {

            auto Points = createPoints(func);
            std::vector<double> xPoints = Points.first;
            std::vector<double> yPoints = Points.second;

            auto start = std::chrono::high_resolution_clock::now();
            double result = integrator.integrate(func);
            auto end = std::chrono::high_resolution_clock::now();
            std::chrono::duration<double> elapsed = end - start;

            // Печать результата
            std::cout << methodName << ": result = " << result
                << ", time = " << elapsed.count() << " seconds\n";

            // Визуализация
            Visualization::visualize(func, methodName, xPoints, yPoints);
        };

    // Пример работы с sin(x)
    std::cout << "Integrating sin(x):\n";
    measureAndVisualize(trapezoidalIntegrator, "Trapezoidal", func1);
    measureAndVisualize(simpsonIntegrator, "Simpson", func1);
    measureAndVisualize(gaussIntegrator, "Gauss", func1);
    measureAndVisualize(monteCarloIntegrator, "Monte Carlo", func1);

    // Пример работы с x^2
    limits.setLimits(0, 1);
    std::cout << "\nIntegrating x^2:\n";
    measureAndVisualize(trapezoidalIntegrator, "Trapezoidal", func2);
    measureAndVisualize(simpsonIntegrator, "Simpson", func2);
    measureAndVisualize(monteCarloIntegrator, "Monte Carlo", func2);

    return 0;
}