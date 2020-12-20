#include <iostream>
#include <random>
#include <ctime>
#include <vector>
#include <algorithm>

const double a = 0.0;
const double b = 3.0;
const double c = -1.0;
const double d = 3.0;
const size_t N = 10;
const double A = 3.0;
const double step = (double)((b - a) / (N - 1));

struct Point {
    double x;
    double y;
};

double func (const double &x) {
    return c * x + d;
}

std::vector<Point> random(const size_t N, const double noise) {
    std::vector<Point> points (N);
    std::random_device rd;
    std::mt19937 gen(rd());
    std::uniform_real_distribution<double> error(-0.5, 0.5);
    for (size_t i = 0; i < N; ++i) {
        points[i].x = a + i * step;
        points[i].y = func( a + i * step) + noise * error(gen);
    }
    return points;
}

std::vector<Point> edge(const double lower, const double upper,
                        const size_t num, const double noise) {
    std::vector<Point> points(num);
    const double step = (upper - lower) / static_cast<double>(num - 1);
    for (size_t i = 0; i < num; ++i) {
        points[i].x = lower + i * step;
        points[i].y = func(points[i].x) -noise/2 + rand() * 1./RAND_MAX * (noise);
    }
    return points;
}

double error(const std::vector<Point>& right, const double c, const double d) {
    double sum = 0.;
    for (auto point : right) {
        sum += pow(point.y - (c * func(point.x)), 2);
    }
    return sum;
}

double golden_ratio(std::vector<Point>& p, double Cmin, double Cmax) {
    double l = std::abs(Cmax - Cmin);
    std::swap(Cmin, Cmax);
    Cmin = std::fabs(Cmin);
    Cmax = std::fabs(Cmax);
    const double e = 0.1;
    const double t = (std::sqrt(5) + 1) / 2;
    double c_k1 = Cmin + (1 - 1/t)*Cmax;
    double c_k2 = Cmin + Cmax / t;
    double f_k1 = func(-c_k1);
    double f_k2 = func(-c_k2);
    while (l > e){
        if (f_k1 < f_k2){
            Cmax = c_k2;
            c_k2 = Cmin + Cmax - c_k1;
            f_k2 = func(-c_k2);
        } else {
            Cmin = c_k1;
            c_k1 = Cmin + Cmax - c_k2;
            f_k1 = func(-c_k1);
        }
        if (c_k1 > c_k2){
            std::swap(c_k1, c_k2);
            std::swap(f_k1, f_k2);
        }
        l = std::abs(Cmax - Cmin);
    }
    return -((Cmax + Cmin) / 2);
}

int F(int f)
{
    if (f==1) return 1;
    else if (f==2) return 1;
    else if (f>2)
        return (F(f-1)+F(f-2));
}

double Fibonacci( std::vector<Point>& p, double Dmin, double Dmax) {
    double ak = Dmin, bk = Dmax, x1, x2, y1, y2;
    int G = 10;
    x1 = ak + (double)F(G - 2) / F(G) * (bk - ak);
    x2 = ak + (double)F(G - 1) / F(G) * (bk - ak);
        y1 = error(p, x1, 0);
        y2 = error(p, x2, 0);
    for (int i=G; i >= 1; --i) {
        if (y1 > y2) {
            ak = x1;
            x1 = x2;
            x2 = bk - (x2 - ak);
            y1 = y2;
            y2 = error(p, x2, 0);
        }
        else {
            bk = x2;
            x2 = x1;
            x1 = ak + (bk - x2);
            y2 = y1;
            y1 = error(p, x1, 0);
        }
    }
    return (x1 + x2) / 2;
}


void print(const double noise)
{
    std::vector<Point> p = random(N, noise);
    double Cmin, Cmax, Dmin, Dmax;
    edge( Cmin, Cmax, Dmin, Dmax);
    std::cout << "Cmin = " << Cmin << "\nCmax = " << Cmax << "\nDmin = "
              << Dmin << "\nDmax = " << Dmax << std::endl;
    double w1 = golden_ratio(p, Cmin, Cmax);
    double w0 = Fibonacci(p, Dmin, Dmax);
    std::cout << "w1 = " << w1 <<std::endl;
    std::cout<< "w0 = " << w0 << std::endl;
}

int main() {
    std::cout << "Function y = -x + 3"<<std::endl;
    std::cout << "Without noise"<<std::endl;
    print(0.0);
    std::cout << "With noise A = " << A <<std::endl;
    print(A);
    return 0;
}
