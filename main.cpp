#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>

using namespace std;

double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double earthRadiusKm = 6371.0;
    double u = sin((lat2 - lat1) / 2);
    double v = sin((lon2 - lon1) / 2);

    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1) * cos(lat2) * v * v));
}

double CFR(double lat1, double lon1, double lat2, double lon2) {
    double ml = (lat1 + lat2) / 2.0 * M_PI / 180.0;
    double kpd_lat = 111132.09 - 566.05 * cos(2 * ml) + 1.20 * cos(4 * ml);
    double kpd_lon = 111415.13 * cos(ml) - 94.55 * cos(3 * ml) + 0.12 * cos(5 * ml);
    double ns = kpd_lat * (lat1 - lat2);
    double ew = kpd_lon * (lon1 - lon2);

    return sqrt(ns * ns + ew * ew);
}

double CFR_optimized(double lat1, double lon1, double lat2, double lon2) {
    double ml = (lat1 + lat2) / 2.0 * M_PI / 180.0;

    double cos1ml = cos(ml);
    double sin1ml = sin(ml);

    double cos2ml = cos1ml * cos1ml - sin1ml * sin1ml;
    double sin2ml = 2.0 * sin1ml * cos1ml;
    double cos3ml = cos1ml * (2.0 * cos2ml - 1.0);
    double cos4ml = cos2ml * cos2ml - sin2ml * sin2ml;
    double cos5ml = cos1ml * (-2.0 * cos2ml + 2.0 * cos4ml + 1.0);

    double kpd_lat = 111132.09 - 566.05 * cos2ml + 1.20 * cos4ml;
    double kpd_lon = 111415.13 * cos1ml - 94.55 * cos3ml + 0.12 * cos5ml;
    double ns = kpd_lat * (lat1 - lat2);
    double ew = kpd_lon * (lon1 - lon2);

    return sqrt(ns * ns + ew * ew);
}

double speed(double lat1, double lon1, double lat2, double lon2, double (*ptr)(double lat1, double lon1, double lat2, double lon2)) {

    auto start = std::chrono::system_clock::now();
    double func = (*ptr)(lat1, lon1, lat2, lon2);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> diff = end - start;
    
    return diff.count();
}

double accuracy(double lat1, double lon1, double lat2, double lon2, double (*ptr)(double lat1, double lon1, double lat2, double lon2)){
    //...
}


int main() {
    double (*ptr)(double lat1, double lon1, double lat2, double lon2) = nullptr;
    ptr = haversine;

    double sum = 0;
    srand(time(0));


    for (int i = 0; i < 100000; ++i) {
        double num1 = 1 + rand() % 100;
        double num2 = 1 + rand() % 100;
        double num3 = 1 + rand() % 100;
        double num4 = 1 + rand() % 100;

        sum += (*ptr)(num1, num2, num3, num4);
    }

    cout << "Haversine speed: " << sum / 100000.0 << endl;

    ptr = CFR;
    sum = 0;

    for (int i = 0; i < 100000; ++i) {
        double num1 = 1 + rand() % 100;
        double num2 = 1 + rand() % 100;
        double num3 = 1 + rand() % 100;
        double num4 = 1 + rand() % 100;

        sum += (*ptr)(num1, num2, num3, num4);
    }

    cout << "CFR speed: " << sum / 100000.0 << endl;

    ptr = CFR_optimized;
    sum = 0;

    for (int i = 0; i < 10000; ++i) {
        double num1 = 1 + rand() % 100;
        double num2 = 1 + rand() % 100;
        double num3 = 1 + rand() % 100;
        double num4 = 1 + rand() % 100;

        sum += (*ptr)(num1, num2, num3, num4);
    }

    cout << "CFR optimized speed: " << sum / 10000.0 << endl;
    return 0;
}
