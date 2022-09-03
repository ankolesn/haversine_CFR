#include <iostream>
#include <cmath>
#include <ctime>
#include <chrono>
#include <fstream>

using namespace std;

double haversine(double lat1, double lon1, double lat2, double lon2) {
    const double earthRadiusKm = 6371.0;
    double u = sin((lat2 - lat1) / 2);
    double v = sin((lon2 - lon1) / 2);

    return 2.0 * earthRadiusKm * asin(sqrt(u * u + cos(lat1) * cos(lat2) * v * v));
}

double CFR(double lat1, double lon1, double lat2, double lon2) {

    double ml = (lat1 + lat2) / 2.0;
    double kpd_lat = 111.13209 - 0.56605 * cos(2 * ml) + 0.00120 * cos(4 * ml);
    double kpd_lon = 111.41513 * cos(ml) - 0.09455 * cos(3 * ml) + 0.00012 * cos(5 * ml);
    double ns = kpd_lat * (lat1 - lat2);
    double ew = kpd_lon * (lon1 - lon2);

    return sqrt(ns * ns + ew * ew);
}

double CFR_optimized(double lat1, double lon1, double lat2, double lon2) {

    double ml = (lat1 + lat2) / 2.0;

    double cos1ml = cos(ml);
    double sin1ml = sin(ml);

    double cos2ml = cos1ml * cos1ml - sin1ml * sin1ml;
    double sin2ml = 2.0 * sin1ml * cos1ml;
    double cos3ml = cos1ml * (2.0 * cos2ml - 1.0);
    double cos4ml = cos2ml * cos2ml - sin2ml * sin2ml;
    double cos5ml = cos1ml * (-2.0 * cos2ml + 2.0 * cos4ml + 1.0);

    double kpd_lat = 111.13209 - 0.56605 * cos2ml + 0.00120 * cos4ml;
    double kpd_lon = 111.41513 * cos1ml - 0.09455 * cos3ml + 0.00012 * cos5ml;
    double ns = kpd_lat * (lat1 - lat2);
    double ew = kpd_lon * (lon1 - lon2);

    return sqrt(ns * ns + ew * ew);
}

double speed(double lat1, double lon1, double lat2, double lon2,
             double (*ptr)(double lat1, double lon1, double lat2, double lon2)) {

    auto start = std::chrono::system_clock::now();
    double func = (*ptr)(lat1, lon1, lat2, lon2);
    auto end = std::chrono::system_clock::now();

    std::chrono::duration<double> diff = end - start;
    auto res = duration_cast<std::chrono::nanoseconds>(diff);

    return res.count();
}

double accuracy(double lat1, double lon1, double lat2, double lon2,
                double (*ptr)(double lat1, double lon1, double lat2, double lon2)) {
    lat1 = lat1 * M_PI / 180;
    lon1 = lon1 * M_PI / 180;
    lat2 = lat2 * M_PI / 180;
    lon2 = lon2 * M_PI / 180;

    double angular_len = acos(sin(lat1) * sin(lat2) + cos(lat1) * cos(lat2) * cos(lon2 - lon1));
    double len = 111.1 * angular_len;

    if (isnan(len)){
        return 0;
    }

    double res = (*ptr)(lat1, lon1, lat2, lon2);
    return abs(len - res);
}


int main() {
    double (*ptr)(double lat1, double lon1, double lat2, double lon2) = nullptr;
    int size = 100000;

    string filename1 = "speed.dat";
    string filename2 = "accuracy.dat";

    ofstream out1(filename1);
    ofstream out2(filename2);

    double sum_speed = 0;
    double sum_accur = 0;
    srand(time(nullptr));


    for (int i = 0; i < size; ++i) {
        ptr = haversine;
        double num1 = double (1 + rand() % 100) / 1000;
        double num2 = double (1 + rand() % 100) / 1000;
        double num3 = double (1 + rand() % 100) / 1000;
        double num4 = double (1 + rand() % 100) / 1000;

        sum_speed += speed(num1, num2, num3, num4, ptr);
        sum_accur += accuracy(num1, num2, num3, num4, ptr);
    }
    int n = 0;
    int k = 0;

    out1 << n++ << " haversine " << sum_speed / 100000.0 << endl;
    out2 << k++ << " haversine " << sum_accur / 100000.0 << endl;

    cout << "Haversine speed: " << sum_speed / 100000.0 << endl;
    cout << "Haversine accuracy: " << sum_accur / 100000.0 << endl;

    ptr = CFR;
    sum_speed = 0;
    sum_accur = 0;

    for (int i = 0; i < size; ++i) {
        double num1 = double (1 + rand() % 100) / 1000;
        double num2 = double (1 + rand() % 100) / 1000;
        double num3 = double (1 + rand() % 100) / 1000;
        double num4 = double (1 + rand() % 100) / 1000;

        sum_speed += speed(num1, num2, num3, num4, ptr);
        sum_accur += accuracy(num1, num2, num3, num4, ptr);
    }

    out1 << n++ << " CFR " << sum_speed / 100000.0 << endl;
    out2 << k++ << " CFR " << sum_accur / 100000.0 << endl;

    cout << "CFR speed: " << sum_speed / 100000.0 << endl;
    cout << "CFR accuracy: " << sum_accur / 100000.0 << endl;

    ptr = CFR_optimized;
    sum_speed = 0;
    sum_accur = 0;

    for (int i = 0; i < size; ++i) {
        double num1 = double (1 + rand() % 100) / 1000;
        double num2 = double (1 + rand() % 100) / 1000;
        double num3 = double (1 + rand() % 100) / 1000;
        double num4 = double (1 + rand() % 100) / 1000;

        sum_speed += speed(num1, num2, num3, num4, ptr);
        sum_accur += accuracy(num1, num2, num3, num4, ptr);
    }

    out1 << n++ << " CFR_optimized " << sum_speed / 100000.0 << endl;
    out2 << k++ << " CFR_optimized " << sum_accur / 100000.0 << endl;

    cout << "CFR optimized speed: " << sum_speed / 100000.0 << endl;
    cout << "CFR accuracy: " << sum_accur / 100000.0 << endl;

    return 0;
}
