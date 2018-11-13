#include <iostream>
#include "Ray.h"

using namespace std;

int main() {

    cout<< "Getting started" << endl;
    double thetasec = 35.0;
    double vx = thetasec*4.84814e-6;

    Ray ray1(vx);
    ray1.runtoz(1e+7);

    cout << ray1.coords_.x_ << endl;

 //   Ray ray2(-0.3);
 //   ray2.runtoz(20.0);

//    cout << ray2.coords_.x_ << endl;

}