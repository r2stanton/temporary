#include <iostream>
#include "Ray.h"

using namespace std;

int main() {

    cout << "Getting started" << endl;
    double thetasecx = 10.0;
    double thetasecy = 10;
    double vy = thetasecy * 4.8414e-6;
    double vx = thetasecx * 4.84814e-6;
    double zval = -1e+6;
    for(int i =0;i<100;i+=5){
        double vx = (thetasecx + i) * 4.84814e-6;
        for(int j = 0;j<100;j+=5){
            double vy = (thetasecy + j) * 4.8414e-6;
            Ray ray1(vx, vy, zval);
            ray1.runtoz(1e+6, 10000, 5000);
            cout << "i = " << i << "|  j  = " << j << endl;
        }
    }
            //Ray ray1(vx);
            //ray1.runtoz(1e+7);
    //cout << ray1.coords_.x_ << endl;

 //   Ray ray2(-0.3);
 //   ray2.runtoz(20.0);

//    cout << ray2.coords_.x_ << endl;

}