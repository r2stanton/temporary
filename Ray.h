//
// Created by Robert on 9/23/2018.
//

#ifndef LENSING_RAY_H
#define LENSING_RAY_H

#include "Coords.h"
#include "derivs.cpp"


class Ray {
public:
    Coords coords_;

    void step(double h);
    void runray();
    void runtoz(double zfin, double xaim, double yaim);

    double defdelta = .00001;
    double M_ = 1.0;  // could replace with a true value 10^30 kg range
    Ray();
    Ray(double vx, double vy, double zval);
    ~Ray();

    double phifunc(Coords temp);
    double goo(Coords temp);
    double g11(Coords temp);
    double g22(Coords temp);
    double g33(Coords temp);
    double vxdotfunc();
    double vydotfunc();
    double vzdotfunc();
    double vtdotfunc();
    double lag(Coords coords);
    double position();
    double gooderiv(Coords temp, int dir);
    double g11deriv(Coords temp, int dir);
    double g22deriv(Coords temp, int dir);
    double g33deriv(Coords temp, int dir);

};


#endif //LENSING_RAY_H
