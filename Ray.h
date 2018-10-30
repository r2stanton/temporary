//
// Created by Robert on 9/23/2018.
//

#ifndef LENSING_RAY_H
#define LENSING_RAY_H

#include "Coords.h"
double M;


class Ray {
public:
    Coords coords_;
    // Maybe Coords temp; ? Depends on RK setup
    void step(Ray ray);
    double M_ = 1.0;
    Ray();
    double phifunc(Coords temp);
    double goo(Coords temp);
    double g11(Coords temp);
    double g22(Coords temp);
    double g33(Coords temp);
    double vxdotfunc(Coords temp);
    double vydotfunc(Coords temp);
    double vzdotfunc(Coords temp);
    double vtdotfunc(Coords temp);
    double lag(Coords coords);
    double position();
    double gooderiv(Coords temp, int dir);
    double g11deriv(Coords temp, int dir);
    double g22deriv(Coords temp, int dir);
    double g33deriv(Coords temp, int dir);
};


#endif //LENSING_RAY_H
