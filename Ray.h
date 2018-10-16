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
    void step(Ray ray);   /* I was not sure here if step would need to have both
                            rays in it, I don't think so, because the trajectory
                           of one does not depend on the trajectory of the other. */

    double M_ = 1.0;
    Ray();
    double phifunc(Coords coords);
    double goo(Coords coords);
    double g11(Coords coords);
    double g22(Coords coords);
    double g33(Coords coords);
    double vxdotfunc(Coords coords); // so I set them to double just incase
    double vydotfunc(Coords coords);
    double vzdotfunc(Coords coords);
    double vtdotfunc(Coords coords);
    double lag(Coords coords);
    double position(); /* Could make this append some values for positions to arrays
                       for each coordinate, or could just do that in step */
    double gooderiv(Coords coords, int dir);
};


#endif //LENSING_RAY_H
