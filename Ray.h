//
// Created by Robert on 9/23/2018.
//

#ifndef LENSING_RAY_H
#define LENSING_RAY_H

#include "Coords.h"
double M;


class Ray {
public:
    Coords coords;
    void step(Ray ray);   /* I was not sure here if step would need to have both
                            rays in it, I don't think so, because the trajectory
                           of one does not depend on the trajectory of the other. */


    Ray();
    double phifunc(Ray ray);
    double goo(Ray ray);
    double g11(Ray ray);
    double g22(Ray ray);
    double g33(Ray ray);
    double vxdotfunc(Ray ray); // so I set them to double just incase
    double vydotfunc(Ray ray);
    double vzdotfunc(Ray ray);
    double vtdotfunc(Ray ray);
    double lag(Ray ray);
    double position(); /* Could make this append some values for positions to arrays
                       for each coordinate, or could just do that in step */
};


#endif //LENSING_RAY_H
