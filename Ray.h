//
// Created by Robert on 9/23/2018.
//

#ifndef LENSING_RAY_H
#define LENSING_RAY_H

#include "Coords.h"
double g00x, g00y, g00z, g00t,
       g11x, g11y, g11z, g11t,
       g22x, g22y, g22z, g22t,
       g33x, g33y, g33z, g33t;

double A, B; // I don't understand the math here or what is going on
             // so I don't get if each A, B is different
             // or if we need for both vx, vydots etc, and x, ydots
             // or just vxdots and then get x, ydots from there


class Ray {
public:
    Coords coords;
    void step(Ray ray);   /* I was not sure here if step would need to have both
                            rays in it, I don't think so, because the trajectory
                           of one does not depend on the trajectory of the other. */

    double xdotfunc(Ray ray);


    Ray();

    double ydotfunc();  // These could all either be void and change variables
    double zdotfunc();  // or they could be this and return values that just
    double tdotfunc();  // we set Ray.Coords.x_, y_, etc = to
    double vxdotfunc(); // so I set them to double just incase
    double vydotfunc();
    double vzdotfunc();
    double vtdotfunc();

    double position(); /* Could make this append some values for positions to arrays
                       for each coordinate, or could just do that in step */
};


#endif //LENSING_RAY_H
