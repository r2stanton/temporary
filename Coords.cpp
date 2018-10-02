//
// Created by Robert on 9/23/2018.
//

#include "Coords.h"

Coords::Coords(double xi, double yi, double zi, double ti,
               double vxi, double vyi, double vzi, double vti) {

    x_ = xi;
    y_ = yi;
    z_ = zi;
    t_ = ti;
    vt_ = vti;
    vx_ = vxi;
    vy_ = vyi;
    vz_ = vzi;
}
