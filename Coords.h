
//
// Created by Robert on 9/23/2018.
//

#ifndef LENSING_COORDS_H
#define LENSING_COORDS_H


class Coords {

    Coords(double x_, double y_, double z_, double t_,
    double vx_, double vy_, double vz_, double vt_);

public:
    Coords();
    double x_, y_, z_, t_, xt_, yt_, zt_, tt_, phi;
    double vx_, vy_, vz_, vt_, vxt_, vyt_, vzt_, vtt_;

};


#endif //LENSING_COORDS_H
