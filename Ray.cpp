//
// Created by Robert on 9/23/2018.
//

#include <cmath>
#include <iostream>
#include "Ray.h"
#include "Coords.h"

using namespace std;

Ray::Ray() {

    coords_.x_ = 10.0;
    coords_.y_ = 0.0;
    coords_.z_ = 0.0;
    coords_.t_ = 0.0;
    coords_.vx_ = 0.1;
    coords_.vy_ = 0.0;
    coords_.vz_ = sqrt(1.0 - coords_.vx_*coords_.vx_ - coords_.vy_*coords_.vy_);

}

void Ray::step(double h){
    double  k1x, k1y, k1z, k1t, k1vx, k1vy, k1vz, k1vt,
            k2x, k2y, k2z, k2t, k2vx, k2vy, k2vz, k2vt,
            k3x, k3y, k3z, k3t, k3vx, k3vy, k3vz, k3vt,
            k4x, k4y, k4z, k4t, k4vx, k4vy, k4vz, k4vt;

    coords_.tt_ = coords_.t_;
    coords_.xt_ = coords_.x_;
    coords_.yt_ = coords_.y_;
    coords_.zt_ = coords_.z_;
    coords_.vtt_ = coords_.vt_;
    coords_.vxt_ = coords_.vx_;
    coords_.vyt_ = coords_.vy_;
    coords_.vzt_ = coords_.vz_;

    k1t = h * coords_.vtt_;
    k1x = h * coords_.vxt_;
    k1y = h * coords_.vyt_;
    k1z = h * coords_.vzt_;


    k1vt = h * vtdotfunc();
    k1vx = h * vxdotfunc();
    k1vy = h * vydotfunc();
    k1vz = h * vzdotfunc();

    // Here I know some stuff has to be reset. Vdotfuncs depend on g__ , g__,i , v_ which need to be updated with temps.
    // Trying to use temp values to compute g__  , g__,i , v_ seems impossible how it is currently set up
    // because they are hardcoded to deal with ray.coords.vx_, x_, etc. even phifunc does the same thing.
    // Need to figure easiest way for




}

// I'm still not sure I get what these functions are going to look like, or
// how the A/B's will differ, if we need an At, Bt, Avt, Bvt, etc. and are they constant?

double  Ray::phifunc(Coords temp){
    return (-1* M)/(pow(pow(temp.xt_, 2) + pow(temp.yt_, 2) + pow(temp.zt_, 2),.5));
}

double Ray::goo(Coords temp){
    return(-1 -2*phifunc(temp));
}
double Ray::g11(Coords temp){
    return(1-2*phifunc(temp));
}
double Ray::g22(Coords temp){
    return(-1 -2*phifunc(temp));
}
double Ray::g33(Coords temp){
    return(-1 -2*phifunc(temp));
}



// Biggest part I'm struggling with is how to do derivatives for the g__'s
// Don't know how I'm going to be able to pass ray.coords.x_ + h, might have to create temporary Ray or something

double Ray::vtdotfunc(Coords temp){   //right below is dL/dt

    double dgoodt = gooderiv(temp, 0);
    return 1.0/goo(temp) * (.5*(dgoodt * pow(temp.vt_, 2_) + dg11dt * pow(temp.vx_, 2) + dg22dt * pow(temp.vy_, 2)
            + dg33dt * pow(temp.vz_, 2)))

                        - 1/ray.goo(ray)*(  dgoodt*ray.coords.vt_*ray.coords.vt_
                                            + dgoodx*ray.coords.vt_*ray.coords.vx_
                                            + dgoody*ray.coords.vt_*ray.coords.vy_
                                            + dgoodz*ray.coords.vt_*ray.coords.vz_);
}

double Ray::gooderiv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = goo(temp);
        temp.tt_ = val - delta;
        minus = goo(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = goo(temp);
        temp.xt_ = val - delta;
        minus = goo(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = goo(temp);
        temp.yt_ = val - delta;
        minus = goo(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = goo(temp);
        temp.zt_ = val - delta;
        minus = goo(temp);
    }
    else{
        cout << "bad integer in gooderiv" << endl;
        break;
    }
    return (plus - minus)/2.0/delta;
}

double::vxdotfunc(Ray ray){
    1/ray.g11(ray) * (.5*(dgoodx * pow(ray.coords.vt_, 2_) + dg11dx * pow(ray.coords.vx_, 2) + dg22dx * pow(ray.coords.vy_, 2) + dg33dx * pow(ray.coords.vz_, 2)))
                            - 1/ray.g11(ray)*(  dgoodt*ray.coords.vx_*ray.coords.vt_
                                              + dgoodx*ray.coords.vx_*ray.coords.vx_
                                              + dgoody*ray.coords.vx_*ray.coords.vy_
                                              + dgoodz*ray.coords.vx_*ray.coords.vz_);
}

double::vydotfunc(Ray ray){
    1/ray.g22(ray) * (.5*(dgoody * pow(ray.coords.vt_, 2_) + dg11dy * pow(ray.coords.vx_, 2) + dg22dy * pow(ray.coords.vy_, 2) + dg33dy * pow(ray.coords.vz_, 2)))
                            - 1/ray.g22(ray)*(  dgoodt*ray.coords.vy_*ray.coords.vt_
                                              + dgoodx*ray.coords.vy_*ray.coords.vx_
                                              + dgoody*ray.coords.vy_*ray.coords.vy_
                                              + dgoodz*ray.coords.vy_*ray.coords.vz_);
}

double::vzdotfunc(Ray ray){
    1/ray.g33(ray) * (.5*(dgoodz * pow(ray.coords.vt_, 2_) + dg11dz * pow(ray.coords.vx_, 2) + dg22dz * pow(ray.coords.vy_, 2) + dg33dz * pow(ray.coords.vz_, 2)))
                            - 1/ray.g33(ray)*(  dgoodt*ray.coords.vz_*ray.coords.vt_
                                              + dgoodx*ray.coords.vz_*ray.coords.vx_
                                              + dgoody*ray.coords.vz_*ray.coords.vy_
                                              + dgoodz*ray.coords.vz_*ray.coords.vz_);
}
