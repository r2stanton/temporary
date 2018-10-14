//
// Created by Robert on 9/23/2018.
//

#include <cmath>
#include "Ray.h"
#include "Coords.h"

Ray::Ray() {

}

void::step(Ray ray){
    double h = .001;
    double M;
    double  k1x, k1y, k1z, k1t, k1vx, k1vy, k1vz, k1vt,
            k2x, k2y, k2z, k2t, k2vx, k2vy, k2vz, k2vt,
            k3x, k3y, k3z, k3t, k3vx, k3vy, k3vz, k3vt,
            k4x, k4y, k4z, k4t, k4vx, k4vy, k4vz, k4vt;

    ray.coords.tt_ = ray.coords.t_;
    ray.coords.xt_ = ray.coords.x_;
    ray.coords.yt_ = ray.coords.y_;
    ray.coords.zt_ = ray.coords.z_;
    ray.coords.vtt_ = ray.coords.vt_;
    ray.coords.vxt_ = ray.coords.vx_;
    ray.coords.vyt_ = ray.coords.vy_;
    ray.coords.vzt_ = ray.coords.vz_;

    k1t = h * ray.coords.vtt_;
    k1x = h * ray.coords.vxt_;
    k1y = h * ray.coords.vyt_;
    k1z = h * ray.coords.vzt_;


    k1vt = h * ray.vtdotfunc(ray);
    k1vx = h * ray.vxdotfunc(ray);
    k1vy = h * ray.vydotfunc(ray);
    k1vz = h * ray.vzdotfunc(ray);

    // Here I know some stuff has to be reset. Vdotfuncs depend on g__ , g__,i , v_ which need to be updated with temps.
    // Trying to use temp values to compute g__  , g__,i , v_ seems impossible how it is currently set up
    // because they are hardcoded to deal with ray.coords.vx_, x_, etc. even phifunc does the same thing.
    // Need to figure easiest way for




}

// I'm still not sure I get what these functions are going to look like, or
// how the A/B's will differ, if we need an At, Bt, Avt, Bvt, etc. and are they constant?

double::phifunc(Ray ray){
    return (-1* M)/(pow(pow(ray.coords.x_, 2) + pow(ray.coords.y_, 2) + pow(ray.coords.z_, 2),.5));
}

double::goo(Ray ray){
    return(-1 -2*ray.phifunc(ray));
}
double::g11(Ray ray){
    return(1-2*ray.phifunc(ray));
}
double::g22(Ray ray){
    return(-1 -2*ray.phifunc(ray));
}
double::g33(Ray ray){
    return(-1 -2*ray.phifunc(ray));
}



// Biggest part I'm struggling with is how to do derivatives for the g__'s
// Don't know how I'm going to be able to pass ray.coords.x_ + h, might have to create temporary Ray or something

double::vtdotfunc(Ray ray){   //right below is dL/dt
    1/ray.goo(ray) * (.5*(dgoodt * pow(ray.coords.vt_, 2_) + dg11dt * pow(ray.coords.vx_, 2) + dg22dt * pow(ray.coords.vy_, 2) + dg33dt * pow(ray.coords.vz_, 2)))
                        - 1/ray.goo(ray)*(  dgoodt*ray.coords.vt_*ray.coords.vt_
                                            + dgoodx*ray.coords.vt_*ray.coords.vx_
                                            + dgoody*ray.coords.vt_*ray.coords.vy_
                                            + dgoodz*ray.coords.vt_*ray.coords.vz_);
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
