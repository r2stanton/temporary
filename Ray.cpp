//
// Created by Robert on 9/23/2018.
//

#include "Ray.h"
#include "Coords.h"

Ray::Ray() {

}

void::step(Ray ray){
    double h = .001;

    double  k1x, k1y, k1z, k1t, k1vx, k1vy, k1vz, k1vt,
            k2x, k2y, k2z, k2t, k2vx, k2vy, k2vz, k2vt,
            k3x, k3y, k3z, k3t, k3vx, k3vy, k3vz, k3vt,
            k4x, k4y, k4z, k4t, k4vx, k4vy, k4vz, k4vt;

    ray.coords.xt_ = ray.coords.x_;
    ray.coords.yt_ = ray.coords.y_;
    ray.coords.zt_ = ray.coords.z_;
    ray.coords.tt_ = ray.coords.t_;
    ray.coords.vxt_ = ray.coords.vx_;
    ray.coords.vyt_ = ray.coords.vy_;
    ray.coords.vzt_ = ray.coords.vz_;
    ray.coords.vtt_ = ray.coords.vt_;

    k1x = h * ray.coords.vxt_;
    k1y = h * ray.coords.vyt_;
    k1z = h * ray.coords.vzt_;
    k1t = h * ray.coords.vtt_;

    /*
    k1vx = h * ?
    k1vy = h * ?
    k1vz = h * ?
    k1vt = h * ?

    update temps with k
    don't understand what xdot and vxdot functions end up doing

     */

}

// I'm still not sure I get what these functions are going to look like, or
// how the A/B's will differ, if we need an At, Bt, Avt, Bvt, etc. and are they constant?

double::vtdotfunc(/* don't think we need to pass anything here? */){
    // Somethinge like 1/g00t * (A-B);
}

double::vxdotfunc(){

}

double::vydotfunc(){

}

double::vzdotfunc(){

}

double::tdotfunc(){

}

double::xdotfunc(){

}

double::ydotfunc(){

}

double::zdotfunc(){

}
