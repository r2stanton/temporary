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

// For RK if use temps
            // Can run 2 rays at once, Coords temp1 Coords temp2 = Ray1.coords. ... ,  Ray2.coords. ...
            // Bad thing is it would look very strange (Ray.step(Coords temp1, Coords temp2) referring to 2 different rays, being called by one ray
            // and might have problems because step is in Ray, would likely want to remove step from Ray entirely

// For RK if use coords, it would make most sense to *only use them for the initial/final setting of Coord temp object,
            // So: temp. ... = coords. ... at beginning -> RK -> coords. ... = temp. ...
            // This would probably work well, and we'd get rid of coords_.(t, x, y, z)t's since temp object would handle that
// I think this is better, but I could be missing something with what you did with temps

// Basically if we don't get rid of one or the other, there will be a lot of unnecessary setting of either
            // coords_.(t, x, y, z)t = temp. ... or vice versa

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





}



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



double Ray::vtdotfunc(Coords temp){   //right below is dL/dt
    return (1.0/goo(temp) *(.5*(gooderiv(temp, 0) * pow(temp.vt_, 2)    +
                                g11deriv(temp, 0) * pow(temp.vx_, 2)    +
                                g22deriv(temp, 0) * pow(temp.vy_, 2)    +
                                g33deriv(temp, 0) * pow(temp.vz_, 2)))  -
            1/goo(temp)   *    (gooderiv(temp, 0) * temp.vt_ * temp.vt_ +
                                gooderiv(temp, 1) * temp.vt_ * temp.vx_ +
                                gooderiv(temp, 2) * temp.vt_ * temp.vy_ +
                                gooderiv(temp, 3) * temp.vt_ * temp.vz_));
}

double Ray::vxdotfunc(Coords temp){   //right below is dL/dt
    return (1.0/g11(temp) *(.5*(gooderiv(temp, 1) * pow(temp.vt_, 2)    +
                                g11deriv(temp, 1) * pow(temp.vx_, 2)    +
                                g22deriv(temp, 1) * pow(temp.vy_, 2)    +
                                g33deriv(temp, 1) * pow(temp.vz_, 2)))  -
            1/g11(temp)   *    (g11deriv(temp, 0) * temp.vx_ * temp.vt_ +
                                g11deriv(temp, 1) * temp.vx_ * temp.vx_ +
                                g11deriv(temp, 2) * temp.vx_ * temp.vy_ +
                                g11deriv(temp, 3) * temp.vx_ * temp.vz_));
}

double Ray::vtdotfunc(Coords temp){   //right below is dL/dt
    return (1.0/g22(temp) *(.5*(gooderiv(temp, 2) * pow(temp.vt_, 2)    +
                                g11deriv(temp, 2) * pow(temp.vx_, 2)    +
                                g22deriv(temp, 2) * pow(temp.vy_, 2)    +
                                g33deriv(temp, 2) * pow(temp.vz_, 2)))  -
            1/g22(temp)   *    (g22deriv(temp, 0) * temp.vy_ * temp.vt_ +
                                g22deriv(temp, 1) * temp.vy_ * temp.vx_ +
                                g22deriv(temp, 2) * temp.vy_ * temp.vy_ +
                                g22deriv(temp, 3) * temp.vy_ * temp.vz_));
}

double Ray::vtdotfunc(Coords temp){   //right below is dL/dt
    return (1.0/g33(temp) *(.5*(gooderiv(temp, 3) * pow(temp.vt_, 2)    +
                                g11deriv(temp, 3) * pow(temp.vx_, 2)    +
                                g22deriv(temp, 3) * pow(temp.vy_, 2)    +
                                g33deriv(temp, 3) * pow(temp.vz_, 2)))  -
            1/g33(temp)   *    (g33deriv(temp, 0) * temp.vz_ * temp.vt_ +
                                g33deriv(temp, 1) * temp.vz_ * temp.vx_ +
                                g33deriv(temp, 2) * temp.vz_ * temp.vy_ +
                                g33deriv(temp, 3) * temp.vz_ * temp.vz_));
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

double Ray::g11deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g11(temp);
        temp.tt_ = val - delta;
        minus = g11(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g11(temp);
        temp.xt_ = val - delta;
        minus = g11(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g11(temp);
        temp.yt_ = val - delta;
        minus = g11(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g11(temp);
        temp.zt_ = val - delta;
        minus = g11(temp);
    }
    else{
        cout << "bad integer in g11deriv" << endl;
        break;
    }
    return (plus - minus)/2.0/delta;
}

double Ray::g22deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g22(temp);
        temp.tt_ = val - delta;
        minus = g22(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g22(temp);
        temp.xt_ = val - delta;
        minus = g22(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g22(temp);
        temp.yt_ = val - delta;
        minus = g22(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g22(temp);
        temp.zt_ = val - delta;
        minus = g22(temp);
    }
    else{
        cout << "bad integer in g22deriv" << endl;
        break;
    }
    return (plus - minus)/2.0/delta;
}

double Ray::g33deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g33(temp);
        temp.tt_ = val - delta;
        minus = g33(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g33(temp);
        temp.xt_ = val - delta;
        minus = g33(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g33(temp);
        temp.yt_ = val - delta;
        minus = g33(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g33(temp);
        temp.zt_ = val - delta;
        minus = g33(temp);
    }
    else{
        cout << "bad integer in g33deriv" << endl;
        break;
    }
    return (plus - minus)/2.0/delta;
}

