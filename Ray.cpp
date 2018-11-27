//
// Created by Robert on 9/23/2018.
//

#include <cmath>
#include <iostream>
#include <fstream>
#include "Ray.h"
#include "Coords.h"
#include "derivs.cpp"

using namespace std;


Ray::Ray() {

    // use to set the initial values
    coords_.x_ = 0.0;
    coords_.y_ = 0.0;
    coords_.z_ = -10.0;
    coords_.t_ = 0.0;
    coords_.vx_ = 0.4;
    coords_.vy_ = 0.0;
    coords_.vz_ = sqrt(1.0 - coords_.vx_*coords_.vx_ - coords_.vy_*coords_.vy_);
    coords_.vt_ = 1.0;
}


Ray::Ray(double vxval, double vyval, double zval) {

    // use to set the initial values
    coords_.x_ = 0.0;
    coords_.y_ = 0.0;
    coords_.z_ = zval;  // seems right for 4 billion years and 10^12 solar masses
    coords_.t_ = 0.0;
    coords_.vx_ = vxval;
    coords_.vy_ = vyval;
    coords_.vz_ = sqrt(1.0 - coords_.vx_*coords_.vx_ - coords_.vy_*coords_.vy_);
    coords_.vt_ = 1.0;
}


Ray::~Ray() {
    // basic empty destructor
}

void Ray::runray() {
    //Commented lines will write to csv. Might be easier to compare the rays w/ a quicker language like python
    //ofstream data;
    //data.open("data.csv");
    double h = 0.005;
    for(int i=0; i<3000; i++){
        step(h);
        if(i%100 == 0){
            cout << i <<"  "<< coords_.t_ <<"  " << coords_.x_ <<"  " << coords_.y_ <<"  "<< coords_.z_ <<"  " << endl;
        }
        //data << i << "," << coords_.t_ << "," << coords_.x_ << "," << coords_.y_ <<","<< coords_.z_<< endl;
    }
}

void Ray::runtoz(double zfin, double xaim, double yaim) {
    double h = 0.005;
    int runsteps = 0;
    ofstream fout;
    fout.open("DVxVyfor3DnewBy5.csv", ios::app);
    double vxi = coords_.vx_;
    double vyi = coords_.vy_;
    while(coords_.z_ < zfin){
        if(abs(coords_.z_)<1000.0){
            h = .1;
        }
        else if(abs(coords_.z_) < 10000){
            h = 300.0;
        }
        else {
            h = 2000.0;
        }
        step(h);
        runsteps++;
        ///Used to check if values are reasonable
        if(runsteps%1000 == 0){
            // cout << runsteps <<"  "<< coords_.t_ <<"  " << coords_.x_ <<"  " << coords_.y_ <<"  "<< coords_.z_ <<"  " << coords_.vz_ << endl;
            // for when vz is acting strangely cout << coords_.vx_ << "  "  <<coords_.vy_ << "  " << coords_.vz_ << "  " << sqrt(1.0 - coords_.vx_*coords_.vx_ - coords_.vy_*coords_.vy_) << endl;
        }


        /* Use for csv to view trajectory in x vs z
        if(runsteps%10000 ==0){
            fout<< coords_.x_ <<"," << coords_.z_<<endl;
        }
        */

        if(runsteps > 1000000){
            cout<<"took too many steps"<<endl;
            break;
        }

        if(pow(coords_.x_, 2) + pow(coords_.y_, 2) + pow(coords_.z_, 2) < 4*M_){
            cout<<"got too close"<<endl;
            break;
        }
    }
    fout << pow( (pow( (coords_.x_ - xaim) , 2) + pow( (coords_.y_ - yaim) , 2)), .5) << "," << vxi << "," << vyi << endl;
    fout.close();
    cout << runsteps <<"  "<< coords_.t_ <<"  " << coords_.x_ <<"  " << coords_.y_ <<"  "<< coords_.z_ <<"  " << endl;
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


    coords_.tt_ = coords_.t_   + k1t/2.0;
    coords_.xt_ = coords_.x_   + k1x/2.0;
    coords_.yt_ = coords_.y_   + k1y/2.0;
    coords_.zt_ = coords_.z_   + k1z/2.0;
    coords_.vtt_ = coords_.vt_ + k1vt/2.0;
    coords_.vxt_ = coords_.vx_ + k1vx/2.0;
    coords_.vyt_ = coords_.vy_ + k1vy/2.0;
    coords_.vzt_ = coords_.vz_ + k1vz/2.0;

    k2t = h * coords_.vtt_;
    k2x = h * coords_.vxt_;
    k2y = h * coords_.vyt_;
    k2z = h * coords_.vzt_;
    k2vt = h * vtdotfunc();
    k2vx = h * vxdotfunc();
    k2vy = h * vydotfunc();
    k2vz = h * vzdotfunc();


    coords_.tt_ = coords_.t_   + k2t/2.0;
    coords_.xt_ = coords_.x_   + k2x/2.0;
    coords_.yt_ = coords_.y_   + k2y/2.0;
    coords_.zt_ = coords_.z_   + k2z/2.0;
    coords_.vtt_ = coords_.vt_ + k2vt/2.0;
    coords_.vxt_ = coords_.vx_ + k2vx/2.0;
    coords_.vyt_ = coords_.vy_ + k2vy/2.0;
    coords_.vzt_ = coords_.vz_ + k2vz/2.0;

    k3t = h * coords_.vtt_;
    k3x = h * coords_.vxt_;
    k3y = h * coords_.vyt_;
    k3z = h * coords_.vzt_;
    k3vt = h * vtdotfunc();
    k3vx = h * vxdotfunc();
    k3vy = h * vydotfunc();
    k3vz = h * vzdotfunc();

    coords_.tt_ = coords_.t_   + k3t;
    coords_.xt_ = coords_.x_   + k3x;
    coords_.yt_ = coords_.y_   + k3y;
    coords_.zt_ = coords_.z_   + k3z;
    coords_.vtt_ = coords_.vt_ + k3vt;
    coords_.vxt_ = coords_.vx_ + k3vx;
    coords_.vyt_ = coords_.vy_ + k3vy;
    coords_.vzt_ = coords_.vz_ + k3vz;

    k4t = h * coords_.vtt_;
    k4x = h * coords_.vxt_;
    k4y = h * coords_.vyt_;
    k4z = h * coords_.vzt_;
    k4vt = h * vtdotfunc();
    k4vx = h * vxdotfunc();
    k4vy = h * vydotfunc();
    k4vz = h * vzdotfunc();

    coords_.t_ = coords_.t_ + (k1t + 2.0*k2t + 2.0*k3t + k4t)/6.0;
    coords_.x_ = coords_.x_ + (k1x + 2.0*k2x + 2.0*k3x + k4x)/6.0;
    coords_.y_ = coords_.y_ + (k1y + 2.0*k2y + 2.0*k3y + k4y)/6.0;
    coords_.z_ = coords_.z_ + (k1z + 2.0*k2z + 2.0*k3z + k4z)/6.0;
    coords_.vt_ = coords_.vt_ + (k1vt + 2.0*k2vt + 2.0*k3vt + k4vt)/6.0;
    coords_.vx_ = coords_.vx_ + (k1vx + 2.0*k2vx + 2.0*k3vx + k4vx)/6.0;
    coords_.vy_ = coords_.vy_ + (k1vy + 2.0*k2vy + 2.0*k3vy + k4vy)/6.0;
    coords_.vz_ = coords_.vz_ + (k1vz + 2.0*k2vz + 2.0*k3vz + k4vz)/6.0;

}



double  Ray::phifunc(Coords temp){
    return (- M_)/(pow(pow(temp.xt_, 2) + pow(temp.yt_, 2) + pow(temp.zt_, 2),.5));
}

double Ray::goo(Coords temp){
    return(-1 -2*phifunc(temp));
}
double Ray::g11(Coords temp){
    return(1-2*phifunc(temp)); // I'm almost sure that the bottom ones were supposed to have a positive 1, so I changed to that.
}
double Ray::g22(Coords temp){
    return(1 -2*phifunc(temp));
}
double Ray::g33(Coords temp){
    return(1 -2*phifunc(temp));
}

double Ray::vtdotfunc(){

    Coords temp;
    temp.tt_ = coords_.tt_;
    temp.xt_ = coords_.xt_;
    temp.yt_ = coords_.yt_;
    temp.zt_ = coords_.zt_;

    double dgoodt = gooderiv(temp, 0);
    double dg11dt = g11deriv(temp, 0);
    double dg22dt = g22deriv(temp, 0);
    double dg33dt = g33deriv(temp, 0);
    double dgoodx = gooderiv(temp, 1);
    double dgoody = gooderiv(temp, 2);
    double dgoodz = gooderiv(temp, 3);


    return 1.0/goo(temp) * (.5*(dgoodt * pow(coords_.vtt_, 2.0) + dg11dt * pow(coords_.vxt_, 2) + dg22dt * pow(coords_.vyt_, 2)
            + dg33dt * pow(coords_.vzt_, 2)))  - 1/goo(temp)*(  dgoodt*coords_.vtt_*coords_.vtt_
            + dgoodx*coords_.vtt_*coords_.vxt_ + dgoody*coords_.vtt_*coords_.vyt_
            + dgoodz*coords_.vtt_*coords_.vzt_);
}

double Ray::vxdotfunc(){

    Coords temp;
    temp.tt_ = coords_.tt_;
    temp.xt_ = coords_.xt_;
    temp.yt_ = coords_.yt_;
    temp.zt_ = coords_.zt_;


    double dgoodx = gooderiv(temp, 1);
    double dg11dx = g11deriv(temp, 1);
    double dg22dx = g22deriv(temp, 1);
    double dg33dx = g33deriv(temp, 1);
    double dg11dt = g11deriv(temp, 0);
    double dg11dy = g11deriv(temp, 2);
    double dg11dz = g11deriv(temp, 3);


    return 1.0/g11(temp) * (.5*(dgoodx * pow(coords_.vtt_, 2.0) + dg11dx * pow(coords_.vxt_, 2) + dg22dx * pow(coords_.vyt_, 2)
            + dg33dx * pow(coords_.vzt_, 2)))
            - 1/g11(temp)*(  dg11dt*coords_.vxt_*coords_.vtt_  + dg11dx*coords_.vxt_*coords_.vxt_
                    + dg11dy*coords_.vxt_*coords_.vyt_  + dg11dz*coords_.vxt_*coords_.vzt_);


}

double Ray::vydotfunc(){

    Coords temp;
    temp.tt_ = coords_.tt_;
    temp.xt_ = coords_.xt_;
    temp.yt_ = coords_.yt_;
    temp.zt_ = coords_.zt_;


    double dgoody = gooderiv(temp, 2);
    double dg11dy = g11deriv(temp, 2);
    double dg22dy = g22deriv(temp, 2);
    double dg33dy = g33deriv(temp, 2);
    double dg22dt = g22deriv(temp, 0);
    double dg22dx = g22deriv(temp, 1);
    double dg22dz = g22deriv(temp, 3);


    return 1.0/g22(temp) * (.5*(dgoody * pow(coords_.vtt_, 2.0) + dg11dy * pow(coords_.vxt_, 2) + dg22dy * pow(coords_.vyt_, 2)
                                + dg33dy * pow(coords_.vzt_, 2)))
           - 1/g22(temp)*(  dg22dt*coords_.vyt_*coords_.vtt_  + dg22dx*coords_.vyt_*coords_.vxt_
                            + dg22dy*coords_.vyt_*coords_.vyt_  + dg22dz*coords_.vyt_*coords_.vzt_);


}

double Ray::vzdotfunc(){

    Coords temp;
    temp.tt_ = coords_.tt_;
    temp.xt_ = coords_.xt_;
    temp.yt_ = coords_.yt_;
    temp.zt_ = coords_.zt_;


    double dgoodz = gooderiv(temp, 3);
    double dg11dz = g11deriv(temp, 3);
    double dg22dz = g22deriv(temp, 3);
    double dg33dz = g33deriv(temp, 3);
    double dg33dt = g33deriv(temp, 0);
    double dg33dx = g33deriv(temp, 1);
    double dg33dy = g33deriv(temp, 2);


    return 1.0/g33(temp) * (.5*(dgoodz * pow(coords_.vtt_, 2.0) + dg11dz * pow(coords_.vxt_, 2) + dg22dz * pow(coords_.vyt_, 2)
                                + dg33dz * pow(coords_.vzt_, 2)))
           - 1/g33(temp)*(  dg33dt*coords_.vzt_*coords_.vtt_  + dg33dx*coords_.vzt_*coords_.vxt_
                            + dg33dy*coords_.vzt_*coords_.vyt_  + dg33dz*coords_.vzt_*coords_.vzt_);


}


double Ray::gooderiv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = goo(temp);
        temp.tt_ = val - delta;
        minus = goo(temp);

    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = goo(temp);
        temp.xt_ = val - delta;
        minus = goo(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = goo(temp);
        temp.yt_ = val - delta;
        minus = goo(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = goo(temp);
        temp.zt_ = val - delta;
        minus = goo(temp);
    }
    else{
        cout << "bad integer in gooderiv" << endl;
    }
    return (plus - minus)/2.0/delta;
}

double Ray::g11deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g11(temp);
        temp.tt_ = val - delta;
        minus = g11(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g11(temp);
        temp.xt_ = val - delta;
        minus = g11(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g11(temp);
        temp.yt_ = val - delta;
        minus = g11(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g11(temp);
        temp.zt_ = val - delta;
        minus = g11(temp);
    }
    else{
        cout << "bad integer in gooderiv" << endl;
    }
    return (plus - minus)/2.0/delta;
}

double Ray::g22deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g22(temp);
        temp.tt_ = val - delta;
        minus = g22(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g22(temp);
        temp.xt_ = val - delta;
        minus = g22(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g22(temp);
        temp.yt_ = val - delta;
        minus = g22(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g22(temp);
        temp.zt_ = val - delta;
        minus = g22(temp);
    }
    else{
        cout << "bad integer in gooderiv" << endl;
    }
    return (plus - minus)/2.0/delta;
}

double Ray::g33deriv(Coords temp, int dir) {

    double delta, plus, minus, val;
    if(dir ==0){
        delta = temp.tt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.tt_;
        temp.tt_ = val + delta;
        plus = g33(temp);
        temp.tt_ = val - delta;
        minus = g33(temp);
    }
    else if(dir == 1){
        delta = temp.xt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.xt_;
        temp.xt_ = val + delta;
        plus = g33(temp);
        temp.xt_ = val - delta;
        minus = g33(temp);
    }

    else if(dir == 2){
        delta = temp.yt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.yt_;
        temp.yt_ = val + delta;
        plus = g33(temp);
        temp.yt_ = val - delta;
        minus = g33(temp);
    }

    else if(dir == 3){
        delta = temp.zt_ / 100.0;
        if (delta == 0){
            delta = defdelta;
        }
        val = temp.zt_;
        temp.zt_ = val + delta;
        plus = g33(temp);
        temp.zt_ = val - delta;
        minus = g33(temp);
    }
    else{
        cout << "bad integer in gooderiv" << endl;
    }
    return (plus - minus)/2.0/delta;
}