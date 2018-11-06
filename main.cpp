#include <iostream>
#include "Ray.h"

using namespace std;

int main() {

    cout<< "Getting started" << endl;
    Ray ray1(0.4);
    ray1.runtoz(20.0);

    cout << ray1.coords_.x_ << endl;

    Ray ray2(-0.3);
    ray2.runtoz(20.0);

    cout << ray2.coords_.x_ << endl;

}