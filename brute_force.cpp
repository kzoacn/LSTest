#include <iostream>
#include <vector>
#include <algorithm>
#include <cassert>
#include "rain.hpp"
#include <ctime>

using namespace std;




int main(){

    Rain rain;
    GF secret_key = 123;
    GF x = 123;
    GF y = rain.enc_2r(secret_key,x);

    double st = clock();
    for(int k=0;k<1000;k++){
        GF key = k;
        GF c = rain.enc_2r(key,x);
        if(c==y){
            cout << "k = " << k << endl; 
        }
    }

    //cout << "time = " << (clock()-st)/CLOCKS_PER_SEC << endl;
    cout<<"average time = "<<(clock()-st)/CLOCKS_PER_SEC<<"ms"<<endl;

    return 0;
}