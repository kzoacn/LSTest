#ifndef RAIN_HPP_
#define RAIN_HPP_

#include <random>
#include <bitset>
#include <vector>
#include "constant.h"
#include "bitbasis.hpp"

using std::vector;
using std::bitset;



typedef bitset<N> GF;

//X^N + X^7 + X^2 + X + 1
GF setup(){
    GF mod;
#if N!=256
    mod[7]=1;
    mod[2]=1;
    mod[1]=1;
    mod[0]=1;
#else
    //X256 + X10 + X5 + X2 + X + 1.
    mod[10]=1;
    mod[5]=1;
    mod[2]=1; 
    mod[0]=1;
#endif
    return mod;
}

const GF GF_mod=setup();
std::mt19937 mt(123);

GF GF_rand(){
    GF res;
    for(int i=0;i<N;i++){
        res[i]=mt()%2;
    }
    return res;
}

GF GF_mul(GF x,GF y){
    GF res = 0;
    for(int i=0;i<N;i++){
        if(y[i]==1)
            res=res^x;
        if (x[N-1]==1)
            x=(x<<1)^GF_mod;
        else
            x=x<<1;
    }
    return res;
}
GF GF_pow(GF x,unsigned long long y){
    GF res = 1;
    //quick
    while(y>0){
        if(y&1)
            res=GF_mul(res,x);
        x=GF_mul(x,x);
        y>>=1;
    }
    return res;
}

GF GF_inv(GF x){
    GF res = 1;
    //2^N-2 = 11111..11110
    for(int i=0;i<N;i++){
        if(i!=0)
            res=GF_mul(res,x);
        x=GF_mul(x,x);
    }
    return res;
}

GF mat_mut(const vector<GF> &M,GF x){
    GF res = 0;
    for(int i=0;i<N;i++){
        res[i] = (M[i]&x).count()%2;
    }
    return res;
}
vector<GF> mat_mut(const vector<GF> &A,const vector<GF> &B){
    vector<GF> C;
    C.resize(N);

    for(int i=0;i<N;i++){
        for(int j=0;j<N;j++){
            C[i][j]=0;
            for(int k=0;k<N;k++){
                if(A[i][k]&B[k][j])
                    C[i][j].flip();
            }
        }
    }
    return C;
}  

vector<GF> mat_inv(vector<GF> A){
    vector<GF> B;
    B.resize(N);
    for(int i=0;i<N;i++){
        B[i][i]=1;
    }

    for(int i=0;i<N;i++){
        if(A[i][i]==0){
            for(int j=i+1;j<N;j++){
                if(A[j][i]==1){
                    swap(A[i],A[j]);
                    swap(B[i],B[j]);
                    break;
                }
            }
        }
        for(int j=0;j<N;j++){
            if(i!=j && A[j][i]==1){
                A[j]=A[j]^A[i];
                B[j]=B[j]^B[i];
            }
        }
    }
    return B;
}
 


vector<GF> rand_mat(){
    int r=0;
    vector<bitset<N> >vec;
    BitBasis<N> basis;
    std::mt19937 mt;
    while(r<N){
        bitset<N>bs;
        for(int i=0;i<N;i++)
            bs[i]=mt()%2;
        if(basis.insert(bs)){
            vec.push_back(bs);
            r++;
        }
    }
    return vec;
}

struct Rain{
    std::mt19937 mt;
    GF c1;
    GF c2;
    vector<GF> M1;
    vector<GF> M2;

    Rain(){
        mt.seed(123);
        c1=GF_rand();
        c2=GF_rand();
        M1=rand_mat();
        M2=rand_mat();
    }

    void simple(){
        c1=1;
        c2=0;
        M1=rand_mat();
        for(int i=0;i<N;i++){
            M2[i]=0; 
            M2[i][i]=1;
        }

    }

    GF enc_1r(GF k,GF P){
        GF x = P;
        x=x^k^c1;
        x=GF_inv(x);
        x=mat_mut(M1,x);
        x=x^k;
        return x;
    }

    GF enc_2r(GF k,GF P){
        GF x = P;
        x=x^k^c1;
        x=GF_inv(x);
        x=mat_mut(M1,x);
        x=x^k^c2;
        x=GF_inv(x);
        x=mat_mut(M2,x);
        x=x^k;
        return x;
    }
};

#endif
