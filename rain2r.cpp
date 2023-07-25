#include <cmath>
#include<iostream>
#include<algorithm>
#include<vector>
#include<set>
#include<map> 
#include "bitbasis.hpp" 
#include <random>
#include "rain.hpp"
using namespace std;

#define RAIN 
//#define AIM

#ifdef RAIN

    #define N 128
    #define LOGD 8

    #define T (2*LOGD)
    #define FREE_VARS (2*LOGD)
    #define QUAD_VARS (FREE_VARS+FREE_VARS*(FREE_VARS-1)/2)
#endif 

const int n=N; 
vector<int> irr_poly;

Rain rain;
GF P,K,C,X=1;
vector<GF> M1,M1_inv;

void set_parameters(){
    irr_poly.resize(n+1); 
 
    irr_poly[n]=1;
    irr_poly[7]=1;
    irr_poly[2]=1;
    irr_poly[1]=1;
    irr_poly[0]=1;

}

vector<int> free_vars;
int mapping[N][N];
pair<int,int> rmapping[QUAD_VARS];

struct Term{
    int var1,var2;
    Term(int a){
        var1=var2=a;
    }
    Term(int a,int b){
        var1=min(a,b);
        var2=max(a,b);
    }
    bool operator==(const Term& t)const{
        return var1==t.var1 && var2==t.var2;
    }
    bool operator<(const Term& t)const{
        if(var1!=t.var1)return var1<t.var1;
        return var2<t.var2;
    }
    Term operator*(const Term& oth)const{

        if(var1==var2){
            if(oth.var1==oth.var2)
                return Term(var1,oth.var1);
            else
                return Term(var1,oth.var1^oth.var2^var1);
        }else{
            return Term(var1,var2);
        } 

    }
};

string to_str(int x){
    string s;
    if(x==0)return "0";
    while(x){
        s+=char('0'+x%10);
        x/=10;
    }
    reverse(s.begin(),s.end());
    return s;
}

string name(int id){
    if(id<n)
        return string("x")+to_str(id);
    return string("y")+to_str(id-n);
}



struct Expression{
    //vector<Term> terms;
    //int constant;

    bitset<QUAD_VARS+1>terms;
    //terms[QUAD_VARS] is constant

    // exp= \sum terms + constant
    Expression(){ } 
    Expression(int t){
        terms[t]=1;
    }
    void set_const(int c){
        terms[QUAD_VARS]=c;
    }
    void add_one(){
        terms[QUAD_VARS].flip();
    }
    int get_const()const{
        return terms[QUAD_VARS];
    }

    void print(bool quad_flag=false){

        vector<string>text;
        if(quad_flag){
            for(int i=0;i<QUAD_VARS;i++)if(terms[i]){
                int v1=rmapping[i].first,v2=rmapping[i].second;
                if(v1==v2){
                    text.push_back(name(v1));
                }else{
                    text.push_back(name(v1)+name(v2));
                }
            }

            return ;
        }else{
            for(int i=0;i<128;i++)if(terms[i]){
                text.push_back(name(i));
            }
        }

        sort(text.begin(),text.end());
        
        if(terms[QUAD_VARS]){
            text.push_back("1");
        }
        for(int i=0;i<text.size();i++){
            if(i)cout<<" + ";
            cout<<text[i];
        }
        cout<<endl;
    }

    void add_term(int term){ 
        terms[term].flip();
    }

    /*void add_terms(const vector<Term>&oth_terms){
        // if v contains pr , remove it
        // else add it
    }*/

    Expression operator+(const Expression& oth)const{
        Expression res=*this;
        //res.add_terms(oth.terms);
        res.terms^=oth.terms;
        return res;
    } 
    
    Expression operator*(const Expression& oth)const{
        Expression res;
        for(auto i : free_vars)if(terms[i]){
            for(auto j : free_vars)if(oth.terms[j]){
                int x=i,y=j;
                if(x>y)swap(x,y);
                res.terms[mapping[x][y]].flip();
            }
        }

        if(get_const()){
            res.terms^=oth.terms;
        }
        
        if(oth.get_const())
            res.terms^=terms;

        res.set_const(get_const()&oth.get_const());

        return res;
    }

    int eval(GF x){
        int ans=0;
        for(int i=0;i<128;i++)
            ans^=terms[i]&x[i];
        ans^=terms[QUAD_VARS];
        return ans;
    }
    /*Expression subs(vector<Expression> repr)const{
        Expression res;
        for(auto term : terms){
            if(term.var1==term.var2)
                res=res+repr[term.var1];
            else
                res=res+repr[term.var1]*repr[term.var2];
        }
        if(constant){
            res.constant^=1;
        }
        return res;
    }*/
};

typedef vector<Expression> Poly;



vector<bitset<n> >rand_mat(){// ensure that the matrix is invertible
    int r=0;
    vector<bitset<n> >vec;
    BitBasis<n> basis;
    std::mt19937 mt;
    while(r<n){
        bitset<n>bs;
        for(int i=0;i<n;i++)
            bs[i]=mt()%2;
        if(basis.insert(bs)){
            vec.push_back(bs);
            r++;
        }
    }
    return vec;
}



Poly multiply_matrix(vector<bitset<n> >mat,Poly a){
    Poly b;
    b.resize(mat.size());
    for(int i=0;i<mat.size();i++){
        for(int j=0;j<n;j++)if(mat[i][j]){
            b[i]=b[i]+a[j];
        }
    }
    return b;
}

/*
Poly multiply(Poly a,Poly b){
    Poly c;
    c.resize(a.size()+b.size()-1);

    for(int i=0;i<a.size();i++){
        for(int j=0;j<b.size();j++){
            Expression exp=a[i]*b[j];
            //add_pairs(c[i+j],terms);
            for(auto term : exp.terms)
                c[i+j].terms.push_back(term);
            c[i+j].constant^=exp.constant;
        }
    } 

    for(auto &cc:c)
        cc.simplify();

    // mod c by irr_poly
    for(int i=(int)c.size()-1;i>=n;i--){
        for(int j=0;j<n;j++)if(irr_poly[j]){
            //add_pairs(c[i-(n-j)],c[i]);
            for(auto term : c[i].terms)
                c[i-(n-j)].terms.push_back(term);
            c[i-(n-j)].constant^=c[i].constant;
        }
    }
    c.resize(n);
    for(auto &cc:c)
        cc.simplify();
    
    return c;
}*/

Poly multiply(Poly a,Poly b){
    Poly c;
    c.resize(a.size()+b.size()-1);

    for(int i=0;i<a.size();i++){
        for(int j=0;j<b.size();j++){
            Expression exp=a[i]*b[j];
            c[i+j].terms^=exp.terms;
        }
    } 
 

    // mod c by irr_poly
    for(int i=(int)c.size()-1;i>=n;i--){
        for(int j=0;j<n;j++)if(irr_poly[j]){
            c[i-(n-j)].terms^=c[i].terms;
        }
    }
    c.resize(n);
    
    return c;
}


Poly multiply_const(Poly a,Poly b){
    Poly c;
    c.resize(a.size()+b.size()-1);

    for(int i=0;i<a.size();i++)if(a[i].get_const()){
        for(int j=0;j<b.size();j++){
            Expression exp=b[j];
            c[i+j].terms^=b[j].terms;
        }
    } 


    // mod c by irr_poly
    for(int i=(int)c.size()-1;i>=n;i--){
        for(int j=0;j<n;j++)if(irr_poly[j]){
            c[i-(n-j)].terms^=c[i].terms;
        }
    }
    c.resize(n);
    
    return c;
}



Poly add(Poly a,Poly b){
    Poly c;
    c.resize(a.size());
    for(int i=0;i<a.size();i++){
        c[i]=a[i]+b[i];
    }
    
    return c;
}

Poly square(Poly a){
    Poly res;
    res.resize(n*2);
    for(int i=0;i<n;i++)
        res[i*2]=a[i];
    
    for(int i=(n-1)*2;i>=n;i--){
        for(int j=0;j<n;j++)if(irr_poly[j]){
            res[i-(n-j)].terms^=res[i].terms;
        }
    }
    res.resize(n);
    return res;
} 

Poly pow(Poly x,int t){
    //d=2^t
    //x^d
    Poly res;
    res.resize(n);
    res[0].set_const(1);

    auto tx=x;
    for(int i=0;i<t;i++)
        tx=square(tx);
    res=tx;

    return res;
}

int ans=0;

GF eval(Poly a,GF x){
    GF res;
    for(int i=0;i<128;i++){
        res[i]= a[i].eval(x);
    }
    return res;
}

void print(GF x){
    for(int i=0;i<128;i++)
        cout<<x[i];
    cout<<endl;
}

void generate(int _i){

    //assume x^d=alpha=a^d
    Poly x;
    Poly power_of_x[n]; 

    Poly a,alpha,one;
    one.resize(n);

    GF _a = _i;
    GF _alpha = GF_pow(_a,(1<<LOGD)+1);

    a.resize(n);
    alpha.resize(n);
    one[0].set_const(1);
    
    
    for(int i=0;i<n;i++){
        a[i].set_const(_a[i]);
        alpha[i].set_const(_alpha[i]);
    }

    x.resize(n);
    for(int i=0;i<n;i++){
        x[i]=Expression(i);
    }
    auto tx=x;
    for(int i=0;i<n;i++){ 
        power_of_x[i]=tx;  
        tx=square(tx); 
    }
    
    vector<Expression> equations;

    for(auto t : {T}){ 
        //auto equation=add(power_of_x[t],multiply(pow(alpha,2*LOGD,-1),x));  
        auto _const = GF_pow(_alpha,(1<<2*LOGD)-1);
        Poly c;
        c.resize(n);
        for(int i=0;i<n;i++)
            c[i].set_const(_const[i]);

        auto equation=add(power_of_x[t],multiply_const(c,x));  
        for(auto eq: equation){
            equations.push_back(eq);
        }
    }

    //cout<<"### Linear equations: ###"<<endl;
    //for(auto eq : equations)
    //    eq.print(); 

    BitBasisWithConstant<n> basis;
    //Finding free variables O(n^3)
    for(auto eq : equations){
        bitset<n+1>bs;
        for(int i=0;i<n;i++)if(eq.terms[i])
            bs[i]=1;
        bs[n]=eq.get_const();
        basis.insert(bs);
    }
    //print basis
    /*for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<basis.bs[i][j];
        }
        cout<<endl;
    }*/

    auto basic_vars=basis.basic_vars();
    free_vars=basis.free_vars();
    int rank=basis.rank();
/*
    cout<<"rank = "<<rank<<endl;
    for(auto var : basic_vars){
        cout<<"basic var: "<<name(var)<<endl;
    }
    for(auto var : free_vars){
        cout<<"free var: "<<name(var)<<endl;
    }*/
    vector<Expression> repr_x;// represent basic vars as linear combination of free vars O(n^2)
    repr_x.resize(n);
    for(int i=n-1;i>=0;i--){
        if(basis.bs[i][i]){
            for(int j=i+1;j<n;j++){
                if(basis.bs[i][j]){
                    repr_x[i]=repr_x[i]+repr_x[j];
                }
            }
            if(basis.bs[i][n])
                repr_x[i].add_one();
        }else{
            repr_x[i]=Expression(i);
        }
    }
    /*cout<<"### represent vars as free vars ###"<<endl;
    for(int i=0;i<n;i++){
        auto r=repr_x[i];
        cout<<name(i)<<" = ";
        r.print();
    }*/
 

    vector<Poly> power_of_repr_x;
    power_of_repr_x.resize(n); 

    tx=repr_x;
    for(int i=0;i<n;i++){ 
        //cerr<<"calc x^2^"<<i<<" which is represented by free vars"<<endl;
        power_of_repr_x[i]=tx; 
        tx=square(tx);
    }


    //cout<<"### finding quadratic equations ###"<<endl; 
    //auto mat = rand_mat(); 


    //output mat
    /*cout<<"M: "<<endl;
    for(int i=0;i<n;i++){
        for(int j=0;j<n;j++){
            cout<<mat[i][j];
        }
        cout<<endl;
    }*/

    
    Poly c1_P,c2,c;
    c1_P.resize(n);
    c2.resize(n);
    c.resize(n);
    for(int i=0;i<n;i++){
        c1_P[i].set_const(rain.c1[i]^P[i]);
        c2[i].set_const(rain.c2[i]);
        c[i].set_const(C[i]);
    }

    //b[0].set_const(1);
    //c[0].set_const(1);

    /*cout<<"b: "<<endl;
    for(int i=0;i<n;i++){
        cout<<b[i].constant;
    }
    cout<<endl;
    cout<<"c: "<<endl;
    for(int i=0;i<n;i++){
        cout<<c[i].constant;
    }
    cout<<endl;*/

    equations.clear();

    vector<Poly> equations_over_GF2n;
  
    
    memset(mapping,-1,sizeof(mapping));
    
    set<int>names;
    for(int i=0;i<QUAD_VARS;i++)    
        names.insert(i);

    for(auto v1 : free_vars){
        mapping[v1][v1]=v1;
        names.erase(v1);
        rmapping[v1]=make_pair(v1,v1);
    }

    for(auto v1 : free_vars)
    for(auto v2 : free_vars){
        if(v1<v2){
            mapping[v1][v2]=*names.begin();
            rmapping[*names.begin()]=make_pair(v1,v2);
            names.erase(names.begin());
        }
    }



    // alpha
    auto K=add(repr_x,c1_P);



    auto before = add( multiply_matrix(M1,pow(repr_x,LOGD)),add(K,c2));  
    auto after = add(K,c);

    
    /*puts("before at X");
    auto eb = eval(before,X);
    print(eb);
    puts("after at X");
    auto ea = eval(after,X);
    print(ea);

    puts("mul");
    print(GF_mul(eb,ea));

    puts("mul inv");
    print(GF_mul(eb,GF_inv(eb)));*/


    auto rx=add(multiply(before,after),one);   
    auto r2x_r=add(multiply(square(before),after),before); 
    auto rx2_x=add(multiply(before,square(after)),after);  

    equations_over_GF2n.push_back(rx); 
    equations_over_GF2n.push_back(r2x_r);
    equations_over_GF2n.push_back(rx2_x); 


    //rx[0].print(true);

    //Finding l
    //Solving final quadratic equations costs O(quad_vars^3)

    for(auto equation : equations_over_GF2n){
        for(auto eq: equation){ 
            equations.push_back(eq);
        }
    }

    for(auto eq : equations){
       // eq.print();
    }


    BitBasisWithConstant<QUAD_VARS> basis2;
    for(auto eq : equations){
        bitset<QUAD_VARS+1>bs;
        for(int i=0;i<QUAD_VARS;i++)if(eq.terms[i])
            bs[i]=1;
        bs[QUAD_VARS]=eq.get_const();
        basis2.insert(bs);
    }

    //ans^=basis2.rank();

    bitset<QUAD_VARS>val;
    
    for(int i=QUAD_VARS-1;i>=0;i--){
        if(basis2.bs[i][i]){
            for(int j=i+1;j<QUAD_VARS;j++){
                if(basis2.bs[i][j]){
                    val[i]=val[i]^val[j];
                }
            }
            if(basis2.bs[i][QUAD_VARS])
                val[i].flip();
        }else{
            throw;
        }
    }


    for(int i=n-1;i>=0;i--){
        if(basis.bs[i][i]){
            for(int j=i+1;j<n;j++)
                if(basis.bs[i][j])
                    val[i]=val[i]^val[j];
            if(basis.bs[i][n])
                val[i].flip();
        }else{
            //val
        }
    }




    GF candidate_key;

    for(int i=0;i<128;i++)
        candidate_key[i]=val[i]^rain.c1[i]^P[i];

    auto _c = rain.enc_2r(candidate_key,P);
    if (_c==C){
        cout<<"found key"<<endl;
        print(candidate_key);
    }

    /*cout<<"candidate key =";
    for(int i=0;i<128;i++)
        cout<<candidate_key[i];
    cout<<endl;*/

    /*cout<<"Field : GF(2^"<<n<<")"<<endl; 
    cout<<"Guess : x^(2^"<<LOGD<<"+1)"<<endl; 
    cout<<"number of free variable : "<<n-rank<<endl;
    cout<<"number of quadratic variables : "<<QUAD_VARS<<endl;
    cout<<"number of linear independent quadratic equations : "<<basis2.rank()<<endl;
    if(basis2.rank() >= QUAD_VARS) 
        cout<<"Complexity : 2^"<<(n-LOGD+3*(log2(max(n,QUAD_VARS))-log2(n)))<<endl; 
    else
        cout<<"Unsolvable"<<endl;*/
}


int main(){
    srand(123);
    set_parameters();
    //for(int i=0;i<100;i++){
    //    cout<<"checking "<<i<<endl;
    //}
    //count();

    
    rain.simple();

    M1=rain.M1;
    M1_inv = mat_inv(M1);
    P=GF_rand();
    K=GF_rand();
    
    //X=K^P^rain.c1;

    K = X^rain.c1^P;
    
    cout<<"K =";
    for(int i=0;i<128;i++)
        cout<<K[i];
    cout<<endl;
      
    C = rain.enc_2r(K,P);
    
    cout<<"C = ";
    for(int i=0;i<128;i++)
        cout<<C[i];
    cout<<endl;
   

    double st=clock();
    for(int i=1;i<=1000;i++)
        generate(i);
    cout<<"time = "<<(clock()-st)/CLOCKS_PER_SEC<<endl;

    return 0;
}