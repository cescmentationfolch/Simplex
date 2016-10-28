#ifndef MATRIX_HH
#define MATRIX_HH

#include<bits/stdc++.h>

using namespace std;

class Matrix {
  
    //classe per operar amb matrius
    
private:
    //vector que marca quines columnes es fan servir 
    vector<int> N;
    
    //elements de la matriu
    vector<vector<long double> > M;
    
    //dimensions de la matriu 
    int n, m;
    
public:
    //Constructores
    
    Matrix();
    //pre: cert
    //post: es crea una matriu amb dimensions (0,0)
    
    Matrix(int _n, int _m, vector<vector<long double> > _M);
    //pre: cert
    //post: es crea una matriu (_n,_m) amb els elements de _M
    
    //Destructores
    
    ~Matrix();
    
    //Modificadores
    
    void modifica_estat_columna(int f, int val);
    //pre: 0 <= f < m, val = 0 o 1
    //post: l'estat de la columna f passa ser val i s'actualitza mr
    
    Matrix operator!();
    //pre: el nombre de columnes actives es n i te rang maxim
    //post: retorna la matriu inversa (n,n)
    
    Matrix operator*(const Matrix& mat);
    //pre: la matriu mat es (m,k)
    //post: retorna la matriu resultant de multiplicar la matriu per mat
};

#endif