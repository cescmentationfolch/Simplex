#include "Matrix.hh"

Matrix::Matrix() {
    n = m = 0;
}

Matrix::Matrix(int _n, int _m, vector<vector<long double> > _M) {
    n = _n;
    m = _m;
    M = _M;
    N = vector<int>(m, 1);
}

Matrix::Matrix(int _n, bool _h, vector<long double> _M) {
  if (_h) {
    n = 1;
    m = _n;
    M = vector<vector<long double> >(1, _M);
  }
  else {
    n = _n;
    m = 1;
    M = vector<vector<long double> >(n, vector<long double>(1));
    for (int i = 0; i < n; ++i) M[i][0] = _M[i];
  }
  N = vector<int>(m, 1);
}

Matrix::~Matrix(){}

void Matrix::modifica_estat_columna(int f, int val) {
    N[f] = val;
}

Matrix Matrix::operator!() {
    vector<vector<long double> > Inv(n, vector<long double>(n, 0));
    for (int i = 0; i < n; ++i) Inv[i][i] = 1;
    vector<vector<long double> > B(n, vector<long double>(n));
    for (int i = 0; i < n; ++i) {
        int k = 0;
        for (int j = 0; j < m; ++j) {
            if (N[j]) B[i][k++] = M[i][j];
        }
    }
    for (int i = 0; i < n; ++i) {
        int k = i;
        long double piv = B[i][i];
        for (int j = i; j < n; ++j) {
            if (abs(B[j][i]) > abs(piv)) {
                piv = B[j][i];
                k = j;
            }
        } 
        for (int j = 0; j < n; ++j) {
	    swap(B[i][j], B[k][j]);
	    swap(Inv[i][j], Inv[k][j]);
	    B[i][j] /= piv;
	    Inv[i][j] /= piv;
	}
        for (int j = 0; j < n; ++j) {
	    if (i == j) continue;
            long double mult = B[j][i];
            for (int w = i; w < n; ++w) B[j][w] -= mult*B[i][w];
            for (int w = 0; w < n; ++w) Inv[j][w] -= mult*Inv[i][w];
        }
    }
    return Matrix(n,n,Inv);
}

Matrix Matrix::operator*(const Matrix& mat) {
    vector<int> VM;
    for (int i = 0; i < m; ++i)
        if (N[i]) VM.push_back(i);
    vector<int> Vmat;
    for (int i = 0; i < mat.m; ++i) 
        if (mat.N[i]) Vmat.push_back(i);
    vector<vector<long double> > Res(n, vector<long double>((int)Vmat.size(), 0));
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < (int)Vmat.size(); ++j) 
            for (int k = 0; k < mat.n; ++k) 
                Res[i][j] += M[i][VM[k]]*mat.M[k][Vmat[j]];
    return Matrix(n, (int)Vmat.size(), Res);
}

Matrix Matrix::operator-(const Matrix& mat) {
     vector<int> VM;
    for (int i = 0; i < m; ++i)
        if (N[i]) VM.push_back(i);
    vector<int> Vmat;
    for (int i = 0; i < mat.m; ++i) 
        if (mat.N[i]) Vmat.push_back(i);
    vector<vector<long double> > Res(n, vector<long double>((int)VM.size()));
    for (int i = 0; i < n; ++i)
	for (int j = 0; j < (int)VM.size(); ++j)
	    Res[i][j] = M[i][VM[j]] - mat.M[i][Vmat[j]];
    return Matrix(n, (int)Vmat.size(), Res);
}

Matrix Matrix::operator+(const Matrix& mat) {
     vector<int> VM;
    for (int i = 0; i < m; ++i)
        if (N[i]) VM.push_back(i);
    vector<int> Vmat;
    for (int i = 0; i < mat.m; ++i) 
        if (mat.N[i]) Vmat.push_back(i);
    vector<vector<long double> > Res(n, vector<long double>((int)VM.size()));
    for (int i = 0; i < n; ++i)
	for (int j = 0; j < (int)VM.size(); ++j)
	    Res[i][j] = M[i][VM[j]] + mat.M[i][Vmat[j]];
    return Matrix(n, (int)Vmat.size(), Res);
}

void Matrix::print() {
    cout << "Matrix:\n";
    for (int i = 0; i < m; ++i) {
	if (i) cout << ' ';
	cout << N[i];
    }
    cout << '\n';
    for (int i = 0; i < n; ++i) {
	for (int j = 0; j < m; ++j) {
	    if (j) cout << ' ';
	    if (abs(M[i][j]) > 1e-15) cout << M[i][j];
	    else cout << 0;
	}
	cout << '\n';
    }
}