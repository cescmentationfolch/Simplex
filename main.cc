#include "Matrix.hh"

typedef long double ld;
typedef vector<ld> vld;
typedef vector<vld> vvld;

int main() {
    vector<Matrix> VM;
    char a;
    while (cin >> a) {
	if (a == 'N') {
	    int n, m;
	    cin >> n >> m;
	    vvld M(n, vld(m));
	    for (int i = 0; i < n; ++i) 
		for (int j = 0; j < m; ++j)
		    cin >> M[i][j];
	    VM.push_back(Matrix(n, m, M));
	}
	else if (a == 'P') {
	    int x;
	    cin >> x;
	    VM[x].print();
	}
	else if (a == 'I') {
	    int x;
	    cin >> x;
	    VM[x].print();
	    (!VM[x]).print();
	    (VM[x]*(!VM[x])).print();
	}
    }
}