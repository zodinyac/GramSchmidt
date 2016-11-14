#include <iostream>
#include <math.h>

using namespace std;

const double eps = 1e-7;

void assignment(double *destination, double *source, int length) {
	while (length > 0) {
		*destination = *source;
		++source;
		++destination;
		--length;
	}
}

void proj_and_minus(double *destination, double *a, double *b, int length) {
	double top(0), bottom(0);
	double *tmp_b = b;
	for (int i = 0;  i < length; ++i, ++a, ++tmp_b) {
		top += *a * *tmp_b;
		bottom += *tmp_b * *tmp_b;
	}
	double s = top / bottom;
	for (int i = 0;  i < length; ++i, ++destination, ++b)
		*destination -= s * *b;
}

bool is_non_zero(double *a, int length) {
	for (int i = 0; i < length; ++i)
		if (fabs(a[i]) > eps)
			return true;
	return false;
}

int main(int argc, char* argv[])
{
	int n, m;
	cin >> n >> m;
	double **a, **b;
	bool *non_zero_b;
	a = new double*[m];
	b = new double*[m];
	non_zero_b = new bool[m];
	for (int i = 0; i < m; ++i) {
		a[i] = new double[n];
		b[i] = new double[n];
	}

	for (int i = 0; i < m; ++i)
		for (int j = 0; j < n; ++j)
			cin >> a[i][j];

	assignment(b[0], a[0], n);
	non_zero_b[0] = is_non_zero(b[0], n);
	for (int i = 1; i < m; ++i) {
		assignment(b[i], a[i], n);
		for (int j = 0; j < i; ++j)
			if (non_zero_b[j])
				//proj_and_minus(b[i], a[i], b[j], n); // default scheme
				proj_and_minus(b[i], b[i], b[j], n); // pro scheme
		non_zero_b[i] = is_non_zero(b[i], n);
	}

	for (int i = 0; i < m; ++i)
		if (non_zero_b[i]) {
			for (int j = 0; j < n; ++j)
				cout << b[i][j] << " ";
			cout << endl;
		}

	for (int i = 0; i < m; ++i) {
		delete[] a[i];
		delete[] b[i];
	}
	delete[] a;
	delete[] b;
	delete[] non_zero_b;

	return 0;
}
