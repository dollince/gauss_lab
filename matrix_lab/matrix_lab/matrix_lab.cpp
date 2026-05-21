#include<iostream>
#include<cmath>
#include<stdexcept>
#include<clocale>

template<typename T>
class Vector {
protected: 
	T* data;
	int size;
public:
	Vector (int n = 0) : size(n) {
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = T(0);
	}
	Vector (const Vector& v) : size(v.size) {
		data = new T[size];
		for (int i = 0; i < size; i++)
			data[i] = v.data[i];
	}
	Vector& operator=(const Vector& v) {
		if (this == &v)
			return *this;
		delete[] data;
		size = v.size;
		data = new T[size];
		for (int i = 0; i < size; ++i)
			data[i] = v.data[i];
		return *this;
	}
	~Vector() {
		delete[] data;
	}
	int getSize() const {
		return size;
	}
	T& operator[](int i) {
		return data[i];
	}
	void print() {
		std::cout << "[";
		for (int i = 0; i < size; i++)
			std::cout << data[i]<< " ";
		std::cout << "]\n";
	}
};

template<typename T>
class Matrix : public Vector<Vector<T>> {
private:
	int lines; //строки
	int cols; //столбцы
public:
	Matrix(int lines, int cols) : Vector<Vector<T>>(lines), lines(lines), cols(cols) {
		for (int i = 0; i < lines; i++)
			(*this)[i] = Vector<T>(cols);
	}
	int getLines() const {
		return lines;
	}
	int getCols() const {
		return cols;
	}
	void print() {
		for (int i = 0; i < lines; i++) {
			std::cout << " ";
			(*this)[i].print();
		}
	}
};


template<typename T>
Vector<T> gauss(Matrix<T> A, Vector<T> b) {
	int n = A.getLines(); //уравнения
	int m = A.getCols(); //неизвестные

	T EPS = 1e-9;

	Vector<int> pivot_col(n);
	for (int i = 0; i < n; i++)
		pivot_col[i] = -1;

	int nowline = 0;

	for (int nowcol = 0; nowcol < m && nowline < n; nowcol++) {
		int pivotline = nowline;
		for (int i = nowline + 1; i < n; i++) {
			if (std::abs(A[i][nowcol]) > std::abs(A[pivotline][nowcol])) {
				pivotline = i;
			}
		}
		if (std::abs(A[pivotline][nowcol]) < EPS)
			continue;

		if (pivotline != nowline) {
			Vector<T> temp = A[nowline];
			A[nowline] = A[pivotline];
			A[pivotline] = temp;

			T tempb = b[nowline];
			b[nowline] = b[pivotline];
			b[pivotline] = tempb;
		}

		T pivot = A[nowline][nowcol];

		for (int j = nowcol; j < m; j++)
			A[nowline][j] /= pivot;
		b[nowline] /= pivot;

		//обнуление
		for (int i = 0; i < n; i++) {
			if (i == nowline)
				continue;
			T factor = A[i][nowcol];

			for (int j = nowcol; j < m; j++)
				A[i][j] -= factor * A[nowline][j];
			b[i] -= factor * b[nowline];
		}
		pivot_col[nowline] = nowcol;
		nowline++;

	}
	int rank = nowline;
	// проверка на несовместную систему
	for (int i = rank; i < n; i++) {
		if (std::abs(b[i]) > EPS)
			throw std::runtime_error("Ошибка: система несовместна.");
	}
	Vector<int> free_index(m);
	for (int j = 0; j < m; j++)
		free_index[j] = -1;
	for (int i = 0; i < rank; i++)
		if (pivot_col[i] != -1)
			free_index[pivot_col[i]] = -2;

	int num_free = 0;
	for (int j = 0; j < m; j++) {
		if (free_index[j] == -1) {  
			free_index[j] = num_free;
			num_free++;
		}
	}

	if (num_free == 0) {
		Vector<T> x(m);
		for (int i = 0; i < rank; i++)
			x[pivot_col[i]] = b[i];
		return x;
	}

	Vector<T> particular(m);       
	Matrix<T> coeffs(m, num_free); 

	
	for (int j = 0; j < m; j++) {
		if (free_index[j] >= 0)  
			coeffs[j][free_index[j]] = T(1);
	}

	
	for (int i = 0; i < rank; i++) {
		int pc = pivot_col[i];
		if (pc == -1) continue;
		particular[pc] = b[i];
		for (int j = 0; j < m; j++) {
			if (free_index[j] >= 0)  
				coeffs[pc][free_index[j]] = -A[i][j];
		}
	}

	std::cout << "\nБесконечно много решений. Параметрическое решение:\n";
	for (int j = 0; j < m; j++) {
		std::cout << "x[" << j + 1 << "] = " << particular[j];
		for (int k = 0; k < num_free; k++) {
			T c = coeffs[j][k];
			if (std::abs(c) < EPS) continue;
			if (c > 0)
				std::cout << " + " << c << "*t" << k + 1;
			else
				std::cout << " - " << -c << "*t" << k + 1;
		}
		std::cout << "\n";
	}
	

	throw std::runtime_error("Бесконечно много решений.");
}


template<typename T>
void inputoutput() {
	int n, m;
	std::cout << "Введите количество строк:";
	std::cin >> n;
	std::cout << "Введите количество столбцов:";
	std::cin >> m;

	if (n <= 0 || m <= 0) {
		std::cout << "Ошибка: размеры должны быть положительнымию\n";
		return;
	}
	Matrix<T> A(n, m);
	Vector<T> b(n);

	std::cout << "\nВведите матрицу А(" << n << "х" << m << ") построчно:\n";
	for (int i = 0; i < n; i++) {
		std::cout << "строка" << i << ":";
		for (int j = 0; j < m; j++)
			std::cin >> A[i][j];

	}
	std::cout << "\nВведите вектор b:\n";
	for (int i = 0; i < n; i++)
		std::cin >> b[i];

	std::cout << "\nМатрица A:\n";
	A.print();
	std::cout << "\nВектор b:\n";
	b.print();

	try {
		Vector<T> x = gauss(A, b);
		std::cout << "\nЕдинственное решение системы:\n";
		for (int i = 0; i < m; i++) {
			std::cout << "x[" << i+1 << "]=" << x[i] << "\n";

		}
	}
	catch (const std::exception& e) {
		std::cout << "\n" << e.what() << "\n";
	}
}

int main() {
	setlocale(LC_ALL, "rus");
	inputoutput<double>();
	return 0;
}