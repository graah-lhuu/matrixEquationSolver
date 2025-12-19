#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <string>

using namespace std;

const double EPS = 1e-10; // 判断是否为0的容差

// 判断一个数是否接近0 (考虑浮点数精度)
bool isZero(double val) {
	return fabs(val) < EPS;
}

// 打印矩阵
void printMatrix(const vector<vector<double>>& matrix, const string& name = "") {
	if (!name.empty()) {
		cout << name << ":" << endl;
	}
	for (const auto& row : matrix) {
		for (double val : row) {
			if (isZero(val)) {
				cout << setw(10) << "0.00";
			} else {
				cout << setw(10) << fixed << setprecision(4) << val;
			}
		}
		cout << endl;
	}
	cout << endl;
}

// 高斯-约当消元法求简化行阶梯形，并记录行变换矩阵E
void gaussJordanElimination(vector<vector<double>>& augmented, vector<vector<double>>& E) {
	int rows = augmented.size();
	if (rows == 0) return;
	int cols = augmented[0].size();
	int rank = 0; // 矩阵的秩

	// 初始化行变换记录矩阵E为单位矩阵
	E.assign(rows, vector<double>(rows, 0.0));
	for (int i = 0; i < rows; i++) {
		E[i][i] = 1.0;
	}

	for (int col = 0, row = 0; col < cols - 1 && row < rows; col++) { // 最后一列是b，不作为主元列
		// 1. 选主元 (当前列中绝对值最大的行)
		int pivotRow = row;
		for (int i = row + 1; i < rows; i++) {
			if (fabs(augmented[i][col]) > fabs(augmented[pivotRow][col])) {
				pivotRow = i;
			}
		}

		if (isZero(augmented[pivotRow][col])) {
			continue; // 当前列全为0，检查下一列
		}

		// 2. 交换当前行和主元行 (在增广矩阵和行变换矩阵E中同时进行)
		if (pivotRow != row) {
			swap(augmented[row], augmented[pivotRow]);
			swap(E[row], E[pivotRow]);
		}

		// 3. 主元归一化
		double pivotVal = augmented[row][col];
		for (int j = col; j < cols; j++) {
			augmented[row][j] /= pivotVal;
		}
		for (int j = 0; j < rows; j++) { // E矩阵的列数与行数相同
			E[row][j] /= pivotVal;
		}

		// 4. 将主元列的其他元素消为0
		for (int i = 0; i < rows; i++) {
			if (i != row && !isZero(augmented[i][col])) {
				double factor = augmented[i][col];
				for (int j = col; j < cols; j++) {
					augmented[i][j] -= factor * augmented[row][j];
				}
				for (int j = 0; j < rows; j++) { // 对E矩阵进行相同的行操作
					E[i][j] -= factor * E[row][j];
				}
			}
		}
		row++;
		rank++;
	}
}

// 分析解的情况并求解
void analyzeAndSolve(const vector<vector<double>>& A, const vector<double>& b) {
	int m = A.size();    // 方程个数
	int n = (m > 0) ? A[0].size() : 0; // 未知数个数

	// 构建增广矩阵 [A | b]
	vector<vector<double>> augmented = A;
	for (int i = 0; i < m; i++) {
		augmented[i].push_back(b[i]);
	}

	cout << "\n原始增广矩阵 [A | b]:" << endl;
	printMatrix(augmented);

	// 记录行变换的矩阵E
	vector<vector<double>> E;
	// 对增广矩阵进行高斯-约当消元，化为简化行阶梯形，并记录行变换E
	gaussJordanElimination(augmented, E);

	cout << "行变换矩阵 E (使得 E * [A|b] = RREF):" << endl;
	printMatrix(E, "E");

	cout << "简化行阶梯形矩阵 RREF([A|b]):" << endl;
	printMatrix(augmented, "RREF");

	// 判断解的情况
	int rankA = 0, rankAB = 0;

	// 系数矩阵A的秩：RREF中前n列非零行的行数
	for (int i = 0; i < m; i++) {
		bool allZeroA = true;
		for (int j = 0; j < n; j++) {
			if (!isZero(augmented[i][j])) {
				allZeroA = false;
				break;
			}
		}
		if (!allZeroA) rankA++;
	}

	// 增广矩阵的秩：RREF中全部n+1列非零行的行数
	for (int i = 0; i < m; i++) {
		bool allZeroAB = true;
		for (int j = 0; j < n + 1; j++) { // 包括常数项列
			if (!isZero(augmented[i][j])) {
				allZeroAB = false;
				break;
			}
		}
		if (!allZeroAB) rankAB++;
	}

	cout << "系数矩阵 A 的秩: " << rankA << endl;
	cout << "增广矩阵 [A|b] 的秩: " << rankAB << endl << endl;

	// 检查无解情况：系数矩阵的秩 < 增广矩阵的秩
	if (rankA < rankAB) {
		cout << "该线性方程组无解。" << endl;
		for (int i = rankA; i < m; i++) {
			if (!isZero(augmented[i][n])) { // 找到矛盾方程，例如 0 = 1
				cout << "矛盾方程出现在第 " << i + 1 << " 行: 0 = " << augmented[i][n] << endl;
				break;
			}
		}
		return;
	}

	// 检查唯一解/无穷多解
	if (rankA == n) {
		cout << "该线性方程组有唯一解。" << endl;
		// 唯一解可以直接从RREF中读出
		vector<double> x_unique(n, 0.0);
		vector<int> pivotCols; // 记录主元列
		int r = 0;
		for (int j = 0; j < n && r < m; j++) {
			if (!isZero(augmented[r][j])) {
				pivotCols.push_back(j);
				x_unique[j] = augmented[r][n]; // 常数项列的值
				r++;
			}
		}
		cout << "解为: ";
		for (int i = 0; i < n; i++) {
			cout << "x" << i + 1 << " = " << x_unique[i];
			if (i < n - 1) cout << ", ";
		}
		cout << endl;

	} else { // rankA < n, 无穷多解
		cout << "该线性方程组有无穷多解。" << endl;
		cout << "自由变量（自由未知量）个数为: " << n - rankA << endl << endl;

		// 寻找主元列和自由列
		vector<bool> isPivotCol(n, false);
		vector<int> pivotRows(n, -1); // 记录每个主元所在的行，-1表示该列不是主元列
		int currentRow = 0;
		for (int j = 0; j < n; j++) {
			if (currentRow < m && !isZero(augmented[currentRow][j])) {
				isPivotCol[j] = true;
				pivotRows[j] = currentRow; // 第j列的主元在第currentRow行
				currentRow++;
			}
		}

		// 1. 求特解 Xp：令所有自由变量为0
		vector<double> x_particular(n, 0.0);
		for (int j = 0; j < n; j++) {
			if (isPivotCol[j]) {
				int row = pivotRows[j];
				x_particular[j] = augmented[row][n]; // 主元变量的值就是对应行的常数项
			}
			// 自由变量默认为0
		}
		cout << "一个特解 Xp (令自由变量为0): " << endl;
		for (int i = 0; i < n; i++) {
			cout << "  x" << i + 1 << " = " << x_particular[i] << endl;
		}

		// 2. 求零空间基向量 Xn：分别令每个自由变量为1，其余自由变量为0，回代求解
		cout << endl << "零空间的基础解系（自由变量依次取1，其余自由变量取0）:" << endl;
		vector<vector<double>> nullspaceBasis;

		for (int freeIndex = 0; freeIndex < n; freeIndex++) {
			if (!isPivotCol[freeIndex]) { // 对于每一个自由变量列
				vector<double> x_null(n, 0.0);
				x_null[freeIndex] = 1.0; // 设置该自由变量为1

				// 回代求解主元变量
				for (int j = 0; j < n; j++) {
					if (isPivotCol[j]) {
						int row = pivotRows[j];
						x_null[j] = -augmented[row][freeIndex]; // 主元变量 = - (自由变量系数)
					}
				}
				nullspaceBasis.push_back(x_null);

				cout << "基础解系向量（对应自由变量 x" << freeIndex + 1 << "=1）: ";
				for (double val : x_null) {
					cout << setw(8) << val << " ";
				}
				cout << endl;
			}
		}

		// 3. 通解形式: X = Xp + c1*Xn1 + c2*Xn2 + ...
		cout << endl << "方程组的通解为: X = Xp + ";
		int freeVarCount = 0;
		for (int i = 0; i < n; i++) {
			if (!isPivotCol[i]) {
				if (freeVarCount > 0) cout << " + ";
				cout << "c" << freeVarCount + 1 << " * (Xn" << freeVarCount + 1 << ")";
				freeVarCount++;
			}
		}
		if (freeVarCount == 0) cout << "0"; // 实际上不会发生，因为rankA<n意味着有自由变量
		cout << endl;
		cout << "其中 c1, c2, ... 为任意常数。" << endl;
	}
}

// 交互式输入函数
void interactiveInput() {
	int m, n;

	cout << "================================" << endl;
	cout << "线性方程组求解器" << endl;
	cout << "求解 A * X = b" << endl;
	cout << "================================" << endl;

	// 输入矩阵大小
	cout << "请输入方程个数 m (矩阵A的行数): ";
	cin >> m;
	cout << "请输入未知数个数 n (矩阵A的列数): ";
	cin >> n;

	if (m <= 0 || n <= 0) {
		cout << "错误：m和n必须是正整数！" << endl;
		return;
	}

	// 输入矩阵A
	vector<vector<double>> A(m, vector<double>(n));
	cout << "\n请输入系数矩阵 A (" << m << "行 × " << n << "列):" << endl;
	for (int i = 0; i < m; i++) {
		cout << "输入第 " << i + 1 << " 行的 " << n << " 个元素: ";
		for (int j = 0; j < n; j++) {
			cin >> A[i][j];
		}
	}

	// 输入向量b
	vector<double> b(m);
	cout << "\n请输入常数向量 b (" << m << "个元素):" << endl;
	for (int i = 0; i < m; i++) {
		cout << "输入 b[" << i + 1 << "]: ";
		cin >> b[i];
	}

	// 显示输入的矩阵和向量
	cout << "\n您输入的系数矩阵 A:" << endl;
	printMatrix(A, "A");

	cout << "您输入的常数向量 b:" << endl;
	cout << "b = [ ";
	for (int i = 0; i < m; i++) {
		cout << b[i];
		if (i < m - 1) cout << ", ";
	}
	cout << " ]^T" << endl;

	// 分析并求解
	cout << "\n" << string(50, '=') << endl;
	cout << "求解结果:" << endl;
	cout << string(50, '=') << endl;

	analyzeAndSolve(A, b);
}

int main() {
	int choice;

	cout << "线性方程组求解器" << endl;
	cout << "==================" << endl;
	cout << "请选择操作模式:" << endl;
	cout << "1. 交互式输入（手动输入矩阵A和向量b）" << endl;
	cout << "2. 运行示例" << endl;
	cout << "请输入选择 (1 或 2): ";
	cin >> choice;

	if (choice == 1) {
		interactiveInput();
	} else if (choice == 2) {
		// 运行预设示例
		cout << "\n运行预设示例..." << endl;

		cout << "\n=== 示例1: 唯一解 ===" << endl;
		vector<vector<double>> A1 = {
			{2, -1, 3, 2},
			{3, -3, 3, 2},
			{3, -1, -1, 2},
			{3, -1, 1, -1}
		};
		vector<double> b1 = { 6, 5, 3, 2 };
		analyzeAndSolve(A1, b1);

		cout << "\n" << string(50, '=') << endl;

		cout << "\n=== 示例2: 无穷多解 ===" << endl;
		vector<vector<double>> A2 = {
			{1, 0, 0, 1, 0, 0},
			{0, 1, 0, 0, 1, 0},
			{0, 0, 1, 0, 0, 1},
			{1, 1, 1, 0, 0, 0},
			{0, 0, 0, 1, 1, 1}
		};
		vector<double> b2 = { 40, 20, 10, 45, 25 };
		analyzeAndSolve(A2, b2);

		cout << "\n" << string(50, '=') << endl;

		cout << "\n=== 示例3: 无解 ===" << endl;
		vector<vector<double>> A3 = {
			{1, 2},
			{1, 1},
			{2, 1}
		};
		vector<double> b3 = { 5, 4, 3 };
		analyzeAndSolve(A3, b3);
	} else {
		cout << "无效选择，程序结束。" << endl;
	}

	// 等待用户查看结果
	cout << "\n" << string(50, '=') << endl;
	cout << "程序运行结束，按Enter键退出...";
	cin.ignore(); // 清除输入缓冲区
	cin.get();    // 等待用户按Enter键

	return 0;
}