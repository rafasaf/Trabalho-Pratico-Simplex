#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#include <fstream>
using namespace std;

class SimplexSolver {
private:
    int numConstraints, numVariables;
    vector<double> costCoefficients;
    vector<vector<double>> constraintsMatrix;

public:
    SimplexSolver(int n, int m, const vector<double>& costs, const vector<vector<double>>& constraints)
        : numConstraints(n), numVariables(m), costCoefficients(costs), constraintsMatrix(constraints) {}

    bool isUnbounded(const vector<vector<double>>& tableau, int col) {
        for (int i = 0; i < tableau.size() - 1; i++) {
            if (tableau[i][col] >= 0)
                return false;
        }
        return true;
    }

    void resetNearZero(double& value) {
        if (value >= -1e-9 && value <= 1e-9) {
            value = 0.0;
        }
    }

    void printMatrix(const vector<vector<double>>& matrix) {
        int numRows = matrix.size();
        int numCols = matrix[0].size();

        cout << fixed << setprecision(3);

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                cout << matrix[i][j] << "\t";
            }
            cout << endl;
        }
        cout << endl;
    }

    vector<vector<double>> secondPhaseTableau(const vector<vector<double>>& firstTableau) {
        vector<vector<double>> secondTableau(numConstraints + 1, vector<double>(numConstraints + numVariables + 1, 0));
        int numRows = secondTableau.size();
        int numCols = secondTableau[0].size();

        for (int i = 0; i < numRows; i++) {
            for (int j = 0; j < numCols; j++) {
                if (i == 0)
                    secondTableau[i][j] = firstTableau[firstTableau.size() - 1][j];
                else
                    secondTableau[i][j] = firstTableau[i][j];
            }
        }

        for (int i = 0; i < numRows; i++)
            secondTableau[i][numCols - 1] = firstTableau[i][firstTableau[0].size() - 1];

        secondTableau[0][numCols - 1] = firstTableau[firstTableau.size() - 1][firstTableau[0].size() - 1];

        return secondTableau;
    }

    void printSolutions(const vector<vector<double>>& tableau) {
        vector<double> solutions(numVariables, 0);
        for (int i = numConstraints; i < numConstraints + numVariables; i++) {
            for (int j = 0; j < numConstraints + 1; j++) {
                if (tableau[j][i] == 1)
                    solutions[j] = tableau[j][tableau[0].size() - 1];
            }
        }

        for (auto sol : solutions)
            cout << sol << " ";
        cout << endl;
    }

    void printCertificate(const vector<vector<double>>& tableau, int col) {
        vector<double> certificate(numVariables, 0);
        int row;
        certificate[col - numConstraints] = 1;

        for (int i = numConstraints; i < numConstraints + numVariables; i++) {
            for (int j = 0; j < numConstraints + 1; j++) {
                if (tableau[j][i] == 1) {
                    row = j;
                    certificate[j] = -tableau[row][col];
                }
            }
        }

        for (auto cert : certificate)
            cout << cert << " ";
        cout << endl;
    }

    void printInfeasible(const vector<vector<double>>& tableau) {
        cout << "inviavel" << endl;
        for (int i = 0; i < tableau.size() - 2; i++)
            cout << tableau[0][i] << " ";
        cout << endl;
    }

    void makePositive() {
        for (int i = 0; i < numConstraints; i++) {
            if (constraintsMatrix[i][numVariables] < 0) {
                for (int j = 0; j < numVariables + 1; j++) {
                    constraintsMatrix[i][j] *= -1;
                }
            }
        }
    }

    vector<vector<double>> auxiliaryProblem() {
        vector<vector<double>> auxTableau(numConstraints + 1, vector<double>(numConstraints + numVariables + 1, 0));

        for (int i = 0; i < numConstraints + numVariables; i++) {
            if (i >= numVariables)
                auxTableau[0][i] = 1;
        }

        for (int i = 1; i < numConstraints + 1; i++) {
            for (int j = 0; j < constraintsMatrix[0].size() - 1; j++) {
                auxTableau[i][j] = constraintsMatrix[i - 1][j];
            }
        }

        for (int i = 1; i < numConstraints + 1; i++) {
            auxTableau[i][numVariables + i - 1] = 1;
        }

        for (int i = 1; i < numConstraints + 1; i++) {
            auxTableau[i][numConstraints + numVariables] = constraintsMatrix[i - 1][constraintsMatrix[0].size() - 1];
        }

        return auxTableau;
    }

    vector<vector<double>> buildTableau(const vector<vector<double>>& auxTableau) {
        vector<vector<double>> tableau(numConstraints + 2, vector<double>(2 * numConstraints + numVariables + 1, 0));
        int numRows = tableau.size();
        int numCols = tableau[0].size();

        for (int i = 0; i < auxTableau.size(); i++) {
            for (int j = numConstraints; j < numCols; j++) {
                tableau[i][j] = auxTableau[i][j - numConstraints];
            }
        }

        for (int i = 1; i < numRows - 1; i++) {
            tableau[i][i - 1] = 1;
        }

        for (int j = 0; j < costCoefficients.size(); j++) {
            tableau[numConstraints + 1][j + numConstraints] = -costCoefficients[j];
        }

        return tableau;
    }

    pair<int, int> findPivot(const vector<vector<double>>& tableau, bool firstPhase) {
        int numRows = firstPhase ? tableau.size() - 1 : tableau.size();
        int numCols = tableau[0].size();

        int pivotCol = -1;
        double minValue = 0;
        for (int j = numConstraints; j < numConstraints + numVariables; j++) {
            if (tableau[0][j] < minValue) {
                minValue = tableau[0][j];
                pivotCol = j;
            }
        }

        int pivotRow = -1;
        double minRatio = 1e+10;

        for (int i = 0; i < numRows; i++) {
            double valueInPivotCol = tableau[i][pivotCol];
            if (valueInPivotCol > 0) {
                double ratio = tableau[i][numCols - 1] / valueInPivotCol;
                if (ratio < minRatio) {
                    minRatio = ratio;
                    pivotRow = i;
                }
            }
        }

        return {pivotRow, pivotCol};
    }

    void completePivoting(vector<vector<double>>& tableau, int pivotRow, int pivotCol) {
        int numRows = tableau.size();
        int numCols = tableau[0].size();

        double pivot = tableau[pivotRow][pivotCol];

        for (int j = 0; j < numCols; ++j) {
            tableau[pivotRow][j] /= pivot;
            resetNearZero(tableau[pivotRow][j]);
        }

        for (int i = 0; i < numRows; ++i) {
            if (i != pivotRow) {
                double factor = tableau[i][pivotCol];
                for (int j = 0; j < numCols; ++j) {
                    tableau[i][j] -= factor * tableau[pivotRow][j];
                    resetNearZero(tableau[i][j]);
                }
            }
        }
    }

    void initialPivoting(vector<vector<double>>& tableau) {
        int col = numVariables + numConstraints;
        for (int i = 1; i < numConstraints + 1; i++) {
            completePivoting(tableau, i, col);
            col++;
        }
    }

    void printUnbounded(vector<vector<double>>& tableau, int col) {
        cout << "ilimitada" << endl;
        printSolutions(tableau);
        printCertificate(tableau, col);
    }

    int simplex(vector<vector<double>>& tableau, bool firstPhase) {
        if (firstPhase)
            initialPivoting(tableau);

        while (true) {
            pair<int, int> pivotCoords = findPivot(tableau, firstPhase);

            if (pivotCoords.first == -1 || pivotCoords.second == -1) {
                break;
            }

            completePivoting(tableau, pivotCoords.first, pivotCoords.second);
        }

        double optimalValue = tableau[0][tableau[0].size() - 1];
        resetNearZero(optimalValue);

        if (optimalValue == 0)
            return 1;
        else
            return 0;
    }

    void printOptimal(vector<vector<double>>& tableau) {
        cout << "otima" << endl << tableau[0][tableau[0].size() - 1] << endl;

        vector<double> solutions;

        for (int i = numConstraints; i < numConstraints + numVariables; i++) {
            bool found = false;
            for (int j = 0; j < numConstraints + 1; j++) {
                if (tableau[j][i] == 1) {
                    solutions.push_back(tableau[j][tableau[0].size() - 1]);
                    found = true;
                }
            }
            if (!found)
                solutions.push_back(0.000);
        }

        for (auto sol : solutions)
            cout << sol << " ";
        cout << endl;

        for (int i = 0; i < numConstraints; i++)
            cout << tableau[0][i] << " ";
        cout << endl;
    }

    void solveFeasible(vector<vector<double>>& tableau) {
        vector<vector<double>> secondTableau = secondPhaseTableau(tableau);

        simplex(secondTableau, false);

        for (int i = 0; i < secondTableau[0].size() - 1; i++) {
            if (secondTableau[0][i] < 0) {
                if (isUnbounded(secondTableau, i)) {
                    printUnbounded(secondTableau, i);
                    return;
                }
            }
        }

        printOptimal(secondTableau);
    }

    void run() {
        makePositive();

        vector<vector<double>> auxTableau = auxiliaryProblem();
        vector<vector<double>> tableau = buildTableau(auxTableau);

        cout << fixed << setprecision(3);

        int result = simplex(tableau, true);
        if (result == 0) {
            printInfeasible(tableau);
        } else if (result == 1) {
            solveFeasible(tableau);
        }
    }
};

int main(int argc, char* argv[]) {
    if (argc != 2) {
        cerr << "Usage: " << argv[0] << " <input_file>" << endl;
        return 1;
    }

    ifstream inputFile(argv[1]);
    if (!inputFile) {
        cerr << "Error opening file." << endl;
        return 1;
    }

    int n, m;
    double input;

    inputFile >> n >> m;

    vector<double> costs(m);
    for (int i = 0; i < m; i++) {
        inputFile >> input;
        costs[i] = input;
    }

    vector<vector<double>> constraints(n, vector<double>(m + 1));
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < m + 1; j++) {
            inputFile >> input;
            constraints[i][j] = input;
        }
    }

    inputFile.close();

    SimplexSolver solver(n, m, costs, constraints);
    solver.run();

    return 0;
}
