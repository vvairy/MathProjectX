#include <vector>
#include <algorithm>
#include <cmath>
#include <map>
#include <string>
#include<set>
#include<array>
#include<iostream>
using std::exception;
using std::cin;
using std::cout;
using std::vector;
using std::map;
using std::string;
using std::endl;
using std::pair;
using std::set;
using std::array;
class Matrix
{
public:
    Matrix(vector <vector<double>> m)
    {
        if (!isMatrixRectangle(m))
            throw std::exception("Not Matrix is Given");
        this->matrix = m;
    }
    Matrix()
    {
    }
    void print()
    {
        for (size_t i = 0; i < size(this->matrix); ++i)
        {
            for (size_t j = 0; j < size(this->matrix[0]); ++j)
                std::cout << this->matrix[i][j] << " ";
            std::cout << std::endl;
        }
    }
    vector<vector<double>> matrix;
    double getDeterminantOfMatrix()
    {
        try { return determinant(this->matrix); }
        catch (std::exception& ex)
        {
            std::cout << ex.what() << std::endl;
        }
    }

    vector <double> solve_SLAU_by_Kramer(vector <vector<double>>& Koefficients, vector <double>& Answers)
    {
        vector<double> determinants;
        for (size_t i = 0; i < Koefficients.size(); i++)
            determinants.push_back(getEveryDeterminant(Koefficients, Answers, i) / determinant(Koefficients));
        return determinants;
    }

    void transpositione()
    {
        for (size_t i = 0; i < this->matrix.size(); i++)
            for (size_t j = i; j < this->matrix[0].size(); j++)
                std::swap(this->matrix[i][j], this->matrix[j][i]);
    }
    
    int rank()
    {
        Matrix m(this->matrix);
        int zeroCount = size(m.matrix) - 2;
        for (size_t j = 0; j < size(m.matrix) - 2; j++)
            for (size_t i = 0; i < zeroCount; i++)
            {
                columnMeltipliedSubtraction(matrix[i], matrix[i + 1], -matrix[i][j] / matrix[i + 1][j]);
                --zeroCount;
            }
        return countNonZeroRows(m);
    }
protected:
    const double PI = 3.14159;
    const double E = 2.71828;
    double determinant(const Matrix& m)
    {
        if (!isMatrixSquare(matrix))
            throw std::exception("Not square Matrix");
        if (size(matrix) == 2)
            return matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
        double det = 0;
        for (size_t idx = 0; idx < size(matrix); idx++)
        {
            vector<vector<double>> NewMatrix(size(matrix) - 1);
            for (size_t i = 1; i < size(matrix); i++)
                for (size_t j = 0; j < size(matrix); j++)
                    if (j != idx)
                        NewMatrix[i - 1].push_back(matrix[i][j]);
            det += pow(-1, idx) * matrix[0][idx] * determinant(NewMatrix);
            NewMatrix.clear();
        }
        return det;
    }
    friend Matrix operator* (const Matrix& M1, const Matrix& M2)
    {
        if (size(M1.matrix) != size(M2.matrix[0]) or size(M2.matrix) != size(M1.matrix[0]))
            throw std::exception("Non-multiplyable matrixes");
        size_t min_len = std::min(size(M1.matrix), size(M1.matrix[0]));
        size_t max_len = std::max(size(M1.matrix), size(M1.matrix[0]));

        vector<vector<double>> NewMatrix(min_len);
        for (size_t i = 0; i < min_len; i++)
            for (size_t j = 0; j < min_len; j++)
            {
                double el = 0;
                for (size_t idx = 0; idx < max_len; idx++)
                    el += M1.matrix[i][idx] * M2.matrix[idx][j];
                NewMatrix[i].push_back(el);
            }
        return NewMatrix;
    }

private:
    template <typename T>
    bool isMatrixSquare(const T& Matrix)
    {
        return size(Matrix) == size(Matrix[0]) and isMatrixRectangle(Matrix);
    }
    template <typename T>
    bool isMatrixRectangle(T Matrix)
    {
        for (size_t i = 1; i < size(Matrix); i++)
            if (Matrix[0].size() != Matrix[i].size())
                return false;
        return true;
    }

    vector<double> columnMeltipliedSubtraction(vector <double>& subtrahend, vector <double>& subtractor, double multiplier)
    {
        if (size(subtrahend) != size(subtractor))
            throw std::exception("vector sizes are not equal");
        for (size_t i = 0; i < size(subtractor); ++i)
            subtrahend[i] -= subtractor[i] * multiplier;
        return subtrahend;
    }
    bool isRowZero(vector<double>& row)
    {
        for (size_t j = 0; j < size(row); ++j)
            if (row[j] != 0)
                return false;
        return true;
    }
    int countNonZeroRows(Matrix& m)
    {
        int count = 0;
        for (size_t i = 0; i < size(m.matrix); ++i)
            if (!isRowZero(matrix[i]))
                count++;
        return count;
    }

    double getEveryDeterminant(vector <vector<double>>& Matrix, vector <double>& Answers, size_t strikethroughСol)
    {
        vector<vector<double>> NewMatrix(size(Matrix));
        for (size_t i = 0; i < size(Matrix); i++)
            for (size_t j = 0; j < size(Matrix[0]); j++)
                if (j != strikethroughСol)
                    NewMatrix[i].push_back(Matrix[i][j]);
                else
                    NewMatrix[i].push_back(Answers[i]);
        return determinant(NewMatrix);
    }
};

class MathVector : Matrix
{
public:
    struct Dot
    {
        array <double, 3> dotCoordinates;
        Dot(double x, double y, double z)
        {
            this->dotCoordinates = { x, y, z };
        }
    };
    MathVector(double x, double y, double z)
    {
        this->coordinates = { x, y, z };
    }
    MathVector(Dot a, Dot b)
    {
        this->coordinates = { b.dotCoordinates[0] - a.dotCoordinates[0], b.dotCoordinates[1] - a.dotCoordinates[1], b.dotCoordinates[2] - a.dotCoordinates[2] };
    }
    double scalarMultiplication(MathVector& A, MathVector& B)
    {
        double res = 0;
        for (size_t i = 0; i < 3; i++)
            res += A.coordinates[i] * B.coordinates[i];
        return res;
    }
    double moduleOfVector(MathVector A)
    {
        double sum = 0;
        for (double el : A.coordinates)
            sum += el * el;
        return sqrt(sum);
    }
    double cornerBetweenVectors(MathVector& A, MathVector& B)
    {
        return acos(scalarMultiplication(A, B) / (moduleOfVector(A) * moduleOfVector(B)));
    }
 
   
    /*double mixedVectorMultiplication(MathVector firstVector, MathVector secondvector)
    {
        vector<vector<double>> matrixOfVectors = {this->coordinates, firstVector.coordinates, secondvector.coordinates};// надо пофиксить тему с конструктором
        return determinant(matrixOfVectors);
    }*/
    /*double distanceBetweenIntersectingLines(Dot& DotA, MathVector& guideVectorA, Dot& DotB, MathVector& guideVectorB)
    {
        MathVector distanceBetweenDots = MathVector(DotA, DotB);
        array<MathVector, 3> volume = {};
        volume[0] = distanceBetweenDots; volume[1] = (guideVectorA); volume[2] = (guideVectorB);
        return abs(mixedVectorMultiplication(volume) / vectorMultiplication(guideVectorA, guideVectorB));
    }*/

    double vectorMultiplication(const MathVector& A, const MathVector& B)
    {
        array<array<double, 3>, 2> matrix = { A.coordinates, B.coordinates };
        MathVector C = {0, 0, 0};
        for (size_t idx = 0; idx < 3; idx++)
        {
            vector<vector<double>> newMatrix;
            for (size_t i = 0; i < 2; i++)
                for (size_t j = 0; j < 3; j++)
                    if (j != idx)
                        newMatrix[i][j] = (matrix[i][j]);
            C.coordinates[idx] = (determinant(newMatrix));
            newMatrix.clear();
        }
        return moduleOfVector(C);
    }
    /*double volumeOfParallelepipedOnVectors(array<array<double, 3>, 3>& mathVectors)
    {
        return abs(mixedVectorMultiplication(mathVectors));
    }*/

    /*double volumeOfPiramidOnVectors(array<array<double, 3>, 3>& mathVectors)
    {
        return volumeOfParallelepipedOnVectors(mathVectors) / 6;
    }*/
private:
    array <double, 3> coordinates;
    array<array<double, 3>, 3> getVectorsFromDots(array<array<double, 3>, 3>& coordinatesOfDots)
    {
        array<array<double, 3>, 3> mathVectors;
        for (size_t i = 0; i < 3; i++)
            for (size_t j = 0; j < 3; j++)
                mathVectors[i][j] = coordinatesOfDots[i + 1][j] - coordinatesOfDots[0][j];
        return mathVectors;
    }
};
void getCanonicEquation(vector<double> Dot, vector<double> guideVector)
{
    std::cout << "(x - " << Dot[0] << ")/" << guideVector[0] << " = (y - " << Dot[1] << ")/" << guideVector[1] << //чет говно какое-то
        " = (z - " << Dot[2] << ")/" << guideVector[2] << " = 0";
}
class ProbabilityTheory
{
protected:
    //ExpectedValue - с англ математическое ожидание
    virtual double getExpectedValue(array<vector<double>, 2>& componentTable) //выборочное среднее составляющЕЙ (1) - аналог математического ожидания для выборки
    {
        double expectedValue = 0;
        for (size_t i = 0; i < size(componentTable[0]); i++)
            expectedValue += componentTable[0][i] * componentTable[1][i];
        return expectedValue;
    }
    double getTwoDimensionalExpectedValue(vector<vector<double>>& correlationTable, int nums_count) //выборочное среднее составляющИХ (2) - аналог математического ожидания двумерной величины для выборки
    {
        double expectedValue = 0;
        for (size_t i = 1; i < size(correlationTable); i++)
            for (size_t j = 1; j < size(correlationTable[0]); j++)
                expectedValue += correlationTable[i][0] * correlationTable[i][j] * correlationTable[0][j];
        return expectedValue / nums_count;
    }
    virtual double getDispersion(array<vector<double>, 2>& componentTable) //дисперсия - величина, показывающая среднее отклонение случайной величины от матожидания
    {
        double Dispersion = 0;
        for (size_t i = 0; i < size(componentTable[0]); i++)
            Dispersion += pow(componentTable[0][i], 2) * componentTable[1][i];
        Dispersion -= pow(getExpectedValue(componentTable), 2); // D[x] = M[x^2] - M^2[x]
        return Dispersion;
    }

    double localLaplassFunction(double arg)
    {
        double value;
        value = (1 / sqrt(PI) * pow(EXP, -1 * pow(arg, 2) / 2));
        return value;
    }

    unsigned NumberOfCombinations(unsigned NumKits, unsigned NumElements)
    {
        return Factorial(NumKits) / Factorial(NumElements) / Factorial(NumKits - NumElements);
    }

    unsigned NumberOfCombinationsRepetitions(unsigned NumKits, unsigned NumElements)
    {
        return Factorial(NumKits) / Factorial(NumKits - NumElements);
    }

    unsigned Factorial(unsigned Number)
    {
        if (Number == 1 or Number == 0)
            return 1;
        return Number * Factorial(Number - 1);
    }

    double getCountInVec(double number, vector<double>& vec)
    {
        double count = 0;
        for (size_t idx = 0; idx < vec.size(); idx++)
            if (vec[idx] == number)
                count++;
        return count;
    }

    const double PI = 3.14159;
    const double EXP = 2.71828;
};
class BinomialDistribution : ProbabilityTheory
{
    BinomialDistribution(double probabiliy, int countOfExperiments)
    {
        this->p = probabiliy;
        this->q = 1 - probabiliy;
        this->n = countOfExperiments;
    }
    double getExpectedValue(array<vector<double>, 2>& componentTable) override
    {
        return this->n * this->p;
    }
    double getDispersion(array <vector<double>, 2>& componentTable) override
    {
        return this->n * this->p * this->q;
    }
    array<vector<double>, 2> fillComponentTable() // p - общепринятое обозначение вероятности (Probability)
    {
        array<vector<double>, 2> componentTable;
        for (int k = 0; k <= this->n; k++)
        {
            componentTable[0].push_back(k);
            componentTable[1].push_back(localLaplassFunction(((k - this->n * this->p) / sqrt(this->n * this->p * this->q))));
            //локальная функция Лапласа принимает отношение отклонения от матождиания к среднему квадр отклонению ((m-M[x])/sqrt(D[x]))
        }
    }
private:
    double n, p, q;
};

class NormalDistribution : protected ProbabilityTheory
{
    //...
};
class PirsonCheck : ProbabilityTheory
{
public:
    PirsonCheck(vector<double>& input)
    {
        this->sample_ = input;
        this->nums_count = size(input);
        this->rows["validated"] = this->equidistant(input);
        this->h = this->SturgessFormula(input);
    }

    double SturgessFormula(vector<double> input)
    {
        return (this->rows["validated"].back() - this->rows["validated"][0]) / (1 + 3.32 * log10(size(input)));//находим шаг равноотстоящей выборки
    }
    bool isPirson()
    {
        cout << "Do some magic: \n h = " << this->h << endl;
        this->rows["counters"] = this->getCounters();
        this->rows["counters_in_range"] = this->getCountersInRange();
        this->rows["summed_with_h"] = this->getSummedWithH();
        this->row_len = size(getSummedWithH());
        this->rows["divided_counters"] = this->getDividedCounters();
        this->selectionTable[0] = this->rows["summed_with_h"]; this->selectionTable[1] = this->rows["counters_in_range"];

        double expectedValue = getExpectedValue(selectionTable);
        this->dispersion = (getDispersion(selectionTable) + pow(expectedValue, 2) - pow(expectedValue, 2) / pow(this->nums_count, 2));
        expectedValue /= this->nums_count;
        double standardDeviation = sqrt(this->dispersion);// 

        double s_2 = (double(this->nums_count) / double(this->nums_count - 1)) * dispersion;
        double s = sqrt(s_2);

        this->rows["laplass_arguments"] = this->getLaplassArguments(expectedValue, standardDeviation);
        this->rows["laplass_results"] = this->getLaplassResults();
        this->rows["theoretical_frequency"] = this->getTheoreticalFrequency(standardDeviation);

        double x_2_q = this->get_x_2_q();
        int degreesOfFreedom = this->getDegreesOfFreedom();

        return x_2_q - degreesOfFreedom * 1.4 >= 2.4; //условие соответсвия нормальному Закону Распредeления
    }
    void getResults()
    {
        for (auto it = rows.begin(); it != rows.end(); ++it)
        {
            cout << it->first << endl;
            for (size_t i = 0; i < (it->second.size()); ++i)
            {
                cout << " " << it->second[i];
            }
            cout << endl;
        }
        cout << "Supporting Value = " << get_x_2_q() << "\n Count degrees of Freedoms = " << getDegreesOfFreedom()
            << "\n expectedValue = " << getExpectedValue(selectionTable) / this->nums_count << "\n Dispersion = " << this->dispersion
            << "\n Average Square Deviation = " << sqrt(this->dispersion);
    }
private:
    vector <double> sample_;
    map<string, vector<double>> rows;
    double h; // обозначение расстояния между равноотстоящими вариантами
    size_t nums_count;
    size_t row_len;
    double dispersion;
    array <vector<double>, 2> selectionTable;
    vector<double> getSummedWithH()
    {
        vector<double> summedWithH;
        double xMin = this->rows["validated"][0], xMax = this->rows["validated"].back();
        while (xMin < xMax)
        {
            summedWithH.push_back(xMin);
            xMin += this->h;
        }
        summedWithH.push_back(xMin);
        return summedWithH;
    }


    vector<double> getCounters()
    {
        vector <double> counters;
        for (size_t idx = 0; idx < this->rows["validated"].size(); idx++)
        {
            counters.push_back(getCountInVec(this->rows["validated"][idx], this->sample_));
        }
        return counters;
    }

    vector<double> getCountersInRange()
    {
        double xMin = this->rows["validated"][0], xMax = this->rows["validated"].back();
        vector<double> counters;
        while (xMin < xMax)
        {
            counters.push_back(this->countInRange(xMin));
            xMin += this->h;
        }
        counters.push_back(this->countInRange(xMin));
        return counters;
    }

    vector<double> getDividedCounters()
    {
        vector<double> dividedCounters(this->row_len);
        for (int idx = 0; idx < this->row_len; idx++)
            dividedCounters[idx] = this->rows["counters_in_range"][idx] / this->nums_count;
        return dividedCounters;
    }

    vector<double> getLaplassArguments(double expectedValue, double standartDeviation)
    {
        vector<double> laplassArguments(this->row_len);
        for (size_t idx = 0; idx < this->row_len; idx++)
            laplassArguments[idx] = (this->rows["summed_with_h"][idx] - expectedValue) / standartDeviation;
        return laplassArguments;
    }

    vector<double> getLaplassResults()
    {
        vector<double> laplassResults(this->row_len);
        for (size_t idx = 0; idx < this->row_len; idx++)
            laplassResults[idx] = localLaplassFunction(this->rows["laplass_arguments"][idx]);
        return laplassResults;
    }

    vector<double> getTheoreticalFrequency(double standartDeviation)
    {
        vector<double> theoreticalFrequency;
        for (size_t idx = 0; idx < this->row_len; idx++)
            theoreticalFrequency.push_back((this->nums_count * this->h) / standartDeviation * this->rows["laplass_results"][idx]);
        return theoreticalFrequency;
    }

    double get_x_2_q()
    {
        double x_2_q = 0;
        for (size_t idx = 0; idx < this->row_len; idx++)
            x_2_q += pow(this->rows["counters_in_range"][idx] - this->rows["theoretical_frequency"][idx], 2) / this->rows["theoretical_frequency"][idx];
        return x_2_q;
    }

    int getDegreesOfFreedom()
    {
        return this->row_len - 3;
    }

    int countInRange(double x1)
    {
        int counter = 0;
        for (int i = 0; i < size(this->sample_); i++)
        {
            if (x1 - this->h / 2 <= this->sample_[i] && this->sample_[i] < x1 + this->h / 2)
                counter++;
        }
        return counter;
    }

    vector<double> equidistant(vector<double>& vec)
    {
        vector<double> res;
        sort(vec.begin(), vec.end());
        for (size_t i = 1; i < size(vec); i++)
            if (vec[i] != vec[i - 1])
                res.push_back(vec[i - 1]);
        res.push_back(vec.back());
        return res;
    }
};

vector<double> setSample()
{
    vector <double> input_data;
    int c;
    cout << "Enter nums count: ";
    cin >> c;
    cout << "Enter nums:\n";
    double number;
    int i = 1;
    while (i <= c)
    {
        cout << i << ". ";
        cin >> number;
        input_data.push_back(number);
        ++i;
    }
    return input_data;
}
class Correlation : ProbabilityTheory
{
public:
    Correlation(array<vector<double>, 2> pairs)
    {
        this->pairs = pairs;
        this->experimentsCount_ = size(pairs[0]);
        this->componentTableX = getComponentTable(pairs[0]);
        this->componentTableY = getComponentTable(pairs[1]);
        this->correlationTable = getCorrelationTable(this->componentTableX, this->componentTableY);
    }
    map <std::string, double> numberCharacteristics;

private:
    size_t experimentsCount_;
    array<vector<double>, 2> pairs;
    vector<vector <double>> correlationTable;
    array <vector <double>, 2> componentTableX;
    array <vector <double>, 2> componentTableY;
    double getExpectedValue(array<vector<double>, 2>& componentTable) override
    {
        return ProbabilityTheory::getExpectedValue(componentTable) / this->experimentsCount_;
    }
    double getDispersion(array <vector<double>, 2>& componentTable) override
    {
        double Dispersion = 0;
        for (size_t i = 0; i < size(componentTable[0]); i++)
            Dispersion += pow(componentTable[0][i], 2) * componentTable[1][i];
        Dispersion -= pow(this->getExpectedValue(componentTable), 2); //D[x] = M[x^2] - M^2{x]
        return Dispersion / this->experimentsCount_;
    }
    void fillNumberCharacteristics()
    {
        this->numberCharacteristics["expectedValueX"] = this->getExpectedValue(this->componentTableX);
        this->numberCharacteristics["expectedValueY"] = this->getExpectedValue(this->componentTableY);
        this->numberCharacteristics["dispersionX"] = getDispersion(this->componentTableX);
        this->numberCharacteristics["dispersionY"] = getDispersion(this->componentTableY);
        this->numberCharacteristics["expectedValueXandY"] = getTwoDimensionalExpectedValue(this->correlationTable, this->experimentsCount_);
        this->numberCharacteristics["correlationMoment"] = getCorrelationMoment(this->numberCharacteristics["expectedValueX"], this->numberCharacteristics["expectedValueY"], this->numberCharacteristics["expectedValueXandY"]);
        this->numberCharacteristics["correlationKoefficient"] = getCorrelationKoefficient(this->numberCharacteristics["correlationMoment"], this->numberCharacteristics["dispersionX"], this->numberCharacteristics["dispersionY"]);
    }
    array <vector<double>, 2> getComponentTable(vector<double>& row)
    {
        vector <double> experements = getUniqueAndSort(row);
        vector <double> countEachExperement;
        for (double el : experements)
            countEachExperement.push_back(getCountInVec(el, row));
        return { countEachExperement, experements };
    }
    vector <double> getUniqueAndSort(vector<double>& row)
    {
        vector <double> componentTable;
        std::sort(row.begin(), row.end());
        for (size_t i = 1; i < this->experimentsCount_; ++i)
            if (row[i] != row[i - 1])
                componentTable.push_back(row[i - 1]);
        componentTable.push_back(row.back());
        return componentTable;
    }

    unsigned countOfSamePairs(array <vector <double>, 2>& pairs, double x, double y)
    {
        unsigned count = 0;
        for (size_t i = 0; i < 4; i++)
            if (pairs[0][i] == x and pairs[1][i] == y)
                count++;
        return count;
    }

    vector <vector <double>> getCorrelationTable(array <vector<double>, 2>& componentTableX, array <vector<double>, 2>& componentTableY)
    {
        vector <vector <double>> correlationTable; correlationTable[0].push_back(NULL);
        correlationTable.push_back(componentTableX[0]);
        for (size_t i = 0; i < size(componentTableY[0]); i++)
            correlationTable[i + 1].push_back(componentTableY[0][i]);
        for (size_t i = 1; i <= size(componentTableY[1]); i++)
            for (size_t j = 1; j <= size(componentTableX[1]); j++)
                correlationTable[i].push_back(countOfSamePairs(this->pairs, componentTableX[1][j - 1], componentTableY[1][i - 1]));
        return correlationTable;
    }

    double getCorrelationMoment(double expectedValueX, double  expectedValueY, double expectedValueXandY)
    {
        return expectedValueXandY - expectedValueX * expectedValueY;
    }

    bool isCorrelate(double correlationMoment)
    {
        return (correlationMoment);
    }

    double getCorrelationKoefficient(double correlationMoment, double dispersionX, double dispersionY)
    {
        return sqrt(pow(correlationMoment, 2) / (dispersionX * dispersionY));
    }

    void strenghtOfCorrelation(double  correlationKoefficient)
    {
        if (correlationKoefficient > 0.7)
        {
            cout << " correlation is strong\n"; makeLineRegression();
        }
        else if (correlationKoefficient >= 0.3)
        {
            cout << " correlation is middle\n";  makeLineRegression();
        }
        else
            cout << " correlation is weak \n";
    }
    vector<double> makeLineRegression()
    {
        //...
    }
};



int main()
{
    /*setlocale(LC_ALL, "russian");
    vector <double> input_data = setSample();
    PirsonCheck pirsonCheck(input_data);
    pirsonCheck.isPirson();
    pirsonCheck.getResults();*/
    Matrix m({ {1, 2, 3,}, {1, 2, 11} , {1,2,-2} });
    std::cout << m.getDeterminantOfMatrix();
    Matrix n({ { 1,4,3 }, { 10, 2, 1 }, { 1,4,5 } });
    Matrix c = m * n; //overloaded operator * to multiply matrixes;
    c.print(); // print matrix to concole
    c.transpositione();
    c.print(); // print matrix to concole
    std::cout << c.rank();
}
/* {
    vector<vector<double>> fillMatrix(size_t row_len, size_t col_len)
    {
        vector<vector<double>> matrix(row_len);
        for (size_t i = 0; i < row_len; i++)
           for (size_t j = 0; j < col_len; j++)
            {
                double el;
                std::cout << "Введите " << i + 1 << " " << j + 1 << " ";
                std::cin >> el;
                matrix[i].push_back(el);
            }
        return matrix;
    }
   }*/
