#include <iostream>
#include <fstream>
#include <cstdlib>
#include <ctime>
#include <random>

using namespace std;

const unsigned short int STATES = 5;
const unsigned short int ROW_METHOD = 5;
const unsigned short int MOVES = 100;
const unsigned short int EXPERIMENTS = 50;

void out_file_matrix(ofstream& file_name, const double& matrix, int size1, int size2, bool key_console_output)
{
    for (int i = 0; i < size1; i++)
    {
        for (int j = 0; j < size2; j++)
            file_name << *(&matrix + i * STATES + j) << "\t";
        file_name << "\n";
    }
    file_name.close();

    if (key_console_output)
        for (int i = 0; i < size1; i++)
        {
            for (int j = 0; j < size2; j++)
                cout << *(&matrix + i * STATES + j) << "\t";
            cout << "\n";
        }
}

void out_file_vector(ofstream& file_name, const double& vector, int size, bool key_console_output)
{
    for (int i = 0; i < size; i++)
        file_name << *(&vector + i) << "\t";
    file_name.close();

    if (key_console_output)
    {
        for (int i = 0; i < size; i++)
            cout << *(&vector + i) << "\t";
        cout << "\n";
    }
}

void show_matrix(const double& matrix)
{
    for (int i = 0; i < STATES; i++)
    {
        for (int j = 0; j < STATES; j++)
            cout << *(&matrix + i * STATES + j) << "\t";
        cout << "\n";
    }
}

void show_vector(const double& vector)
{
    for (int i = 0; i < STATES; i++)
        cout << *(&vector + i) << "\n";
}


float RNG_generator_float()
{
    static mt19937 rng(std::random_device{}());

    return generate_canonical< float, 128 >(rng);
}


int main()
{
    ifstream in_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_in_data.txt");

    if (!in_file.is_open())
    {
        cerr << "Error: Files for input or output could not be opened." << endl;
        exit(2);
    }
    double table0[STATES][STATES];

    double table[STATES][STATES];

    double moves[EXPERIMENTS][STATES] = {};

    for (int i = 0; i < STATES; i++)
        for (int j = 0; j < STATES; j++)
        {
            in_file >> table0[i][j];
            table[i][j] = table0[i][j];
        }
    in_file.close();
    show_matrix(**table);

    // Транспонируем
    cout << "\nTransposed Matrix:\n";
        double a;
        for (int i = 0; i < STATES - 1; i++)
            for (int j = i + 1; j < STATES; j++)
            {
                a = table[i][j];
                table[i][j] = table[j][i];
                table[j][i] = a;
            }
    show_matrix(**table);

    // Готовим матрицу для расчта вектора предельных вероятностей
    cout << "\npT - E Matrix\n";
    for (int i = 0; i < STATES; i++)
        table[i][i] -= 1;
    show_matrix(**table);

    // Модифицируем матрицу
    cout << "\nModified Matrix:\n";
    for (int i = 0; i < STATES; i++)
        table[STATES - 1][i] = 1;
    show_matrix(**table);

    double Y_table[STATES];
    for (int i = 0; i < STATES; i++)
        Y_table[i] = 0;
    Y_table[ROW_METHOD - 1] = 1;

    //Gauss Method
    cout << "\nGauss Method Matrix Conversion:\n";
    double mod;
    for (int i = 0; i < STATES; i++)
        for (int j = i + 1; j < STATES; j++)
        {
            if (fabs(table[i][i]) < 1e-8)
                continue;
            mod = table[j][i] / table[i][i];
            table[j][i] = 0;
            for (int k = i + 1; k < STATES; k++)
            {
                table[j][k] -= mod * table[i][k];
            }
        }
    show_matrix(**table);

    // Считаем итоговый вектор
    double res[STATES];
    for (int i = STATES - 1; i >= 0; i--)
    {
        res[i] = Y_table[i];
        for (int j = STATES - 1; j > i; j--)
            res[i] -= table[i][j] * res[j];
        res[i] = res[i] / table[i][i];
    }
    ofstream out_res_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_res_data.txt");
    cout << "\nResult Vector:\n";
    out_file_vector(out_res_file, *res, STATES, true);

    cout << "\nFinal Transition Matrix:\n";
    for (int i = 0; i < STATES; i++)
        for (int j = 0; j < STATES; j++)
            table[j][i] = res[i];
    show_matrix(**table);


    for (int i = 0; i < STATES; i++)
        for (int j = 0; j < STATES; j++)
        {
            table[i][j] = table0[i][j];
        }


    cout << "\nTable for RNG gen:\n";
    for (int i = 0; i < STATES; i++)
        for (int j = 1; j < STATES; j++)
            table[i][j] += table[i][j - 1];
    show_matrix(**table);

    {
        ofstream out_moves_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_moves_data.txt");
        double value;
        int state;

        for (int i = 0; i < EXPERIMENTS; i++)
        {
            state = static_cast<int>(RNG_generator_float() * 5 + 1);

            out_moves_file << state << "\t";
            for (int j = 1; j < MOVES; j++)
            {
                value = RNG_generator_float();
                for (int k = 0; k < STATES; k++)
                    if (value < table[state - 1][k])
                    {
                        state = k + 1;
                        moves[i][state - 1]++;
                        break;
                    }
                out_moves_file << state << "\t";
            }
            out_moves_file << "\n";
        }
        out_moves_file.close();
    }

    ofstream out_end_moves_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_end_moves_data.txt");
    out_file_matrix(out_end_moves_file, **moves, EXPERIMENTS, STATES, false);

    // Расчет вероятностей перехода
    for (int i = 0; i < EXPERIMENTS; i++)
        for (int j = 0; j < STATES; j++)
            moves[i][j] = moves[i][j] / MOVES;
    ofstream out_end_chances_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_end_chances_data.txt");
    cout << "\nEnd Chances Table:\n";
    out_file_matrix(out_end_chances_file, **moves, EXPERIMENTS, STATES, true);


    double sr_ar[STATES] = {};
    for (int j = 0; j < STATES; j++)
    {
        for (int i = 0; i < EXPERIMENTS; i++)
            sr_ar[j] += moves[i][j];
        sr_ar[j] = sr_ar[j] / EXPERIMENTS;
    }
    ofstream out_sr_ar_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_sr_ar_data.txt");
    cout << "\nSR-AR Vector:\n";
    out_file_vector(out_sr_ar_file, *sr_ar, STATES, true);

    double sr_kv[STATES] = {};
    for (int i = 0; i < EXPERIMENTS; i++)
        for (int j = 0; j < STATES; j++)
            sr_kv[j] += static_cast<double>(pow(moves[i][j] - sr_ar[j], 2));
    for (int i = 0; i < STATES; i++)
        sr_kv[i] = sqrt(sr_kv[i] / (EXPERIMENTS - 1));

    ofstream out_sr_kv_file("../Source/AMIMSU_NikolaychukDS_Task1_v108_out_sr_kv_data.txt");
    cout << "\nSR-KV Matrix:\n";
    out_file_vector(out_sr_kv_file, *sr_kv, STATES, true);


    return 0;
}