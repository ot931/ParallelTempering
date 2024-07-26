#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <numeric>
#include <random>
#include <fstream>
#include <time.h>
#include <algorithm>
#include <cstdlib>

// container's structure

float probability();

struct _TEMPERING_CONTAINER
{
    std::vector<double> CAPACITY;
    std::vector<double> ENERGY;
    std::vector<std::vector<int>> SUCCESSFUL_COPIES;
    std::vector<float> TEMPERATURES;
};

class _SPIN_SYSTEM
{
private:
    int SIZE, NUM_OF_COPIES, REPEAT, NUM_OF_STEPS, TN_STEPS, ACCURACY;
    float T_MIN, T_MAX;
    std::vector<std::vector<int>> CAPACITY;
    std::vector<std::vector<int>> ENERGY;
public:
    //  Class constructor
    _SPIN_SYSTEM(int _SIZE, int _NUM_OF_COPIES, int _REPEAT, int _NUM_OF_STEPS, int _TN_STEPS, int _ACCURACY, float _T_MIN, float _T_MAX)
        : SIZE(_SIZE),
        NUM_OF_COPIES(_NUM_OF_COPIES),
        REPEAT(_REPEAT),
        NUM_OF_STEPS(_NUM_OF_STEPS),
        TN_STEPS(_TN_STEPS),
        ACCURACY(_ACCURACY),
        T_MIN(_T_MIN),
        T_MAX(_T_MAX)
    {}

    // functions of class

    // generate the random state of spins
    std::vector<std::vector<int>> random_states()
    {
        std::vector<int> ones = { -1, 1 };
        std::vector<std::vector<int>> states(NUM_OF_COPIES, std::vector<int>(SIZE * SIZE));

        srand(time(0));

        for (int i = 0; i < NUM_OF_COPIES; ++i)
        {
            for (int j = 0; j < SIZE * SIZE; ++j)
                states[i][j] = ones[rand() % 2];
        }

        return states;
    };

    // generate the vector with temperatures 
    std::vector<float> temperatures_generate()
    {
        std::vector<float> temperatures(NUM_OF_COPIES);
        float delta_T = (T_MAX - T_MIN) / (NUM_OF_COPIES - 1);
        float plus_delta_T = T_MIN;

        for (int i = 0; i < NUM_OF_COPIES; ++i)
        {

            temperatures[i] = plus_delta_T;
            plus_delta_T += delta_T;
        }

        return temperatures;
    };

    // Make probability of exchange to 20%
    float energy_approx(float T1, float T2, double E1, double E2, float delta)
    {
        return ((T2 + delta - T1) / (T2 - T1)) * (E2 - E1) + E1;
    };

    float approximation(std::vector<double>& energies, std::vector<float>& temperatures, int tn_steps, double e_0, double e_1, float t_0, float t_1, float alpha, int i)
    {
        double en = e_1;
        float delta = 0;
        int stop = 0;

        for (int k = 0; k < tn_steps; k++)
        {
            double p = std::pow(2.718282, (((double)e_0 - (double)en) * ((1 / (double)t_0) - 1 / (double)(t_1 + delta))));

            if (p > 0.2)
            {
                delta += alpha;
                en = energy_approx(t_0, t_1, e_0, e_1, delta);
            }

            if (p < 0.2 || t_0 > temperatures[i+1])
            {
                if (stop == 30) break;

                delta -= alpha;
                en = energy_approx(t_0, t_1, e_0, e_1, delta);

                stop += 1;
            }
        }

        temperatures[i + 1] = t_1 + delta;
        energies[i + 1] = en;

        return 0;
    };

    std::vector<float> temperature_normalize(std::vector<float> temperatures, std::vector<double> energies, float alpha, int tn_steps, int m)
    {
        std::vector<float> temperatures_0 = temperatures;
        std::vector<double> energies_0 = energies;
        float t_0 = 0; float t_1 = 0;
        double e_0 = 0; double e_1 = 0;

        for (int i = 0; i < NUM_OF_COPIES - 1; i++)
        {
            t_0 = temperatures[i];  t_1 = temperatures_0[i + 1];
            e_0 = energies[i];      e_1 = energies_0[i + 1];

            approximation(energies, temperatures, tn_steps, e_0, e_1, t_0, t_1, alpha, i);

            //std::sort(energies.begin(), energies.end());
        }


        std::cout << "energies: ";
        for (int i = 0; i < NUM_OF_COPIES; i++)
            std::cout << energies[i] << " ";
        std::cout << std::endl;

        std::ofstream out;

        out.open("Temperatures_bf.txt", std::fstream::app);
        if (out.is_open())
        {
            for (double x : temperatures_0)
                out << x << " ";
            out << " \n";
        }
        out.close();

        out.open("Temperatures_af.txt", std::fstream::app);
        if (out.is_open())
        {
            for (double x : temperatures)
                out << x << " ";
            out << " \n";
        }
        out.close();


        out.open("Energies_bf.txt", std::fstream::app);
        if (out.is_open())
        {
            for (double x : energies_0)
                out << x << " ";
            out << " \n";
        }
        out.close();

        out.open("Energies_af.txt", std::fstream::app);
        if (out.is_open())
        {
            for (double x : energies)
                out << x << " ";
            out << " \n";
        }
        out.close();

        std::sort(temperatures.begin(), temperatures.end());
        std::sort(energies.begin(), energies.end());

        return temperatures;
    };

    double energy_mean(std::vector<double> e_arr, int REPEAT)
    {
        double e_sum = 0;
        for (int i = 0; i < REPEAT; i++) e_sum += e_arr[i];

        return e_sum / (double)REPEAT;
    };

    int energy_metropolis(std::vector<int> state)
    {
        std::vector<std::vector<int>> energy_array(SIZE, std::vector<int>(SIZE));

        for (int i = 0; i < SIZE; ++i)
        {
            for (int j = 0; j < SIZE; ++j)
            {
                if (i + 1 == SIZE && j + 1 == SIZE)
                    energy_array[i][j] = (-1) * state[i * SIZE + j] * (state[i * SIZE] + state[j]);
                else if (j + 1 == SIZE)
                    energy_array[i][j] = (-1) * state[i * SIZE + j] * (state[i * SIZE] + state[(i + 1) * SIZE + j]);
                else if (i + 1 == SIZE)
                    energy_array[i][j] = (-1) * state[i * SIZE + j] * (state[i * SIZE + (j + 1)] + state[j]);
                else
                    energy_array[i][j] = (-1) * state[i * SIZE + j] * (state[i * SIZE + (j + 1)] + state[(i + 1) * SIZE + j]);
            }
        }

        int sum = 0;

        for (int i = 0; i < SIZE * SIZE; ++i)
            sum += energy_array[i / SIZE][i % SIZE];

        return sum;
    };


    _TEMPERING_CONTAINER tempering(std::vector<float> temperatures, std::vector<std::vector<int>> states)
    {
        std::vector<std::vector<int>> successful_copies(NUM_OF_STEPS, std::vector<int>(NUM_OF_COPIES - 1));

        std::vector<double> capacity(NUM_OF_COPIES);
        std::vector<double> energy(NUM_OF_COPIES);

        std::vector<std::vector<double>> E_ARRAY(NUM_OF_COPIES, std::vector<double>(REPEAT));
        std::vector<std::vector<double>> E_ARRAY_SQUARED(NUM_OF_COPIES, std::vector<double>(REPEAT));

        std::vector<double> E_MEAN_SQUARED(NUM_OF_COPIES);
        std::vector<double> E_MEAN(NUM_OF_COPIES);

        std::mt19937 gen((int)time(0));
        std::uniform_real_distribution<> urd(0., 1.);


        // feelling zeroes in parametres
        for (int i = 0; i < NUM_OF_COPIES; ++i)
        {
            E_MEAN_SQUARED[i] = 0;
            E_MEAN[i] = 0;
            for (int j = 0; j < REPEAT; ++j)
            {
                E_ARRAY[i][j] = 0;
                E_ARRAY_SQUARED[i][j] = 0;
            }
        }

        for (int i = 0; i < NUM_OF_STEPS; ++i)
        {
            for (int j = 0; j < NUM_OF_COPIES - 1; ++j)
                successful_copies[i][j] = 0;
        }

        for (int m = 0; m < NUM_OF_STEPS; ++m)
        {
            for (int i = 0; i < NUM_OF_COPIES; ++i)
            {

                for (int j = 0; j < REPEAT; ++j)
                {
                    int k = rand() % (SIZE * SIZE);

                    std::vector<int> copied_state = states[i];

                    copied_state[k] = (-1) * copied_state[k];

                    int E1 = energy_metropolis(states[i]);
                    int E2 = energy_metropolis(copied_state);
                    int delta_E = E2 - E1;

                    double one = urd(gen);
                    double p = std::pow(2.718282, (-1) * delta_E / temperatures[i]);


                    if (one < p)
                    {
                        states[i] = copied_state;
                        // if (delta_E > 0) std::cout << "ERROR";
                    }

                    E_ARRAY[i][j] = energy_metropolis(states[i]);
                    E_ARRAY_SQUARED[i][j] = std::pow(E_ARRAY[i][j], 2);
                }
            }

            std::vector<double> E_MEAN_1(NUM_OF_COPIES);

            if (m > ACCURACY)
            {
                for (int i = 0; i < NUM_OF_COPIES; i++)
                {
                    E_MEAN[i] += energy_mean(E_ARRAY[i], REPEAT);
                    E_MEAN_SQUARED[i] += energy_mean(E_ARRAY_SQUARED[i], REPEAT);
                }

                for (int i = 0; i < NUM_OF_COPIES; i++)
                    E_MEAN_1[i] = E_MEAN[i] / (m - ACCURACY);

                temperatures = temperature_normalize(temperatures, E_MEAN_1, 0.45, TN_STEPS, m);

                for (int i = 0; i < NUM_OF_COPIES; i++)
                    std::cout << temperatures[i] << " ";
                std::cout << std::endl;

            }
            for (int i = NUM_OF_COPIES - 2; i > -1; i--)
            {
                double p = std::pow(2.718282, (((double)energy_mean(E_ARRAY[i + 1], REPEAT) / (double)(m + 1) - (double)energy_mean(E_ARRAY[i], REPEAT) / (double)(m + 1)) * (1 / (double)temperatures[i + 1] - 1 / (double)temperatures[i])));

                if (urd(gen) < p)
                {
                    std::vector<int> tmp = states[i + 1];
                    states[i + 1] = states[i];
                    states[i] = tmp;

                    successful_copies[m][i] += 1;
                }
            }
        }

        for (int i = 0; i < NUM_OF_COPIES; i++)
        {
            E_MEAN[i] /= (NUM_OF_STEPS - (ACCURACY + 1));
            E_MEAN_SQUARED[i] /= (NUM_OF_STEPS - (ACCURACY + 1));

        }

        for (int i = 0; i < NUM_OF_COPIES; i++)
        {

            capacity[i] = (E_MEAN_SQUARED[i] - E_MEAN[i] * E_MEAN[i]) / (temperatures[i] * temperatures[i] * SIZE * SIZE);

        }
        _TEMPERING_CONTAINER tempering_container;

        tempering_container.CAPACITY = capacity;
        tempering_container.ENERGY = E_MEAN;
        tempering_container.SUCCESSFUL_COPIES = successful_copies;
        tempering_container.TEMPERATURES = temperatures;

        return tempering_container;
    };

};

float probability(std::vector<std::vector<int>> probabilities, int num_of_copies, int NUM_OF_STEPS)
{
    int num = num_of_copies - 1;
    int probabilities_size = 0;
    std::vector<float> probabilities_1d(num);

    float sum;

    for (int i = 0; i < NUM_OF_STEPS; i++)
        probabilities_size += probabilities[i].size();

    for (int i = 0; i < num; i++)
    {
        sum = 0;

        for (int j = 0; j < NUM_OF_STEPS; j++)
            sum += probabilities[j][i];
        probabilities_1d[i] = (float)sum / (float)NUM_OF_STEPS;
    }

    sum = 0;

    for (float x : probabilities_1d)
        sum += x;



    return ((float)sum / (float)(num)) * 100;

};



int main()
{
    int num_of_copies = 6;
    int num_of_steps = 20;
    int tn_steps = 10000;
    int size = 4;
    float percent = 0.1;
    int accuracy = num_of_steps * percent;


    _SPIN_SYSTEM spin_system{ size, num_of_copies, 1500 * (int)std::pow(4,2), num_of_steps, tn_steps, accuracy, 0.1, 5 };
    //std::cout << 1500*(int)std::pow(4,2);
    std::vector<std::vector<int>> initial_states = spin_system.random_states();
    std::vector<float> temperatures = spin_system.temperatures_generate();

    std::cout << "START:" << "\n";

    clock_t tStart = clock();

    _TEMPERING_CONTAINER tempering_container = spin_system.tempering(temperatures, initial_states);

    std::vector<double> capacity = tempering_container.CAPACITY;
    std::vector<double> energy = tempering_container.ENERGY;
    std::vector<std::vector<int>> successful_copies = tempering_container.SUCCESSFUL_COPIES;

    temperatures = tempering_container.TEMPERATURES;


    std::cout << "Executing time: " << (double)(clock() - tStart) / CLOCKS_PER_SEC << " sec. \n";

    for (double x : temperatures)
        std::cout << x << " ";
    std::cout << " \n";


    for (double x : capacity)
        std::cout << x << " ";
    std::cout << " \n";



    for (double x : energy)
        std::cout << x << " ";
    std::cout << " \n";


    std::cout << probability(successful_copies, num_of_copies, num_of_steps) << " %" << "\n";

    std::ofstream out;

    // Temperatures
    out.open("Temperatures.txt");

    if (out.is_open())
    {
        for (double x : temperatures)
            out << x << " ";
        out << " \n";
    }
    out.close();


    // Capacity
    out.open("Capacity.txt");

    if (out.is_open())
    {

        for (double x : capacity)
            out << x << " ";
        out << " \n";

    }
    out.close();


    // Energy
    out.open("Energy.txt");

    if (out.is_open())
    {

        for (double x : energy)
            out << x << " ";
        out << " \n";

    }
    out.close();

    return 0;
}
