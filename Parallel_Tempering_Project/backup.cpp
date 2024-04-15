#include <iostream>
#include <vector>
#include <cmath>
#include <ctime>
#include <numeric>
#include <random>

// container's structures
struct _METROPOLIS_CONTAINER
{
    std::vector<int> STATE;
    double CAPACITY;
    double ENERGY;
};

struct _TEMPERING_CONTAINER
{
    std::vector<std::vector<double>> CAPACITY;
    std::vector<std::vector<double>> ENERGY;
    std::vector<std::vector<int>> SUCCESSFUL_COPIES;
};

class _SPIN_SYSTEM
{
private:
    int SIZE, NUM_OF_COPIES, REPEAT, NUM_OF_STEPS;
    float T_MIN, T_MAX;
    std::vector<std::vector<int>> CAPACITY;
    std::vector<std::vector<int>> ENERGY;
public:
    //  Class constructor
    _SPIN_SYSTEM(int _SIZE, int _NUM_OF_COPIES, int _REPEAT, int _NUM_OF_STEPS, float _T_MIN, float _T_MAX)
        : SIZE(_SIZE),
        NUM_OF_COPIES(_NUM_OF_COPIES),
        REPEAT(_REPEAT),
        NUM_OF_STEPS(_NUM_OF_STEPS),
        T_MIN(_T_MIN),
        T_MAX(_T_MAX) {}

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
    std::vector<std::vector<float>> temperatures_generate()
    {
        std::vector<std::vector<float>> temperatures(NUM_OF_STEPS, std::vector<float>(NUM_OF_COPIES));
        float delta_T = (T_MAX - T_MIN) / (NUM_OF_COPIES - 1);
        float plus_delta_T;

        for (int i = 0; i < NUM_OF_STEPS; ++i)
        {
            plus_delta_T = T_MIN;
            for (int j = 0; j < NUM_OF_COPIES; ++j)
            {
                temperatures[i][j] = plus_delta_T;
                plus_delta_T += delta_T;
            }
        }
        return temperatures;
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

    _METROPOLIS_CONTAINER metropolis(float T, std::vector<int> state, int REPEAT, int SIZE)
    {
        std::vector<int> E_VECT(REPEAT);

        std::mt19937 gen ((int)time(0));
        std::uniform_real_distribution<> urd (0., 1.);

        for (int j = 0; j < REPEAT; ++j)
        {
            int k = rand() % (SIZE * SIZE);
            std::vector<int> copied_state = state;

            copied_state[k] = (-1) * copied_state[k];

            int E1 = energy_metropolis(state);
            int E2 = energy_metropolis(copied_state);
            int delta_E = E2 - E1;

            double one = urd(gen);
            double p = std::pow(2.718282, (-1) * delta_E / T);


            if (one < p)
            {
                state = copied_state;
                // if (delta_E > 0) std::cout << "ERROR";
            }
            else
                delta_E = 0;

            if (j == 0)
                E_VECT[j] = energy_metropolis(state);
            else
                E_VECT[j] = E_VECT[j - 1] - delta_E;

            E_VECT[j] = energy_metropolis(state);
        }
        // mean energy
        int E_SUM = 0;
        for (int x : E_VECT) E_SUM += x;


        double MEAN_E = (double)E_SUM / (double)E_VECT.size();

        // mean squared_energy
        E_SUM = 0;
        for (int x : E_VECT) E_SUM += x * x;


        double MEAN_E_SQUARED = (double)E_SUM / (double)E_VECT.size();


        double capacity = (MEAN_E_SQUARED - MEAN_E * MEAN_E) / (T * T * SIZE * SIZE);
        _METROPOLIS_CONTAINER _container;

        _container.STATE = state;
        _container.ENERGY = MEAN_E;
        _container.CAPACITY = capacity;

        return _container;
    };

    _TEMPERING_CONTAINER tempering(std::vector<std::vector<float>> temperatures, std::vector<std::vector<int>> states)
    {
        std::vector<std::vector<int>> successful_copies(NUM_OF_STEPS, std::vector<int>(NUM_OF_COPIES - 1));
        std::vector<std::vector<double>> capacity(NUM_OF_STEPS, std::vector<double>(NUM_OF_COPIES));
        std::vector<std::vector<double>> energy(NUM_OF_STEPS, std::vector<double>(NUM_OF_COPIES));

        std::mt19937 gen ((int)time(0));
        std::uniform_real_distribution<> urd (0., 1.);

        for (int i = 0; i < NUM_OF_STEPS; ++i)
        {
            for (int j = 0; j < NUM_OF_COPIES - 1; ++j)
                successful_copies[i][j] = 0;
        }


        for (int m = 0; m < NUM_OF_STEPS; ++m)
        {
            for (int i = 0; i < NUM_OF_COPIES; ++i)
            {
                _METROPOLIS_CONTAINER metropolis_container = metropolis(temperatures[m][i], states[i], REPEAT, SIZE);

                states[i] = metropolis_container.STATE;
                capacity[m][i] = metropolis_container.CAPACITY;
                energy[m][i] = metropolis_container.ENERGY;
            }

            for (int i = NUM_OF_COPIES - 2; i > -1; i--)
            {
                double p = std::pow(2.718282, (energy[m][i + 1] - energy[m][i])*(1 /(double)temperatures[m][i + 1] - 1 /(double)temperatures[m][i]));

                if (urd(gen) < p)
                {
                    std::vector<int> tmp = states[i + 1];
                    states[i + 1] = states[i];
                    states[i] = tmp;

                    successful_copies[m][i] += 1;
                }
            }
        }
        _TEMPERING_CONTAINER tempering_container;

        tempering_container.CAPACITY = capacity;
        tempering_container.ENERGY = energy;
        tempering_container.SUCCESSFUL_COPIES = successful_copies;

        return tempering_container;
    };

};

float probability(std::vector<std::vector<int>> probabilities, int num_of_copies, int NUM_OF_STEPS)
{
    int num = num_of_copies-1;
    int probabilities_size = 0;
    std::vector<float> probabilities_1d (num*num);

    float sum;

    for (int i = 0; i < NUM_OF_STEPS; i++)
        probabilities_size += probabilities[i].size();

    for (int i = 0; i < num*num; i++)
    {
        sum = 0;

        for (int j = 0; j < NUM_OF_STEPS; j++)
            sum += probabilities[j][i];
        probabilities_1d[i] = (float)sum/(float)NUM_OF_STEPS;
    }

    sum = 0;

    for (float x : probabilities_1d)
        sum += x;
    


    return ((float)sum/(float)(num*num))*100;

};

int main()
{

    int num_of_copies = 2;
    int num_of_steps = 8;
    int size = 4;

    _SPIN_SYSTEM spin_system{size, num_of_copies, 1500 * (int)std::pow(4,2), num_of_steps, 0.1, 1.2 };
    //std::cout << 1500*(int)std::pow(4,2);
    std::vector<std::vector<int>> initial_states = spin_system.random_states();
    std::vector<std::vector<float>> temperatures = spin_system.temperatures_generate();

    std::cout << "TEST" << "\n";

    _TEMPERING_CONTAINER tempering_container = spin_system.tempering(temperatures, initial_states);

    std::vector<std::vector<double>> capacity = tempering_container.CAPACITY;
    std::vector<std::vector<double>> energy = tempering_container.ENERGY;
    std::vector<std::vector<int>> successful_copies = tempering_container.SUCCESSFUL_COPIES;

    for (double x : temperatures[0])
        std::cout << x << " ";
    std::cout << " \n";

    for (int i = 0; i < num_of_steps; ++i)
    {
        for (double x : capacity[i])
            std::cout << x << " ";
        std::cout << " \n";
    }

    for (int i = 0; i < num_of_steps; ++i)
    {
        for (double x : energy[i])
            std::cout << x << " ";
        std::cout << " \n";
    }

    for (int i = 0; i < num_of_steps; ++i)
    {
        for (int x : successful_copies[i])
            std::cout << x << " ";
        std::cout << " \n";
    }

    std::cout << probability(successful_copies, num_of_copies, num_of_steps) << " \n";

    return 0;
}