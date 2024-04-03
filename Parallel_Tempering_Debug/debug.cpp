#include <iostream>
#include <vector>
#include <ctime>

#define SIZE 4

double r2()
{
    return (double)rand() / (double)RAND_MAX ;
}

std::vector<int> random_state()
{
    std::vector<int> ones = { -1, 1 };
    std::vector<int> state(SIZE * SIZE);

    srand(time(0));

    for (int i = 0; i < SIZE * SIZE; i++)
        state[i] = ones[rand() % 2];

    return state;
};


int energy_metropolis(std::vector<int> state)
{
    std::vector<std::vector<int>> energy_array(SIZE, std::vector<int>(SIZE));

    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
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

    for (int i = 0; i < SIZE * SIZE; i++)
    
        sum += energy_array[i / SIZE][i % SIZE];

    return sum;
};

float random_to_one()
{

    float random = 0.1 * (rand() % 11);

    if (random == 0) return 0.001; 
    else return random;

};


int main()
{
    std::vector<int> state = {1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1};
    int energy = energy_metropolis(state);

    for (int x : state)
        std::cout << x << " ";
    std::cout << " \n";


    std::cout << energy  << "\n";



    srand(time(0));
    std::cout << random_to_one() <<"\n";
    std::cout << random_to_one() <<"\n";
    std::cout << random_to_one() <<"\n";
    std::cout << r2();


    return 0;
}
