#include <iostream>
#include <string.h>
#include <ctime>
#include <vector>

class demo
{
    private:
        uint32_t m_var_1 = 0;
        bool m_var_2 = false;
        std::string m_var_3 = "";
        float m_var_4 = 0.0;

    public:
        demo(uint32_t var_1, bool var_2, std::string var_3, float var_4)
            : m_var_1(var_1),
            m_var_2(var_2),
            m_var_3(var_3),
            m_var_4(var_4) {}
};

demo obj{123, true, "lol", 1.1};

int main()
{
    //srand(time(0));
   // for (int i = 0; i < 10; ++i) std::cout << rand()%2;

    int NUM_OF_COPIES = 2;
    float T_MAX = 1.2;
    float T_MIN = 0.1;

    std::vector<float> temperatures (NUM_OF_COPIES);
    float delta_T = (T_MAX - T_MIN)/NUM_OF_COPIES;
    std::cout << delta_T;
    float plus_delta_T = delta_T;
    std:: cout << "\n";
    for (int i=0; i < NUM_OF_COPIES; ++i)
    {   
        temperatures[i] = T_MIN + plus_delta_T;
        std::cout << temperatures[i] << " "; 
        plus_delta_T += delta_T; 
    }
    
    for (float x : temperatures)
        std::cout << x << " ";
    
    int temp_sum = -476300;
    int size_a = 24000;

    double E =  (double)(temp_sum) / (double)(size_a);
    std::cout << E << "\n";

    std::cout << 5 % 4;

    return 0;
}