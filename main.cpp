#include "problem.h"
#include "tsfc_kernel.hpp"

#include <iostream>
#include <ufc.h>
#include <chrono>

int main(int argc, char *argv[])
{
    int ncells = 100'000;

    ufc_form a = *form_problem_a;
    int ndofs_cell = a.finite_elements[1]->space_dimension;
    auto &&kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

    const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0};

    std::size_t type = 0;
    if (argc == 2)
        type = std::stoi(argv[1]);

    if (type == 0)
    {
        double Ae[DOFS * DOFS];
        double coefficients[DOFS];

        for (int i = 0; i < DOFS; i++)
            coefficients[i] = 1.0;
        
        auto start = std::chrono::steady_clock::now();
        for (int c = 0; c < ncells; c++)
        {
            kernel(Ae, coefficients, nullptr, coordinate_dofs, 0, 0);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    else
    {
        double Ae[DOFS][DOFS];
        auto start = std::chrono::steady_clock::now();
        for (int c = 0; c < ncells; c++)
        {
            form_cell_integral_otherwise(Ae, coordinate_dofs);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    return 0;
}