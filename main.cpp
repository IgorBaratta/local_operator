#include "problem.h"
#include "tsfc_kernel.cpp"
#include <iostream>
#include <ufc.h>
#include <chrono>
#include <vector>
#include <algorithm>
#include <xsimd/xsimd.hpp>

int main(int argc, char *argv[])
{
    int ncells = 100'000'000;

    ufc_form a = *form_problem_a;
    int ndofs_cell = a.finite_elements[1]->space_dimension;
    ufc_tabulate_tensor *kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

    const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0};

    std::size_t type = 0;
    if (argc == 2)
        type = std::stoi(argv[1]);

    constexpr std::size_t dofs = DOFS;
    std::vector<double> A(ncells * ndofs_cell * ndofs_cell);

    double coefficients[dofs];
    std::fill_n(coefficients, dofs, 1);

    if (type == 0)
    {
        auto start = std::chrono::steady_clock::now();
        std::vector<double, xsimd::aligned_allocator<double, 256>> Ae(dofs * dofs);
        for (int c = 0; c < ncells; c++)
        {
            std::size_t offset = c * dofs * dofs;
            std::fill(Ae.begin(), Ae.end(), 0);
            kernel(Ae.data(), nullptr, nullptr, coordinate_dofs, 0, 0);
            auto result = std::next(A.begin(), offset);
            std::copy(Ae.begin(), Ae.end(), result);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    else
    {
        auto start = std::chrono::steady_clock::now();
        for (int c = 0; c < ncells; c++)
        {
            double Ae[DOFS][DOFS] = {0};
            form_cell_integral_otherwise(Ae, coordinate_dofs);
            double *data = A.data() + c * DOFS * DOFS;
            std::copy_n(&Ae[0][0], dofs * dofs, data);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    return 0;
}