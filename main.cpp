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
    int ncells = 1'000'000;

    ufc_form a = *form_problem_a;
    int ndofs = a.finite_elements[1]->space_dimension;
    ufc_tabulate_tensor *kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

    const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                        0.0, 0.0, 0.0, 1.0};

    std::size_t type = 0;
    if (argc == 2)
        type = std::stoi(argv[1]);

    std::vector<double> A(ncells * ndofs * (ndofs + 1));

    std::vector<double, xsimd::aligned_allocator<double, 256>> coefficients(ndofs + 1);
    std::fill(coefficients.begin(), coefficients.end(), 0.5);

    if (type == 0)
    {
        // FFCx
        auto start = std::chrono::steady_clock::now();
        std::vector<double, xsimd::aligned_allocator<double, 256>> Ae(ndofs * ndofs);
        for (int c = 0; c < ncells; c++)
        {
            std::size_t offset = c * ndofs * ndofs;
            std::fill(Ae.begin(), Ae.end(), 0);
            kernel(Ae.data(), coefficients.data(), nullptr, coordinate_dofs, 0, 0);
            auto result = std::next(A.begin(), offset);
            std::copy(Ae.begin(), Ae.end(), result);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    else
    {
        // TSFC
        auto start = std::chrono::steady_clock::now();
        std::vector<double> Ae(ndofs * ndofs);
        for (int c = 0; c < ncells; c++)
        {
            std::size_t offset = c * ndofs * ndofs;
            std::fill(Ae.begin(), Ae.end(), 0);
            form_cell_integral_otherwise(Ae.data(), coordinate_dofs);
            auto result = std::next(A.begin(), offset);
            std::copy(Ae.begin(), Ae.end(), result);
        }
        auto end = std::chrono::steady_clock::now();
        double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
        std::cout << ncells << ", " << duration;
    }
    return 0;
}