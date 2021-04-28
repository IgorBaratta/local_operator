// #define XTENSOR_USE_XSIMD 1
#include "lagrange.h"
#include <iostream>
#include <xtensor/xtensor.hpp>
#include <xtensor/xio.hpp>
#include <ufc.h>
#include <chrono>

int main(int argc, char *argv[])
{
    int ncells = 100'000;

    ufc_form a = *form_lagrange_a;
    int ndofs_cell = a.finite_elements[1]->space_dimension;
    auto &&kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor;

    xt::xtensor<double, 2> coordinate_dofs = {{0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0}, {0.0, 0.0, 1.0}};
    // xt::xtensor<double, 2> Ae = xt::empty<double>({ndofs_cell, ndofs_cell});

    double *Ae = new double[ndofs_cell * ndofs_cell];

    auto start = std::chrono::steady_clock::now();
    for (int c = 0; c < ncells; c++)
    {
        kernel(Ae, nullptr, nullptr, coordinate_dofs.data(), 0, 0, 0);
    }
    auto end = std::chrono::steady_clock::now();
    double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
    std::cout << ncells << ", " << duration;
    return 0;
}