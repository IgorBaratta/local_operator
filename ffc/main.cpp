#include "data.hpp"
#include "problem.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <ufc.h>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

int main(int argc, char *argv[]) {

  problem_cell_integral_0_otherwise a;
  problem_cell_integral_1_otherwise L;

  // Read input
  constexpr int ndofs = dim;
  constexpr int ncoeffs = 2;
  constexpr int rank = kernel_rank;
  int ncells = rank == 1 ? 500000 : 100000;


  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::size_t local_size = rank == 1 ? ndofs : ndofs * ndofs;

  // Allocate and initialize data
  std::vector<double> A(ncells * local_size);
  std::vector<double> Ae(local_size);
  std::size_t stride = ndofs * ncoeffs;
  std::vector<double> coefficients(ncells * stride);
  std::uniform_real_distribution<double> dis(0, 1);
  std::generate(coefficients.begin(), coefficients.end(),
                [&]() { return dis(gen); });

  if (rank == 2) {
    // FFC
    auto start = std::chrono::steady_clock::now();
    for (int cell = 0; cell < ncells; cell++) {
      std::fill(Ae.begin(), Ae.end(), 0);
      double *coeffs = coefficients.data() + cell * stride;
      a.tabulate_tensor(Ae.data(), &coeffs, coordinate_dofs, 0);
      auto result = std::next(A.begin(), cell * local_size);
      std::copy(Ae.begin(), Ae.end(), result);
    }
    auto end = std::chrono::steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << ncells << ", " << t.count() / 1e6;
  } else {

    double **coeffs = new double *[2];
    coeffs[0] = new double[ndofs];
    coeffs[1] = new double[ndofs];

    auto start = std::chrono::steady_clock::now();
    for (int cell = 0; cell < ncells; cell++) {
      std::fill(Ae.begin(), Ae.end(), 0);
      coeffs[0] = coefficients.data() + cell * stride;
      coeffs[1] = coeffs[0] + ndofs;
      L.tabulate_tensor(Ae.data(), coeffs, coordinate_dofs, 0);
      auto result = std::next(A.begin(), cell * local_size);
      std::copy(Ae.begin(), Ae.end(), result);
    }
    auto end = std::chrono::steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    std::cout << ncells << ", " << t.count() / 1e6;
  }

  return 0;
}