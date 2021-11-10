#include "problem.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <numeric>
#include <random>
#include <vector>

std::random_device rd;
std::mt19937 gen(rd());

int main(int argc, char *argv[]) {

  // Read input
  int ncells = 1'000'000;
  constexpr int ndofs = dim;
  constexpr int ncoeffs = 2;
  constexpr int rank = kernel_rank;

  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::size_t local_size = rank == 1 ? dim : dim * dim;

  // Allocate and initialize data
  std::vector<double> A(ncells * local_size);
  std::vector<double> Ae(local_size);
  std::size_t stride = ndofs * ncoeffs;
  std::vector<double> coefficients(ncells * stride);
  std::uniform_real_distribution<double> dis(0, 1);
  std::generate(coefficients.begin(), coefficients.end(),
                [&]() { return dis(gen); });

  // tsfc
  auto start = std::chrono::steady_clock::now();
  for (int cell = 0; cell < ncells; cell++) {
    std::fill(Ae.begin(), Ae.end(), 0);
    double *coeffs = coefficients.data() + cell * stride;
#if MF == 1
    form_cell_integral_otherwise(Ae.data(), coordinate_dofs, coeffs,
                                 coeffs + ndofs);
#else
    form_cell_integral_otherwise(Ae.data(), coordinate_dofs, coeffs);
#endif
    auto result = std::next(A.begin(), cell * local_size);
    std::copy(Ae.begin(), Ae.end(), result);
  }
  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}