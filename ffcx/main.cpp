#include "problem.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {

  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  // Read input
  constexpr int ndofs = dim;
  constexpr int ncoeffs = 2;
  constexpr int rank = kernel_rank;
  int ncells = rank == 1 ? 500000 : 100000;

  std::size_t local_size = rank == 1 ? dim : dim * dim;

  // Allocate and initialize data
  std::vector<double> A(ncells * local_size);
  std::vector<double> Ae(local_size);
  std::size_t stride = ndofs * ncoeffs;
  std::vector<double> coefficients(ncells * stride, 1.0);

  // FFCx
  auto start = std::chrono::steady_clock::now();
  for (int cell = 0; cell < ncells; cell++) {
    std::fill(Ae.begin(), Ae.end(), 0);
    double *coeffs = coefficients.data() + cell * stride;
    kernel(Ae.data(), coeffs, nullptr, coordinate_dofs, 0, 0);
    auto result = std::next(A.begin(), cell * local_size);
    std::copy(Ae.begin(), Ae.end(), result);
  }
  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}