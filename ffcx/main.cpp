#include "problem.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {

  const scalar_type coordinate_dofs[24] = {
      0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0,
      0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  // Read input
  constexpr int ndofs = dim;
  constexpr int ncoeffs = 2;
  constexpr int rank = kernel_rank;
  constexpr int local_size = rank == 1 ? dim : dim * dim;
  constexpr int stride = ndofs * ncoeffs;
  constexpr int ncells = global_size / ndofs;

  // Allocate and initialize data
  std::vector<scalar_type> A(ncells * local_size);
  std::vector<scalar_type> Ae(local_size);
  std::vector<scalar_type> coefficients(ncells * stride);
  std::fill(coefficients.begin(), coefficients.end(), 1.0);

  auto start = std::chrono::steady_clock::now();
  for (int cell = 0; cell < ncells; cell++) {
    std::fill(Ae.begin(), Ae.end(), 0);
    scalar_type *coeffs = coefficients.data() + cell * stride;
    kernel(Ae.data(), coeffs, nullptr, coordinate_dofs, 0, 0);
    auto result = std::next(A.begin(), cell * local_size);
    std::copy(Ae.begin(), Ae.end(), result);
  }
  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}