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

  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  std::vector<scalar_type> coords(24);
  for (std::size_t i = 0; i < coords.size(); i++)
    for (int j = 0; j < batch_size; j++)
      coords[i][j] = coordinate_dofs[i];

  // Read input
  constexpr int ndofs = dim;
  constexpr int ncoeffs = 2;
  constexpr int rank = kernel_rank;
  int ncells = rank == 1 ? 500000 : 100000;

  std::size_t local_size = rank == 1 ? dim : dim * dim;

  // Allocate and initialize data
  std::vector<scalar_type> A(ncells/batch_size * local_size);
  std::vector<scalar_type> Ae(local_size);
  std::size_t stride = ndofs * ncoeffs;
  std::vector<scalar_type> coefficients(ncells/batch_size * stride);

  for (std::size_t i = 0; i < coefficients.size(); i++)
    for (int j = 0; j < batch_size; j++)
      coefficients[i][j] = 1.;

  scalar_type zero = {0.0};

  // FFCx
  auto start = std::chrono::steady_clock::now();
  for (int cell = 0; cell < ncells / batch_size; cell++) {
    std::fill(Ae.begin(), Ae.end(), zero);
    scalar_type *coeffs = coefficients.data() + cell * stride;
    kernel(Ae.data(), coeffs, nullptr, coords.data(), 0, 0);
    auto result = std::next(A.begin(), cell * local_size);
    std::copy(Ae.begin(), Ae.end(), result);
  }
  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}