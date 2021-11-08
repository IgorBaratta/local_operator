#include "problem.h"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <ufc.h>
#include <vector>

template <typename T>
using kernel_fn = std::function<void(T *, const T *, const T *, const T *,
                                     const int *, const std::uint8_t *)>;

int main(int argc, char *argv[]) {

  // Read input
  int ncells = 1'000'000;

  ufc_form a = *form_problem_a;
  int ndofs = a.finite_elements[1]->space_dimension;

  const kernel_fn<double> &kernel =
      a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor_float64;

  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  // Allocate data
  int local_size = (ndofs) * (ndofs);
  std::vector<double> A(ncells * local_size);
  std::vector<double> Ae(local_size);
  std::vector<double> coefficients(ndofs);
  std::fill(coefficients.begin(), coefficients.end(), 0.5);

  // FFCx
  auto start = std::chrono::steady_clock::now();

  for (int c = 0; c < ncells; c++) {
    std::size_t offset = c * local_size;
    std::fill(Ae.begin(), Ae.end(), 0);
    kernel(Ae.data(), coefficients.data(), nullptr, coordinate_dofs, 0, 0);
    auto result = std::next(A.begin(), offset);
    std::copy(Ae.begin(), Ae.end(), result);
  }

  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}