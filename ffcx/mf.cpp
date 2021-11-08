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

template <typename T>
using kernel_fn = std::function<void(T *, const T *, const T *, const T *,
                                     const int *, const std::uint8_t *)>;

int main(int argc, char *argv[]) {

  // Read input
  int ncells = 1'000'000;

  ufc_form L = *form_problem_L;
  int ndofs = L.finite_elements[0]->space_dimension;
  int ncoeffs = L.num_coefficients;

  const kernel_fn<double> &kernel =
      L.integrals(ufc_integral_type::cell)[0]->tabulate_tensor_float64;

  if (L.rank != 1)
    std::runtime_error("Form rank should be 1 for matrix free application");

  const double coordinate_dofs[24] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0, 0.1, 0.0, 0.1, 1.0,
                                      0.0, 0.1, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0};

  // Allocate and initialize data
  std::vector<double> A(ncells * ndofs);
  std::vector<double> Ae(ndofs + 1);
  std::vector<double> coefficients(ncells * (ndofs * ncoeffs));
  std::uniform_real_distribution<double> dis(0, 1);
  std::generate(coefficients.begin(), coefficients.end(),
                [&]() { return dis(gen); });

  // FFCx
  auto start = std::chrono::steady_clock::now();
  for (int c = 0; c < ncells; c++) {
    std::size_t offset = c * ndofs;
    std::fill(Ae.begin(), Ae.end(), 0);
    double *coeffs = coefficients.data() + offset;
    kernel(Ae.data(), coeffs, nullptr, coordinate_dofs, 0, 0);
    auto result = std::next(A.begin(), offset);
    std::copy(Ae.begin(), Ae.end(), result);
  }

  auto end = std::chrono::steady_clock::now();
  auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
  std::cout << ncells << ", " << t.count() / 1e6;

  return 0;
}