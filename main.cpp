#include "problem.h"
#include <iostream>
#include <ufc.h>
#include <chrono>
#include <vector>
#include <algorithm>

int main(int argc, char *argv[])
{
  int ncells = 100'000;

  ufc_form a = *form_problem_a;
  int ndofs = a.finite_elements[1]->space_dimension;
  auto kernel = a.integrals(ufc_integral_type::cell)[0]->tabulate_tensor_float64;

  int rank = a.rank;

  const double coordinate_dofs[12] = {0.1, 0.0, 0.1, 1.0, 0.0, 0.1, 0.0, 1.0,
                                      0.0, 0.0, 0.0, 1.0};

  std::string type = "ffcx";
  if (argc == 2)
    type = argv[1];

  // Allocate data
  int local_size = (rank == 2) ? ndofs * ndofs : ndofs;
  std::vector<double> A(ncells * local_size);
  std::vector<double> Ae(local_size);
  std::vector<double> coefficients(ndofs + 1);
  std::fill(coefficients.begin(), coefficients.end(), 0.5);

  // FFCx
  auto start = std::chrono::steady_clock::now();
  for (int c = 0; c < ncells; c++)
  {
    std::size_t offset = c * local_size;
    std::fill(Ae.begin(), Ae.end(), 0);
    kernel(Ae.data(), coefficients.data(), nullptr, coordinate_dofs, 0, 0);
    auto result = std::next(A.begin(), offset);
    std::copy(Ae.begin(), Ae.end(), result);
  }
  auto end = std::chrono::steady_clock::now();
  double duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start).count() / 1.e6;
  std::cout << ncells << ", " << duration;
  return 0;
}