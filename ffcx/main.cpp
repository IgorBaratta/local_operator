#include "problem.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char *argv[]) {

  MPI_Init(&argc, &argv);
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // coordinate of a single cell
    std::array<scalar_type, 24> coords = {
        0.1, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.1, 1.0, 1.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0};

    // Read input
    constexpr int ndofs = dim;
    constexpr int ncoeffs = 2;
    constexpr int rank = kernel_rank;
    constexpr int local_size = rank == 1 ? dim : dim * dim;
    constexpr int stride = ndofs * ncoeffs;
    constexpr int ncells = global_size / ndofs;

    std::vector<scalar_type> geometry(ncells * 24);

    for (int cell = 0; cell < ncells; cell++) {
      auto it = std::next(geometry.begin(), 24 * cell);
      std::copy_n(coords.begin(), 24, it);
    }

    // Allocate and initialize data
    std::vector<scalar_type> A(ncells * local_size);
    std::vector<scalar_type> Ae(local_size);
    std::vector<scalar_type> coefficients(ncells * stride);
    std::fill(coefficients.begin(), coefficients.end(), 1.0);

    auto start = std::chrono::steady_clock::now();
    for (int cell = 0; cell < ncells; cell++) {
      std::fill(Ae.begin(), Ae.end(), 0);
      scalar_type *coeffs = coefficients.data() + cell * stride;
      auto geo = geometry.data() + 24 * cell;
      kernel(Ae.data(), coeffs, nullptr, geo, 0, 0);
      auto result = std::next(A.begin(), cell * local_size);
      std::copy(Ae.begin(), Ae.end(), result);
    }
    auto end = std::chrono::steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    double local_time = t.count();
    double max_time = 0;
    MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (mpi_rank == 0)
      std::cout << ncells << ", " << max_time / 1e6;
  }
  MPI_Finalize();

  return 0;
}