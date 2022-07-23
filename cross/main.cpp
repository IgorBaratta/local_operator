#include "problem.hpp"
#include <algorithm>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <numeric>
#include <vector>

int main(int argc, char *argv[]) {
  MPI_Init(&argc, &argv);
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    constexpr int ndofs = dim;
    constexpr int ncoeffs = 2;
    constexpr int rank = kernel_rank;
    constexpr int local_size = rank == 1 ? dim : dim * dim;
    constexpr int stride = ndofs * ncoeffs;
    constexpr int ncells = global_size / ndofs;
    constexpr int num_batches = ncells / batch_size;

    std::array<double, 24> coordinate_dofs = {
        0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 1.0, 0.0, 0.0, 1.0, 1.0,
        1.0, 0.0, 0.0, 1.0, 0.0, 1.0, 1.0, 1.0, 0.0, 1.0, 1.0, 1.0};

    std::vector<scalar_type> coords(24 * num_batches);
    for (std::size_t c = 0; c < num_batches; c++)
      for (std::size_t i = 0; i < 24; i++)
        for (int j = 0; j < batch_size; j++)
          coords[c * 24 + i][j] = coordinate_dofs[i];

    // Allocate and initialize data
    std::vector<scalar_type> A(num_batches * local_size);
    std::vector<scalar_type> Ae(local_size);
    std::vector<scalar_type> coefficients(num_batches * stride);

    for (std::size_t i = 0; i < coefficients.size(); i++)
      for (int j = 0; j < batch_size; j++)
        coefficients[i][j] = 1.;

    scalar_type zero = {0.0};

    // FFCx
    auto start = std::chrono::steady_clock::now();
    for (int cell = 0; cell < num_batches; cell++) {
      std::fill(Ae.begin(), Ae.end(), zero);
      scalar_type *coeffs = coefficients.data() + cell * stride;
      scalar_type *geo = coords.data() + 24 * cell;
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