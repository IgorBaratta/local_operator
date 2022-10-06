#include "problem.hpp"
#include "geometry.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <vector>

int main(int argc, char *argv[])
{

  MPI_Init(&argc, &argv);
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // Read input
    constexpr int num_dofs = dim;
    constexpr int local_size = kernel_rank == 1 ? dim : dim * dim;
    constexpr int stride = num_dofs * num_coefficients;
    constexpr int num_cells = global_size / num_dofs;
    constexpr int num_batches = num_cells / batch_size;
    constexpr int geom_size = num_nodes * 3;

    // Allocate and initialize data
    std::vector<scalar_type> A(num_batches * local_size);

    scalar_type one = {1.};
    std::vector<geom_type> geometry = create_geometry<geom_type>(num_batches, batch_size, geom_size);
    std::vector<scalar_type> coefficients(num_batches * stride);
    std::for_each(
        coefficients.begin(), coefficients.end(), [=](auto &e)
        { e = one; });

    std::vector<scalar_type> Ae(local_size);
    scalar_type zero = {0.};

    auto start = std::chrono::steady_clock::now();
    for (int batch = 0; batch < num_batches; batch++)
    {
      std::fill(Ae.begin(), Ae.end(), zero);
      scalar_type *coeffs = coefficients.data() + batch * stride;
      geom_type *geo = geometry.data() + geom_size * batch;
      kernel(Ae.data(), coeffs, nullptr, geo, 0, 0);
      scalar_type *result = A.data() + batch * local_size;
      std::copy_n(Ae.begin(), num_dofs, result);
    }
    auto end = std::chrono::steady_clock::now();
    auto t = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    double local_time = t.count();
    double max_time = 0;
    MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX,
                  MPI_COMM_WORLD);

    if (mpi_rank == 0)
      std::cout << num_cells << ", " << max_time / 1e6;
  }
  MPI_Finalize();

  return 0;
}