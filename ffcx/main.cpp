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

    // Const data from kernel
    constexpr int num_dofs = dim;
    constexpr int local_size = kernel_rank == 1 ? dim : dim * dim;
    constexpr int stride = num_dofs * num_coefficients;
    constexpr int num_cells = global_size / num_dofs;
    constexpr int num_batches = num_cells / batch_size;
    constexpr int geom_size = num_nodes * 3;

    // Allocate and initialize data
    std::vector<scalar_type> A(num_batches * local_size);

    // Constants for cross element vectorization
    scalar_type one = {1.};
    scalar_type zero = {0.};

    // Create geometry and coefficients
    std::vector<geom_type> geometry = create_geometry<geom_type>(num_batches, batch_size, geom_size);
    std::vector<scalar_type> coefficients(num_batches * stride);
    auto set_ = [one](auto &e)
    { e = one; };
    std::for_each(coefficients.begin(), coefficients.end(), set_);

    std::array<scalar_type, local_size> Ae;

    double start = MPI_Wtime();
    for (int batch = 0; batch < num_batches; batch++)
    {
      std::fill(Ae.begin(), Ae.end(), zero);
      scalar_type *coeffs = coefficients.data() + batch * stride;
      geom_type *geo = geometry.data() + batch * geom_size;
      kernel(Ae.data(), coeffs, nullptr, geo, 0, 0);
      std::vector<scalar_type>::iterator result = std::next(A.begin(), batch * local_size);
      std::transform(Ae.begin(), Ae.end(), result, result, std::plus<scalar_type>());
    }
    double end = MPI_Wtime();
    double local_time = end - start;

    double max_time = 0;
    MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);

    if (mpi_rank == 0)
      std::cout << num_cells << ", " << max_time;
  }
  MPI_Finalize();

  return 0;
}