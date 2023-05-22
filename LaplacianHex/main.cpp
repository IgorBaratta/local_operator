#include "types.hpp"
#include "geometry.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <cassert>
#include <any>

#ifndef PRECISION
#error Floating point precision not defined.
#endif

#ifndef BATCH_SIZE
#error Batch size not defined.
#endif

#ifndef DEGREE
#error Polynomial degree not defined.
#endif

#ifndef BLOCK_SIZE
#error Block size not defined.
#endif

// This block enables to compile the code with and without the likwid header in place
#ifdef LIKWID_PERFMON
#include <likwid-marker.h>
#else
#define LIKWID_MARKER_INIT
#define LIKWID_MARKER_THREADINIT
#define LIKWID_MARKER_REGISTER(regionTag)
#define LIKWID_MARKER_START(regionTag)
#define LIKWID_MARKER_STOP(regionTag)
#define LIKWID_MARKER_CLOSE
#endif

int main(int argc, char *argv[])
{

  constexpr int precision = PRECISION;
  constexpr int batch_size = BATCH_SIZE;
  constexpr int P = DEGREE;
  constexpr int cub_nq = (P + 3) * (P + 3) * (P + 3);
  constexpr int bs = BLOCK_SIZE > 0 ? BLOCK_SIZE : cub_nq;

  constexpr bool precompute = PRECOMPUTE;

  using T = VectorExtensions<precision, batch_size>::T;
  using S = VectorExtensions<precision, batch_size>::S;
  constexpr int global_size = precompute ? 5000000 : 10000000;

  MPI_Init(&argc, &argv);
  {

    // Make sure that sizeof(T) == batch_size * sizeof(S)
    if (sizeof(T) != batch_size * sizeof(S))
    {
      std::string error_msg = "Size of T should be " + std::to_string(batch_size * sizeof(S));
      error_msg += " but current size is " + std::to_string(sizeof(T));
      throw std::runtime_error(error_msg);
    }

    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    // Create finite element operator
    Operator<T, S, P, bs, precompute> op;

    // Const data from kernel
    constexpr int local_size = op.num_dofs;
    constexpr int num_nodes = 8; // 1st order hexahedron
    constexpr int stride = op.num_dofs + num_nodes;
    constexpr int num_cells = global_size / op.num_dofs;
    constexpr int num_batches = num_cells / batch_size;
    constexpr int geom_size = precompute ? cub_nq * 6 : num_nodes * 3;

    // Allocate and initialize data
    std::vector<T> A(num_batches * local_size);

    // Constants for cross element vectorization
    T zero = {0};

    // Create geometry and coefficients
    std::vector<T> geometry = create_geometry<T, S>(num_batches, batch_size, cub_nq, precompute);
    std::vector<T> coefficients(num_batches * stride);

    std::size_t count = {0};
    std::for_each(coefficients.begin(), coefficients.end(), [&count, nd = op.num_dofs](auto &e)
                  { e = count < num_nodes? T(1) : T(S(count-num_nodes));
                    count = (count + 1) % (nd + num_nodes); });

    std::array<T, op.num_dofs> Ae;

    MPI_Barrier(comm);

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_REGISTER("kernel");

    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("kernel");

    double start = MPI_Wtime();
    for (int batch = 0; batch < num_batches; batch++)
    {
      std::fill(Ae.begin(), Ae.end(), zero);
      T *coeffs = coefficients.data() + batch * stride;
      T *geo = geometry.data() + batch * geom_size;
      op.apply(Ae.data(), coeffs, geo);
      std::vector<T>::iterator result = std::next(A.begin(), batch * local_size);
      std::transform(Ae.begin(), Ae.end(), result, result, std::plus<T>());
    }
    double end = MPI_Wtime();
    double local_time = end - start;

    LIKWID_MARKER_STOP("kernel");
    LIKWID_MARKER_CLOSE;

    MPI_Barrier(comm);

    double max_time = 0;
    double min_time = 0;
    MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);
    MPI_Allreduce(&local_time, &min_time, 1, MPI_DOUBLE, MPI_MIN, comm);

    if (mpi_rank == 0)
    {
      std::cout << precision << ", " << batch_size << ", ";
      std::cout << num_cells << ", " << P << ", " << max_time << ", ";
      std::cout << min_time << ", " << bs << ", " << precompute << ", ";
      std::cout << OPTIMIZE_SUM_FACTORIZATION;
    }
  }
  MPI_Finalize();

  return 0;
}