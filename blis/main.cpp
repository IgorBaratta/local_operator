#include "types.hpp"

#include "blis/blis.h"

#include "problem.hpp"
#include "geometry.hpp"
#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <numeric>
#include <mpi.h>
#include <vector>
#include <cassert>
#include <any>

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

#ifndef PRECISION
#define PRECISION 8
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 1
#endif

#ifndef DEGREE
#define DEGREE 1
#endif

int main(int argc, char *argv[])
{

  using T = VectorExtensions<PRECISION, BATCH_SIZE>::S;
  using S = VectorExtensions<PRECISION, BATCH_SIZE>::S;
  constexpr int global_size = 50000000;

  MPI_Init(&argc, &argv);
  {
    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    constexpr int P = DEGREE;
    Operator<T, P, BATCH_SIZE> op;

    rntm_t rntm;
    bli_rntm_init(&rntm);
    bli_rntm_set_num_threads(1, &rntm);

    // Const data from kernel
    constexpr int local_size = op.num_dofs;
    constexpr int stride = op.num_dofs;
    constexpr int num_cells = global_size / op.num_dofs;
    constexpr int num_batches = num_cells / BATCH_SIZE;
    constexpr int geom_size = 4 * 3;

    // Allocate and initialize data
    std::vector<T> A(num_batches * BATCH_SIZE * local_size);

    // Constants for cross element vectorization
    T zero = {0};

    // Create geometry and coefficients
    std::vector<T> geometry = create_geometry<T, S>(num_batches * BATCH_SIZE, 1, geom_size);
    std::vector<T> coefficients(num_batches * BATCH_SIZE * stride);
    std::fill(coefficients.begin(), coefficients.end(), T(1.0));

    // Sanity check: Are we computing the correct values?
    std::array<T, local_size * BATCH_SIZE> Ae;

    MPI_Barrier(comm);

    LIKWID_MARKER_INIT;
    LIKWID_MARKER_REGISTER("kernel");

    LIKWID_MARKER_THREADINIT;
    LIKWID_MARKER_START("kernel");

    double start = MPI_Wtime();
    for (int batch = 0; batch < num_batches; batch++)
    {
      std::fill(Ae.begin(), Ae.end(), zero);
      T *coeffs = coefficients.data() + batch * stride * BATCH_SIZE;
      T *geo = geometry.data() + batch * geom_size * BATCH_SIZE;
      op.apply(Ae.data(), coeffs, geo);
      std::vector<T>::iterator result = std::next(A.begin(), batch * local_size * BATCH_SIZE);
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
      std::cout << PRECISION << ", " << BATCH_SIZE << ", " << num_cells << ", " << DEGREE << ", " << max_time;
      std::cout << ", " << min_time;
    }
  }
  MPI_Finalize();

  return 0;
}
