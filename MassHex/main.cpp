#include "types.hpp"
#include "problem.hpp"

#include <algorithm>
#include <array>
#include <chrono>
#include <iostream>
#include <mpi.h>
#include <vector>
#include <cassert>

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
#define PRECISION 0
#endif

#ifndef OPTIMIZE_SUM_FACTORIZATION
#define OPTIMIZE_SUM_FACTORIZATION 0
#endif

#ifndef BATCH_SIZE
#define BATCH_SIZE 1
#endif

#ifndef DEGREE
#define DEGREE 0
#endif

int main(int argc, char *argv[])
{

  using T = VectorExtensions<PRECISION, BATCH_SIZE>::T;
  using S = VectorExtensions<PRECISION, BATCH_SIZE>::S;
  constexpr int global_size = 50000000;

  MPI_Init(&argc, &argv);
  {
    if (sizeof(T) != BATCH_SIZE * sizeof(S))
    {
      std::string error_msg = "Size of T should be " + std::to_string(BATCH_SIZE * sizeof(S));
      error_msg += " but current size is " + std::to_string(sizeof(T));
      throw std::runtime_error(error_msg);
    }

    MPI_Comm comm = MPI_COMM_WORLD;
    int mpi_rank;
    MPI_Comm_rank(comm, &mpi_rank);

    constexpr int P = DEGREE;
    Operator<T, S, P> op;

    // Const data from kernel
    constexpr int local_size = op.num_dofs;
    constexpr int stride = op.num_dofs;
    constexpr int num_cells = global_size / op.num_dofs;
    constexpr int num_batches = num_cells / BATCH_SIZE;

    constexpr int geom_size = (P + 2) * (P + 2) * (P + 2);

    // Allocate and initialize data
    std::vector<T> A(num_batches * local_size);

    // Constants for cross element vectorization
    T one = {1};
    T zero = {0};

    std::vector<T> geometry(geom_size * num_batches);
    std::vector<T> coefficients(num_batches * stride);
    auto set_ = [one](auto &e)
    { e = one; };
    std::for_each(coefficients.begin(), coefficients.end(), set_);
    std::for_each(geometry.begin(), geometry.end(), set_);
    std::array<T, op.num_dofs> Ae;

    // Sanity check: Are we computing the correct values?
    // for (int batch = 0; batch < 100; batch++)
    // {
    //   std::array<T, op.num_dofs> Ae = {0};
    //   T *coeffs = coefficients.data() + batch * stride;
    //   T *geo = geometry.data() + batch * geom_size;
    //   op.apply(Ae.data(), coeffs, geo);

    //   // Compute area of cell times number of quadrature points
    //   // T acc = 0;
    //   // for (std::size_t i = 0; i < Ae.size(); i++)
    //   //   acc += Ae[i];

    //   // check_solution<T, S>(acc, reference);
    // }
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
      std::cout << PRECISION << ", " << BATCH_SIZE << ", ";
      std::cout << BATCH_SIZE * num_batches << ", " << DEGREE << ", " << max_time;
      std::cout << ", " << min_time << ", " << OPTIMIZE_SUM_FACTORIZATION;
    }
  }

  MPI_Finalize();

  return 0;
}