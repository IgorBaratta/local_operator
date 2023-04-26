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

template <typename T, typename S>
void check_solution(T &acc, T &reference)
{
  if constexpr (std::is_same<T, S>::value)
  {
    if ((acc - reference) * (acc - reference) > T(0.001))
    {
      throw std::runtime_error("Please verify solution.");
    }
  }
  else
  {
    if ((acc[0] - reference[0]) * (acc[0] - reference[0]) > S(0.001))
    {
      throw std::runtime_error("Please verify solution.");
    }
  }
}

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
    T reference = {0};

    // Create geometry and coefficients
    std::vector<T> geometry = create_geometry<T, S>(num_batches, batch_size, cub_nq, precompute);

    std::vector<T> coefficients(num_batches * stride);
    // std::fill(coefficients.begin(), coefficients.end(), one);

    std::size_t count = {0};
    std::for_each(coefficients.begin(), coefficients.end(), [&count, nd = op.num_dofs](auto &e)
                  { e = count < num_nodes? T(1) : T(S(count-num_nodes));
                    count = (count + 1) % (nd + num_nodes); });

    MPI_Barrier(comm);

    std::array<T, op.num_dofs> Ae;
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

    double max_time = 0;
    MPI_Allreduce(&local_time, &max_time, 1, MPI_DOUBLE, MPI_MAX, comm);

    // Sanity check: Are we computing the correct values?
    for (int batch = 0; batch < num_batches; batch++)
    {
      T acc = 0;
      for (int i = 0; i < local_size; i++)
        acc += A[i + batch * local_size];
      check_solution<T, S>(acc, reference);
    }

    if (mpi_rank == 0)
    {
      std::cout << precision << ", " << batch_size << ", ";
      std::cout << num_cells << ", " << P << ", " << max_time << ", ";
      std::cout << bs << ", " << precompute;
    }
  }
  MPI_Finalize();

  return 0;
}