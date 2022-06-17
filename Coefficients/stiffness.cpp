#include "utils.hpp"
#include <basix/cell.h>

using duration = std::chrono::duration<double>;
int main() {
  constexpr int ncells = 10000;
  constexpr int max_degree = 10;
  constexpr int nrepeats = 5;
  auto celldofs = cell_dofs(ncells);

  // Two cell types are available (Hexes and Tets)
  for (auto [cell_type, dofs] : celldofs) {
    auto kernels = stiffness_coefficient_kernels<double>(cell_type);
    // Select Polynomial Degree
    for (int p = 1; p <= max_degree; p++) {
      auto [dphi0, dphi1, dphi2, u, w] =
          init_vectors_stiffness(p, ncells, cell_type);
      auto kernel_p = kernels[p - 1];
      // For each cell_type/degree multiple implementations are available
      for (auto [method, compute] : kernel_p) {
        int nquads = dphi0.shape(0);
        int ndofs = dphi0.shape(1);
        // Repeat computation - statistical purposes
        for (int r = 0; r < nrepeats; r++) {
          auto start = std::chrono::high_resolution_clock::now();
          // Actual Computation
          for (int c = 0; c < ncells; c++) {
            double *_u = u.data() + c * ndofs;
            double *_w = w.data() + c * nquads;
            compute(dphi0.data(), dphi1.data(), dphi2.data(), _u, _w);
          }
          auto end = std::chrono::high_resolution_clock::now();
          auto t = std::chrono::duration_cast<duration>(end - start);
          // Print output to terminal
          output(p, method, t.count() / 10, dofs[p - 1], cell_type);
        }
      }
    }
  }
}