#include "mass.h"
#include "kernel.hpp"
#include <dolfinx.h>
#include <dolfinx/fem/assemble_matrix_impl.h>

int main(int argc, char *argv[])
{
    common::subsystem::init_logging(argc, argv);
    common::subsystem::init_mpi(argc, argv);

    MPI_Comm mpi_comm{MPI_COMM_WORLD};

    std::shared_ptr<mesh::Mesh> mesh = std::make_shared<mesh::Mesh>(
        generation::BoxMesh::create(mpi_comm, {{{0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}}}, {20, 20, 20},
                                    mesh::CellType::tetrahedron, mesh::GhostMode::none));

    mesh->topology().create_entity_permutations();
    const std::vector<std::uint32_t> &cell_info = mesh->topology().get_cell_permutation_info();

    auto constants = std::map<std::string, std::shared_ptr<const fem::Constant<double>>>();
    const std::shared_ptr<fem::FunctionSpace> &V = fem::create_functionspace(functionspace_form_mass_a, "u", mesh);
    std::shared_ptr<fem::Form<double>> a = std::make_shared<fem::Form<double>>(fem::create_form<double>(
        *form_mass_a, {V, V}, {}, constants, {}));

    int tdim = mesh->topology().dim();
    const auto &topology = mesh->topology();
    int ncells = topology.index_map(tdim)->size_global();
    int ndofs_cell = V->element()->space_dimension();

    std::vector<std::int32_t> active_cells(ncells);
    std::iota(active_cells.begin(), active_cells.end(), 0);

    std::function<int(std::int32_t, const std::int32_t *, std::int32_t,
                      const std::int32_t *, const double *)>
        insert_block = [](std::int32_t nr, const std::int32_t *rows, const std::int32_t nc,
                          const std::int32_t *cols, const double *data) {
            // do nothing
            return 0;
        };

    // Call std function from kernel.hpp
    //=====================================================================================
    std::function<void(double *, const double *, const double *, const double *, const int *,
                       const std::uint8_t *, const std::uint32_t)>
        kernel = [](double *A, const double *w, const double *c,
                    const double *coordinate_dofs, const int *entity_local_index,
                    const uint8_t *quadrature_permutation, const uint32_t cell_permutation) {
            tabulate_tensor(A, w, c, coordinate_dofs, entity_local_index, quadrature_permutation,
                            cell_permutation);
        };

    const auto &geometry = mesh->geometry();
    auto dofmap0 = V->dofmap()->list();
    auto dofmap1 = V->dofmap()->list();
    int bs0 = 1;
    int bs1 = 1;
    std::vector<bool> bc0(0);
    std::vector<bool> bc1(0);
    array2d<double> coeffs(0, 0);

    double t_func = MPI_Wtime();
    dolfinx::fem::impl::assemble_cells(insert_block, geometry, active_cells, dofmap0, bs0, dofmap1, bs1, bc0,
                                       bc1, kernel, coeffs, {}, cell_info);
    t_func = MPI_Wtime() - t_func;
    std::cout << "\n std::function \t " << t_func;


    // Call std function from form
    //=====================================================================================
    auto ids = a->integral_ids(dolfinx::fem::IntegralType::cell);
    auto a_kernel = a->kernel(dolfinx::fem::IntegralType::cell, ids[0]);
    double t_form = MPI_Wtime();
    dolfinx::fem::impl::assemble_cells(insert_block, geometry, active_cells, dofmap0, bs0, dofmap1, bs1,
                                       bc0, bc1, a_kernel, coeffs, {}, cell_info);

    t_form = MPI_Wtime() - t_form;
    std::cout << "\n assemble form \t " << t_form;

    return 0;
}