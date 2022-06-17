#include "src.hpp"

#include <basix/finite-element.h>
#include <basix/quadrature.h>
#include <functional>
#include <libxsmm_source.h>
#include <random>
#include <vector>
#include <xtensor/xio.hpp>
#include <xtensor/xtensor.hpp>
#include <xtensor/xview.hpp>

enum Type {
  Auto,
  OpenMP,
  Blis,
  libxsmm,
  SumFactorization,
  SumFactorizationGEMM
};

// --------------------------------------------------------
const std::string &type_to_str(basix::cell::type type) {
  static const std::map<basix::cell::type, std::string> type_to_name = {
      {basix::cell::type::point, "point"},
      {basix::cell::type::interval, "interval"},
      {basix::cell::type::triangle, "triangle"},
      {basix::cell::type::tetrahedron, "tetrahedron"},
      {basix::cell::type::quadrilateral, "quadrilateral"},
      {basix::cell::type::pyramid, "pyramid"},
      {basix::cell::type::prism, "prism"},
      {basix::cell::type::hexahedron, "hexahedron"}};

  auto it = type_to_name.find(type);
  if (it == type_to_name.end())
    throw std::runtime_error("Can't find type");

  return it->second;
}

// --------------------------------------------------------
const std::string &type_to_str(Type type) {
  static const std::map<Type, std::string> type_to_name = {
      {Type::Auto, "Auto Vectorization"},
      {Type::OpenMP, "OpenMP SIMD"},
      {Type::Blis, "Blis"},
      {Type::libxsmm, "libxsmm"},
      {Type::SumFactorization, "Sum Factorization"},
      {Type::SumFactorizationGEMM, "Sum Factorization GEMM"}};

  auto it = type_to_name.find(type);
  if (it == type_to_name.end())
    throw std::runtime_error("Can't find type");

  return it->second;
}

// --------------------------------------------------------
template <int p> constexpr int num_points() {
  if constexpr (p == 1)
    return 4;
  if constexpr (p == 2)
    return 14;
  if constexpr (p == 3)
    return 24;
  if constexpr (p == 4)
    return 45;
  if constexpr (p == 5)
    return 74;
  if constexpr (p == 6)
    return 122;
  if constexpr (p == 7)
    return 177;
  if constexpr (p == 8)
    return 729;
  if constexpr (p == 9)
    return 1000;
  if constexpr (p == 10)
    return 1331;
}

// --------------------------------------------------------
// Compute the number of dofs.
// Note: Regular mesh.
inline auto cell_dofs(int num_cells) {
  std::map<basix::cell::type, std::vector<int>> celldofs;
  celldofs[basix::cell::type::tetrahedron] = {
      294, 1859, 5776, 13125, 24986, 42439, 66564, 98441, 139150, 189771};
  // celldofs[basix::cell::type::hexahedron] = {1331,   9261,   29791,  68921,
  //                                            132651, 226981, 357911, 531441,
  //                                            753571, 1030301};
  return celldofs;
}

// --------------------------------------------------------
// Tabulate order P basis functions on a cell of type cell_type
xt::xtensor<double, 2> tabulate_basis(int p, int q, basix::cell::type cell) {
  auto family = basix::element::family::P;
  auto quad = basix::quadrature::type::Default;
  auto [points, weights] = basix::quadrature::make_quadrature(quad, cell, q);
  auto variant = basix::element::lagrange_variant::equispaced;
  auto element = basix::create_element(family, cell, p, variant);
  xt::xtensor<double, 4> basis = element.tabulate(0, points);
  return xt::view(basis, 0, xt::all(), xt::all(), 0);
}

// --------------------------------------------------------
// Tabulate order P basis functions on a Hexahedron
auto tabulate_basis_derivatives(int p, int q, basix::cell::type cell) {
  auto family = basix::element::family::P;
  auto quad = basix::quadrature::type::Default;
  auto [points, weights] = basix::quadrature::make_quadrature(quad, cell, q);
  auto variant = basix::element::lagrange_variant::equispaced;
  auto element = basix::create_element(family, cell, p, variant);
  xt::xtensor<double, 4> basis = element.tabulate(1, points);
  xt::xtensor<double, 2> dphi0 = xt::view(basis, 1, xt::all(), xt::all(), 0);
  xt::xtensor<double, 2> dphi1 = xt::view(basis, 2, xt::all(), xt::all(), 0);
  xt::xtensor<double, 2> dphi2 = xt::view(basis, 3, xt::all(), xt::all(), 0);
  return std::tuple{dphi0, dphi1, dphi2};
}

// --------------------------------------------------------
auto init_vectors(int p, int ncells, basix::cell::type cell) {
  xt::xtensor<double, 2> phi = tabulate_basis(p, 2 * p, cell);
  std::vector<double> u(ncells * phi.shape(1), 0);
  std::vector<double> w(ncells * phi.shape(0), 0);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0, 1);
  std::generate(u.begin(), u.end(), [&]() { return dis(gen); });
  return std::tuple{phi, u, w};
}

// --------------------------------------------------------
auto init_vectors_stiffness(int p, int ncells, basix::cell::type cell) {
  auto [dphi0, dphi1, dphi2] = tabulate_basis_derivatives(p, 2 * p, cell);
  std::vector<double> u(ncells * dphi0.shape(1), 0);
  std::vector<double> w(3 * ncells * dphi0.shape(0), 0);
  std::random_device rd;
  std::mt19937 gen(rd());
  std::uniform_real_distribution<double> dis(0, 1);
  std::generate(u.begin(), u.end(), [&]() { return dis(gen); });
  return std::tuple{dphi0, dphi1, dphi2, u, w};
}

// --------------------------------------------------------
template <typename T, int p> auto generate_mass(basix::cell::type cell) {
  std::map<Type, std::function<void(T *, T *, T *)>> kernels;
  if (cell == basix::cell::type::hexahedron) {
    constexpr int nq = (p + 1) * (p + 1) * (p + 1);
    constexpr int nd = (p + 1) * (p + 1) * (p + 1);
    kernels[Type::Auto] = [=](T *phi, T *u, T *w) {
      compute_coefficients_auto<T, nq, nd>(phi, u, w);
    };
    kernels[Type::OpenMP] = [=](T *phi, T *u, T *w) {
      compute_coefficients_openmp<T, nq, nd>(phi, u, w);
    };
    kernels[Type::Blis] = [=](T *phi, T *u, T *w) {
      compute_coefficients_blis<T, nq, nd>(phi, u, w);
    };
    kernels[Type::SumFactorization] = [=](T *phi, T *u, T *w) {
      compute_coefficients_sf<T, p + 1, p + 1>(phi, u, w);
    };
    kernels[Type::SumFactorizationGEMM] = [=](T *phi, T *u, T *w) {
      compute_coefficients_sfgemm<T, p + 1, p + 1>(phi, u, w);
    };

    typedef libxsmm_mmfunction<T> kernel_type;
    kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, nq, 1, nd, 1.0, 0.0);
    kernels[Type::libxsmm] = [=](T *phi, T *u, T *w) { kernel(phi, u, w); };

  } else if (cell == basix::cell::type::tetrahedron) {
    constexpr int nd = (p + 1) * (p + 2) * (p + 3) / 6;
    constexpr int nq = num_points<p>();
    kernels[Type::Auto] = [=](T *phi, T *u, T *w) {
      compute_coefficients_auto<T, nq, nd>(phi, u, w);
    };
    // kernels[Type::OpenMP] = [=](T *phi, T *u, T *w) {
    //   compute_coefficients_openmp<T, nq, nd>(phi, u, w);
    // };
    kernels[Type::Blis] = [=](T *phi, T *u, T *w) {
      compute_coefficients_blis<T, nq, nd>(phi, u, w);
    };

    // typedef libxsmm_mmfunction<T> kernel_type;
    // kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, nq, 1, nd, 1.0, 0.0);
    // kernels[Type::libxsmm] = [=](T *phi, T *u, T *w) { kernel(phi, u, w); };
  }
  return kernels;
}

template <typename T, int p> auto generate_stiffness(basix::cell::type cell) {
  std::map<Type, std::function<void(T *, T *, T *, T *, T *)>> kernels;
  if (cell == basix::cell::type::hexahedron) {
    constexpr int nq = (p + 1) * (p + 1) * (p + 1);
    constexpr int nd = (p + 1) * (p + 1) * (p + 1);
    kernels[Type::Auto] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_auto<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };
    kernels[Type::OpenMP] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_openmp<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };

    kernels[Type::Blis] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_blis<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };

    typedef libxsmm_mmfunction<T> kernel_type;
    kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, nq, 1, nd, 1.0, 0.0);
    kernels[Type::libxsmm] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      kernel(dphi0, u, w);
      kernel(dphi1, u, w + nq);
      kernel(dphi2, u, w + 2 * nq);
    };
  } else if (cell == basix::cell::type::tetrahedron) {
    constexpr int nd = (p + 1) * (p + 2) * (p + 3) / 6;
    constexpr int nq = num_points<p>();
    kernels[Type::Auto] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_auto<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };
    kernels[Type::OpenMP] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_openmp<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };

    kernels[Type::Blis] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      compute_coefficients_blis<T, nq, nd>(dphi0, dphi1, dphi2, u, w);
    };

    typedef libxsmm_mmfunction<T> kernel_type;
    kernel_type kernel(LIBXSMM_GEMM_FLAG_NONE, nq, 1, nd, 1.0, 0.0);
    kernels[Type::libxsmm] = [=](T *dphi0, T *dphi1, T *dphi2, T *u, T *w) {
      kernel(dphi0, u, w);
      kernel(dphi1, u, w + nq);
      kernel(dphi2, u, w + 2 * nq);
    };
  }
  return kernels;
}

// --------------------------------------------------------
template <typename T> auto mass_coefficient_kernels(basix::cell::type cell) {
  std::vector<std::map<Type, std::function<void(T *, T *, T *)>>> kernels;
  kernels.push_back(generate_mass<double, 1>(cell));
  kernels.push_back(generate_mass<double, 2>(cell));
  kernels.push_back(generate_mass<double, 3>(cell));
  kernels.push_back(generate_mass<double, 4>(cell));
  kernels.push_back(generate_mass<double, 5>(cell));
  kernels.push_back(generate_mass<double, 6>(cell));
  kernels.push_back(generate_mass<double, 7>(cell));
  kernels.push_back(generate_mass<double, 8>(cell));
  kernels.push_back(generate_mass<double, 9>(cell));
  kernels.push_back(generate_mass<double, 10>(cell));
  return kernels;
}

// --------------------------------------------------------
template <typename T>
auto stiffness_coefficient_kernels(basix::cell::type cell) {
  std::vector<std::map<Type, std::function<void(T *, T *, T *, T *, T *)>>>
      kernels;
  kernels.push_back(generate_stiffness<double, 1>(cell));
  kernels.push_back(generate_stiffness<double, 2>(cell));
  kernels.push_back(generate_stiffness<double, 3>(cell));
  kernels.push_back(generate_stiffness<double, 4>(cell));
  kernels.push_back(generate_stiffness<double, 5>(cell));
  kernels.push_back(generate_stiffness<double, 6>(cell));
  kernels.push_back(generate_stiffness<double, 7>(cell));
  kernels.push_back(generate_stiffness<double, 8>(cell));
  kernels.push_back(generate_stiffness<double, 9>(cell));
  kernels.push_back(generate_stiffness<double, 10>(cell));
  return kernels;
}

void output(int p, Type method, double t, double num_dofs,
            basix::cell::type cell_type) {
  std::cout << p << ", ";
  std::cout << type_to_str(method) << ", ";
  std::cout << t << ", ";
  std::cout << num_dofs << ", ";
  std::cout << type_to_str(cell_type) << std::endl;
}