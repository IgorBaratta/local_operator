
import typing
from importlib import reload
import scipy.special

_headers = """
#include <cmath>
#include <cstdint>
#define restrict __restrict__

constexpr int dim = {dim};
constexpr int global_size = {global_size};
constexpr int kernel_rank = {rank};

using scalar_type={scalar_type};
using namespace std;
"""


def generate_code(n, scalar_type, global_size):

    ndofs = (n + 1)*(n + 2)*(n + 3)//6

    headers = _headers.format(dim=ndofs, global_size=global_size,
                              scalar_type=scalar_type, rank=1)

    q = n + 1
    rule0 = scipy.special.roots_jacobi(q, 0, 0)
    rule1 = scipy.special.roots_jacobi(q, 1, 0)
    rule2 = scipy.special.roots_jacobi(q, 2, 0)
    rule0 = ((rule0[0] + 1) / 2, rule0[1] / 2)
    rule1 = ((rule1[0] + 1) / 2, rule1[1] / 4)
    rule2 = ((rule2[0] + 1) / 2, rule2[1] / 8)

    code= f"""
    void kernel({scalar_type}* restrict f3,
                   const {scalar_type}* restrict c0,
                   const {scalar_type}* restrict c,
                   const {scalar_type}* restrict coordinate_dofs,
                   const int* restrict entity_local_index,
                   const uint8_t* restrict quadrature_permutation)
    {{

  // Input: c0 ({(n)*(n+1)*(n+2)//6}) - dofs
  // Output: c3 (3, {q*q*q}) - gradient values at quadrature points
  double rule2p[{q}] = {{{', '.join([str(p) for p in rule2[0]])}}};
  double rule1p[{q}] = {{{', '.join([str(p) for p in rule1[0]])}}};
  double rule0p[{q}] = {{{', '.join([str(p) for p in rule0[0]])}}};
  double rule2w[{q}] = {{{', '.join([str(p) for p in rule2[1]])}}};
  double rule1w[{q}] = {{{', '.join([str(p) for p in rule1[1]])}}};
  double rule0w[{q}] = {{{', '.join([str(p) for p in rule0[1]])}}};

  double c1[3][{n}][{n}][{q}];

  for (int i2 = 0; i2 < {q}; ++i2)
  {{
    double s = 1.0 - rule0p[i2];
    double r = rule0p[i2] / s;
    int c = 0;
    int dz = {(n+2)*(n+1)//2};
    for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
    {{
      int dy = {n + 1} - alpha1;
      for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
      {{
        double w = 1.0;
        for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
          w *= s;
        double c1vx = 0.0;
        double c1vy = 0.0;
        double c1vz = 0.0;
        for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
        {{
          c1vx += w * (c0[c + dz] - c0[c]);
          c1vy += w * (c0[c + dy] - c0[c]);
          c1vz += w * (c0[c + 1] - c0[c]);
          ++c;
          w *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
        }}
        c1[0][alpha1][alpha2][i2] = {n} * c1vx;
        c1[1][alpha1][alpha2][i2] = {n} * c1vy;
        c1[2][alpha1][alpha2][i2] = {n} * c1vz;
        ++c;
        --dy;
        --dz;
      }}
      ++c;
      --dz;
    }}
  }}

  double c2[3][{n}][{q}][{q}] = {{}};

  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i1 = 0; i1 < {q}; ++i1)
    {{
      double s = 1.0 - rule1p[i1];
      double r = rule1p[i1]/s;
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        double w = 1.0;
        for (int j = 0; j < {n - 1} - alpha1; ++j)
          w *= s;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          for (int i2 = 0; i2 < {q}; ++i2)
            c2[dim][alpha1][i1][i2] += w * c1[dim][alpha1][alpha2][i2];
          w *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}
  }}

  double c3[3][{q}][{q}][{q}] = {{0}};

  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i0 = 0; i0 < {q}; ++i0)
    {{
      double s = 1.0 - rule2p[i0];
      double r = rule2p[i0] / s;
      double w = 1.0;
      for (int j = 0; j < {n - 1}; ++j)
        w *= s;
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        for (int i1 = 0; i1 < {q}; ++i1)
          for (int i2 = 0; i2 < {q}; ++i2)
            c3[dim][i0][i1][i2] += w * c2[dim][alpha1][i1][i2];
         w *= r * ({n - 1} - alpha1) / (1.0 + alpha1);
      }}
    }}
  }}


  double f1[3][{n + 1}][{q}][{q}] = {{}};
  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i0 = 0; i0 < {q}; ++i0)
    {{
      double s = 1.0 - rule2p[i0];
      double r = rule2p[i0] / s;
      double ww = rule2w[i0];
      for (int j = 0; j < {n}; ++j)
        ww *= s;
      for (int alpha1 = 0; alpha1 < {n + 1}; ++alpha1)
      {{
        for (int i1 = 0; i1 < {q}; ++i1)
          for (int i2 = 0; i2 < {q}; ++i2)
            f1[dim][alpha1][i1][i2] += ww * c3[dim][i0][i1][i2];
        ww *= r * ({n} - alpha1) / (1.0 + alpha1);
      }}
    }}
  }}

  double f2[3][{n+1}][{n+1}][{q}] = {{}};
  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i1 = 0; i1 < {q}; ++i1)
    {{
      double s = 1 - rule1p[i1];
      double r = rule1p[i1] / s;
      double w = rule1w[i1];
      for (int alpha1 = 0; alpha1 < {n + 1}; ++alpha1)
      {{
        double ww = w;
        for (int j = 0; j < {n} - alpha1; ++j)
          ww *= s;
        for (int alpha2 = 0; alpha2 < {n+1} - alpha1; ++alpha2)
        {{
          for (int i2 = 0; i2 < {q}; ++i2)
            f2[dim][alpha1][alpha2][i2] += ww * f1[dim][alpha1][i1][i2];
          ww *= r * ({n} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}
  }}

    // double f3[{(n + 3) * (n + 2) * (n + 1) // 6}]
    for (int dim = 0; dim < 3; ++dim)
    {{
      for (int i2 = 0; i2 < {q}; ++i2)
      {{
        double s = 1.0 - rule0p[i2];
        double r = rule0p[i2] / s;
        double w = rule0w[i2];
        int c = 0;
        for (int alpha1 = 0; alpha1 < {n + 1}; ++alpha1)
          for (int alpha2 = 0; alpha2 < {n + 1} - alpha1; ++alpha2)
          {{
            double ww = w;
            for (int j = 0; j < {n} - alpha1 - alpha2; ++j)
              ww *= s;
            for (int alpha3 = 0; alpha3 < {n + 1} - alpha1 - alpha2; ++alpha3)
            {{
               f3[c++] += ww * f2[dim][alpha1][alpha2][i2];
               ww *= r * ({n} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
            }}
          }}
        }}
      }}

    }}

    """


    with open("bernstein/problem.hpp", "w") as file:
        file.write(headers)
        file.write(code)
