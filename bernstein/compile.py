
import typing
from importlib import reload
import scipy.special

_headers = """
#include <cmath>
#include <cstring>
#include <cstdint>
#include <iostream>
#define restrict __restrict__
#include<time.h>

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

    code1 = f"""
    void kernel({scalar_type}* restrict f3,
                const {scalar_type}* restrict c0,
                const {scalar_type}* restrict consts,
                const {scalar_type}* restrict coordinate_dofs,
                const int* restrict entity_local_index,
                const uint8_t* restrict quadrature_permutation)
    {{

  // Input: c0 ({(n)*(n+1)*(n+2)//6}) - dofs
  // Output: f3 ({(n)*(n+1)*(n+2)//6}) - dofs

  double J[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
      J[i][j] = coordinate_dofs[(i+1)*3 + j] - coordinate_dofs[j];

  // Cofactors
  double C[3][3];
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
    {{
      int i1 = (i + 1) % 3;
      int i2 = (i + 2) % 3;
      int j1 = (j + 1) % 3;
      int j2 = (j + 2) % 3;
      C[i][j] = J[i1][j1] * J[i2][j2] - J[i1][j2] * J[i2][j1];
    }}

  // FIXME: scale by detJ
  double JtJ[3][3] = {{}};
  for (int i = 0; i < 3; ++i)
    for (int j = 0; j < 3; ++j)
        for (int k = 0; k < 3; ++k)
          JtJ[i][k] += C[i][j] * C[k][j];


  double rule2w[{q}] = {{{', '.join([str(p) for p in rule2[1]])}}};
  double rule1w[{q}] = {{{', '.join([str(p) for p in rule1[1]])}}};
  double rule0w[{q}] = {{{', '.join([str(p) for p in rule0[1]])}}};
  double rule2p[{q}] = {{{', '.join([str(p) for p in rule2[0]])}}};
  double rule1p[{q}] = {{{', '.join([str(p) for p in rule1[0]])}}};
  double rule0p[{q}] = {{{', '.join([str(p) for p in rule0[0]])}}};

  double c1[3][{n}][{n}][{q}] = {{}};
  int dx[3] = {{1, 1, 1}};

  for (int i2 = 0; i2 < {q}; ++i2)
  {{
    double s = 1.0 - rule0p[i2];
    double r = rule0p[i2] / s;
    int c = 0;
    dx[0] = {(n+2)*(n+1)//2};
    for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
    {{
      dx[1] = {n + 1} - alpha1;
      for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
      {{
        double w = 1.0;
        for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
          w *= s;
        for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
        {{
          for (int dim = 0; dim < 3; ++dim)
            c1[dim][alpha1][alpha2][i2] += {n} * w * c0[c + dx[dim]];
          for (int dim = 0; dim < 3; ++dim)
            c1[dim][alpha1][alpha2][i2] -= {n} * w * c0[c];
          ++c;
          w *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
        }}
        ++c;
        --dx[1];
        --dx[0];
      }}
      ++c;
      --dx[0];
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

  double f1[3][{n}][{q}][{q}] = {{}};
  for (int i0 = 0; i0 < {q}; ++i0)
  {{
    double s = 1.0 - rule2p[i0];
    double r = rule2p[i0] / s;
    double w[{n}] = {{1.0}};
    for (int j = 0; j < {n - 1}; ++j)
      w[0] *= s;
    for (int j = 0; j < {n - 1}; ++j)
      w[j + 1] = w[j] * r * ({n - 1} - j) / (1.0 + j);

    for (int i1 = 0; i1 < {q}; ++i1)
      for (int i2 = 0; i2 < {q}; ++i2)
      {{
        double c3i[3] = {{}};
        for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
          for (int dim = 0; dim < 3; ++dim)
            c3i[dim] += w[alpha1] * c2[dim][alpha1][i1][i2];

        // Transform by Jacobian
        double f3i[3] = {{}};
        for (int i = 0; i < 3; ++i)
          for (int j = 0; j < 3; ++j)
            f3i[i] += JtJ[i][j] * c3i[j];

        for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
          for (int dim = 0; dim < 3; ++dim)
            f1[dim][alpha1][i1][i2] += w[alpha1] * rule2w[i0] * f3i[dim];
      }}
  }}


  memset(c1, 0, {3*n*n*q}*sizeof(double));

  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i1 = 0; i1 < {q}; ++i1)
    {{
      double s = 1 - rule1p[i1];
      double r = rule1p[i1] / s;
      double w = rule1w[i1];
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        double ww = w;
        for (int j = 0; j < {n - 1} - alpha1; ++j)
          ww *= s;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          for (int i2 = 0; i2 < {q}; ++i2)
            c1[dim][alpha1][alpha2][i2] += ww * f1[dim][alpha1][i1][i2];
          ww *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}
  }}

  // double f3[{(n + 3) * (n + 2) * (n + 1) // 6}]

    for (int i2 = 0; i2 < {q}; ++i2)
    {{
      double s = 1.0 - rule0p[i2];
      double r = rule0p[i2] / s;
      double w = rule0w[i2];
      int c = 0;
      dx[0] = {(n+2)*(n+1)//2};
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        dx[1] = {n + 1} - alpha1;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          double ww = {n} * w;
          for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
            ww *= s;
          for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
          {{
             for (int dim = 0; dim < 3; ++dim)
               f3[c] -= ww * c1[dim][alpha1][alpha2][i2];
             for (int dim = 0; dim < 3; ++dim)
               f3[c + dx[dim]] += ww * c1[dim][alpha1][alpha2][i2];
             ww *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
             ++c;
          }}
          ++c;
          --dx[1];
          --dx[0];
        }}
        ++c;
        --dx[0];
      }}
    }}

}}

    """

    code2 = f"""
    void kernel({scalar_type}* restrict f3,
                   const {scalar_type}* restrict c0,
                   const {scalar_type}* restrict consts,
                   const {scalar_type}* restrict coordinate_dofs,
                   const int* restrict entity_local_index,
                   const uint8_t* restrict quadrature_permutation)
    {{

  // Input: c0 ({(n+1)*(n+2)*(n+3)//6}) - dofs
  // Output: f3 ({(n+1)*(n+2)*(n+3)//6}) - dofs

  double rule2w[{q}] = {{{', '.join([str(p) for p in rule2[1]])}}};
  double rule1w[{q}] = {{{', '.join([str(p) for p in rule1[1]])}}};
  double rule0w[{q}] = {{{', '.join([str(p) for p in rule0[1]])}}};
  double rule2p[{q}] = {{{', '.join([str(p) for p in rule2[0]])}}};
  double rule1p[{q}] = {{{', '.join([str(p) for p in rule1[0]])}}};
  double rule0p[{q}] = {{{', '.join([str(p) for p in rule0[0]])}}};

  // Copy of dofs in B(n-1) for dx dy and dz
  double dofs[{(n+2)*n*(n+1)//6}][4] = {{}};
  int dx = {(n+2)*(n+1)//2};
  int c = 0;
  int k = 0;
  for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
  {{
    int dy = {n + 1} - alpha1;
    for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
    {{
      for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
      {{
        dofs[k][0] = c0[c + dx];
        dofs[k][1] = c0[c + dy];
        dofs[k][2] = c0[c + 1];
        dofs[k][3] = c0[c];
        ++c;
        ++k;
      }}
      ++c;
      --dy;
      --dx;
    }}
    ++c;
    --dx;
  }}

  for (int dim = 0; dim < 4; ++dim)
  {{
  double c1[{n}][{n}][{q}] = {{}};
  for (int i2 = 0; i2 < {q}; ++i2)
  {{
    double s = 1.0 - rule0p[i2];
    double r = rule0p[i2] / s;
    int c = 0;
    for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
    {{
      for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
      {{
        double w = {n};
        for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
          w *= s;
        for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
        {{
          c1[alpha1][alpha2][i2] += w * dofs[c][dim];
          ++c;
          w *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
        }}
      }}
    }}
  }}

  double c2[{n}][{q}][{q}] = {{}};
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
              c2[alpha1][i1][i2] += w * c1[alpha1][alpha2][i2];
          w *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}


    double f1[{n}][{q}][{q}] = {{}};
    for (int i0 = 0; i0 < {q}; ++i0)
    {{
      double c3[{q}][{q}] = {{}};
      double s = 1.0 - rule2p[i0];
      double r = rule2p[i0] / s;
      double w = 1.0;
      for (int j = 0; j < {n - 1}; ++j)
        w *= s;
      double ww = w * rule2w[i0];
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        for (int i1 = 0; i1 < {q}; ++i1)
          for (int i2 = 0; i2 < {q}; ++i2)
            c3[i1][i2] += w * c2[alpha1][i1][i2];
        w *= r * ({n - 1} - alpha1) / (1.0 + alpha1);
      }}
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        for (int i1 = 0; i1 < {q}; ++i1)
          for (int i2 = 0; i2 < {q}; ++i2)
            f1[alpha1][i1][i2] += ww * c3[i1][i2];
        ww *= r * ({n - 1} - alpha1) / (1.0 + alpha1);
      }}
    }}


    for (int i = 0; i < {n}; ++i)
      for (int j = 0; j < {n}; ++j)
        for (int k = 0; k < {q}; ++k)
            c1[i][j][k] = 0;

  for (int i1 = 0; i1 < {q}; ++i1)
  {{
    double s = 1 - rule1p[i1];
    double r = rule1p[i1] / s;
    double w = rule1w[i1];
    for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
    {{
      double ww = w;
      for (int j = 0; j < {n - 1} - alpha1; ++j)
        ww *= s;
      for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
      {{
        for (int i2 = 0; i2 < {q}; ++i2)
          c1[alpha1][alpha2][i2] += ww * f1[alpha1][i1][i2];
        ww *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
      }}
    }}
  }}


  // double f3[{(n + 3) * (n + 2) * (n + 1) // 6}]

    for (int i2 = 0; i2 < {q}; ++i2)
    {{
      double s = 1.0 - rule0p[i2];
      double r = rule0p[i2] / s;
      double w = rule0w[i2];
      int c = 0;
      int dx = {(n+2)*(n+1)//2};
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        int dy = {n + 1} - alpha1;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          double ww = {n} * w;
          for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
            ww *= s;
          for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
          {{
             int d = dx;
             if (dim == 1)
                d = dy;
             else if (dim == 2)
                d = 1;
             f3[c] -= ww * c1[alpha1][alpha2][i2];
             f3[c + d] += ww * c1[alpha1][alpha2][i2];
             ww *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
             ++c;
          }}
          ++c;
          --dy;
          --dx;
        }}
        ++c;
        --dx;
      }}
    }}
  }}


}}
    """

    code3 = f"""
        void kernel({scalar_type}* restrict f3,
                   const {scalar_type}* restrict c0,
                   const {scalar_type}* restrict consts,
                   const {scalar_type}* restrict coordinate_dofs,
                   const int* restrict entity_local_index,
                   const uint8_t* restrict quadrature_permutation)
{{
  // Input: c0 ({(n)*(n+1)*(n+2)//6}) - dofs
  // Output: f3 ({(n)*(n+1)*(n+2)//6}) - dofs

  double rule2w[{q}] = {{{', '.join([str(p) for p in rule2[1]])}}};
  double rule1w[{q}] = {{{', '.join([str(p) for p in rule1[1]])}}};
  double rule0w[{q}] = {{{', '.join([str(p) for p in rule0[1]])}}};
  double rule2p[{q}] = {{{', '.join([str(p) for p in rule2[0]])}}};
  double rule1p[{q}] = {{{', '.join([str(p) for p in rule1[0]])}}};
  double rule0p[{q}] = {{{', '.join([str(p) for p in rule0[0]])}}};

  double c1[3][{n}][{n}][{q}];

  for (int i2 = 0; i2 < {q}; ++i2)
  {{
    double s = 1.0 - rule0p[i2];
    double r = rule0p[i2] / s;
    int c = 0;
    int dx = {(n+2)*(n+1)//2};
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
          c1vx += w * (c0[c + dx] - c0[c]);
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
        --dx;
      }}
      ++c;
      --dx;
    }}
  }}

  std::cout << "t0 = " <<  clock() - time << std::endl;
  time = clock();

  double c2[3][{n}][{q*q}] = {{}};
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
            c2[dim][alpha1][i1 * {q} + i2] += w * c1[dim][alpha1][alpha2][i2];
          w *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}
  }}

  std::cout << "t1 = " <<  clock() - time << std::endl;
  time = clock();

  double f1[3][{n}][{q*q}] = {{}};
  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i0 = 0; i0 < {q}; ++i0)
    {{
      double s = 1.0 - rule2p[i0];
      double r = rule2p[i0] / s;
      double w[{n}] = {{1.0}};
      for (int j = 0; j < {n - 1}; ++j)
        w[0] *= s;
      for (int j = 0; j < {n - 1}; ++j)
        w[j + 1] = w[j] * r * ({n - 1} - j) / (1.0 + j);

      for (int i = 0; i < {q*q}; ++i)
      {{
        double c3i = 0;
        for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
          c3i += w[alpha1] * c2[dim][alpha1][i];
        for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
          f1[dim][alpha1][i] += w[alpha1] * rule2w[i0] * c3i;
      }}
    }}
  }}

    std::cout << "t2 = " <<  clock() - time << std::endl;
  time = clock();

  memset(c1, 0, {3*n*n*q}*sizeof(double));

  for (int dim = 0; dim < 3; ++dim)
  {{
    for (int i1 = 0; i1 < {q}; ++i1)
    {{
      double s = 1 - rule1p[i1];
      double r = rule1p[i1] / s;
      double w = rule1w[i1];
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        double ww = w;
        for (int j = 0; j < {n - 1} - alpha1; ++j)
          ww *= s;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          for (int i2 = 0; i2 < {q}; ++i2)
            c1[dim][alpha1][alpha2][i2] += ww * f1[dim][alpha1][i1*{q} + i2];
          ww *= r * ({n - 1} - alpha1 - alpha2) / (1.0 + alpha2);
        }}
      }}
    }}
  }}

  std::cout << "t4 = " <<  clock() - time << std::endl;
  time = clock();

  // double f3[{(n + 3) * (n + 2) * (n + 1) // 6}]

    for (int i2 = 0; i2 < {q}; ++i2)
    {{
      double s = 1.0 - rule0p[i2];
      double r = rule0p[i2] / s;
      double w = rule0w[i2];
      int c = 0;
      int dx = {(n+2)*(n+1)//2};
      for (int alpha1 = 0; alpha1 < {n}; ++alpha1)
      {{
        int dy = {n + 1} - alpha1;
        for (int alpha2 = 0; alpha2 < {n} - alpha1; ++alpha2)
        {{
          double ww = {n} * w;
          for (int j = 0; j < {n - 1} - alpha1 - alpha2; ++j)
            ww *= s;
          for (int alpha3 = 0; alpha3 < {n} - alpha1 - alpha2; ++alpha3)
          {{
             f3[c] -= ww * (c1[0][alpha1][alpha2][i2] + c1[1][alpha1][alpha2][i2] + c1[2][alpha1][alpha2][i2]);
             f3[c + dx] += ww * c1[0][alpha1][alpha2][i2];
             f3[c + dy] += ww * c1[1][alpha1][alpha2][i2];
             f3[c + 1] += ww * c1[2][alpha1][alpha2][i2];
             ww *= r * ({n - 1} - alpha1 - alpha2 - alpha3)/(1.0 + alpha3);
             ++c;
          }}
          ++c;
          --dy;
          --dx;
        }}
        ++c;
        --dx;
      }}
    }}

  std::cout << "t5 = " <<  clock() - time << std::endl;


}}

    """

    with open("bernstein/problem.hpp", "w") as file:
        file.write(headers)
        file.write(code1)
