
#include <cmath>
#include <cstdint>
#define restrict __restrict__

template <typename T>
void compute_geometry_tensor(const T *restrict coords, T K[3][3])
{
  // Quadrature independent code
  T J[3][3] = {0};
  J[0][0] = coords[1] - coords[0];  // 1 additions
  J[0][1] = coords[2] - coords[0];  // 1 additions
  J[0][2] = coords[3] - coords[0];  // 1 additions
  J[1][0] = coords[5] - coords[4];  // 1 additions
  J[1][1] = coords[6] - coords[4];  // 1 additions
  J[1][2] = coords[7] - coords[4];  // 1 additions
  J[2][0] = coords[9] - coords[8];  // 1 additions
  J[2][1] = coords[10] - coords[8]; // 1 additions
  J[2][2] = coords[11] - coords[8]; // 1 additions

  // Determinant of jacobian
  // FLOPS: 3 * 4
  T detJ = J[0][0] * (J[1][1] * J[2][2] - J[1][2] * J[2][1]) - // 2 additions, 3 multiplication
           J[0][1] * (J[1][0] * J[2][2] - J[1][2] * J[2][0]) + // 2 additions, 3 multiplication
           J[0][2] * (J[1][0] * J[2][1] - J[1][1] * J[2][0]);  // 1 additions, 3 multiplication

  // FLOPS : 9 * 4
  K[0][0] = (J[1][1] * J[2][2] - J[1][2] * J[2][1]) / detJ; // 1 additions, 3 multiplication
  K[0][1] = (J[0][2] * J[2][1] - J[0][1] * J[2][2]) / detJ; // 1 additions, 3 multiplication
  K[0][2] = (J[0][1] * J[1][2] - J[0][2] * J[1][1]) / detJ; // 1 additions, 3 multiplication
  K[1][0] = (J[1][2] * J[2][0] - J[1][0] * J[2][2]) / detJ; // 1 additions, 3 multiplication
  K[1][1] = (J[0][0] * J[2][2] - J[0][2] * J[2][0]) / detJ; // 1 additions, 3 multiplication
  K[1][2] = (J[0][2] * J[1][0] - J[0][0] * J[1][2]) / detJ; // 1 additions, 3 multiplication
  K[2][0] = (J[1][0] * J[2][1] - J[1][1] * J[2][0]) / detJ; // 1 additions, 3 multiplication
  K[2][1] = (J[0][1] * J[2][0] - J[0][0] * J[2][1]) / detJ; // 1 additions, 3 multiplication
  K[2][2] = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) / detJ; // 1 additions, 3 multiplication
}

template <typename T, typename S, int P>
struct Operator
{
  constexpr static int num_dofs = (P + 1) * (P + 2) * (P + 3) / 2;
  constexpr static int Nc = (P + 1) * (P + 2) * (P + 3) / 6;
  constexpr static int Nq = (P + 1) * (P + 2) * (P + 3) / 2;

  template <int Np>
  inline void kernel(T *restrict A, const T *restrict w, const T K[3][3], const T &detJ, int offset)
  {
    T w1[3][3][Np] = {0};
    T w0[Np] = {0};

    // Restrict finite element coefficients to values values at Np quadrature points
    // Note: Auto-vectorize across degrees of freedom.
    for (int ip = 0; ip < Np; ip++)
    {
      for (int ic = 0; ic < Nc; ++ic)
      {
        w1[0][0][ip] += w[4 + ic] * D0[offset + ip][ic];          // 1 additions, 1 multiplications
        w1[0][1][ip] += w[4 + ic + Nc] * D0[offset + ip][ic];     // 1 additions, 1 multiplications
        w1[0][2][ip] += w[4 + ic + 2 * Nc] * D0[offset + ip][ic]; // 1 additions, 1 multiplications
      }
      for (int ic = 0; ic < Nc; ++ic)
      {
        w1[1][0][ip] += w[4 + ic] * D1[offset + ip][ic];          // 1 additions, 1 multiplications
        w1[1][1][ip] += w[4 + ic + Nc] * D1[offset + ip][ic];     // 1 additions, 1 multiplications
        w1[1][2][ip] += w[4 + ic + 2 * Nc] * D1[offset + ip][ic]; // 1 additions, 1 multiplications
      }
      for (int ic = 0; ic < Nc; ++ic)
      {
        w1[2][0][ip] += w[4 + ic] * D2[offset + ip][ic];          // 1 additions, 1 multiplications
        w1[2][1][ip] += w[4 + ic + Nc] * D2[offset + ip][ic];     // 1 additions, 1 multiplications
        w1[2][2][ip] += w[4 + ic + 2 * Nc] * D2[offset + ip][ic]; // 1 additions, 1 multiplications
      }
      for (int ic = 0; ic < 4; ++ic)
        w0[ip] += w[ic] * Phi[offset + ip][ic]; // 1 additions, 1 multiplications
    }

    T t0[3][3][Np] = {0};
    for (int ip = 0; ip < Np; ip++)
    {
      T temp[3][3] = {{0}};
      // Compute K^T * w1
      temp[0][0] = K[0][0] * w1[0][0][ip] + K[1][0] * w1[1][0][ip] + K[2][0] * w1[2][0][ip]; // 2 additions, 3 multiplications
      temp[0][1] = K[0][0] * w1[0][1][ip] + K[1][0] * w1[1][1][ip] + K[2][0] * w1[2][1][ip]; // 2 additions, 3 multiplications
      temp[0][2] = K[0][0] * w1[0][2][ip] + K[1][0] * w1[1][2][ip] + K[2][0] * w1[2][2][ip]; // 2 additions, 3 multiplications
      temp[1][0] = K[0][1] * w1[0][0][ip] + K[1][1] * w1[1][0][ip] + K[2][1] * w1[2][0][ip]; // 2 additions, 3 multiplications
      temp[1][1] = K[0][1] * w1[0][1][ip] + K[1][1] * w1[1][1][ip] + K[2][1] * w1[2][1][ip]; // 2 additions, 3 multiplications
      temp[1][2] = K[0][1] * w1[0][2][ip] + K[1][1] * w1[1][2][ip] + K[2][1] * w1[2][2][ip]; // 2 additions, 3 multiplications
      temp[2][0] = K[0][2] * w1[0][0][ip] + K[1][2] * w1[1][0][ip] + K[2][2] * w1[2][0][ip]; // 2 additions, 3 multiplications
      temp[2][1] = K[0][2] * w1[0][1][ip] + K[1][2] * w1[1][1][ip] + K[2][2] * w1[2][1][ip]; // 2 additions, 3 multiplications
      temp[2][2] = K[0][2] * w1[0][2][ip] + K[1][2] * w1[1][2][ip] + K[2][2] * w1[2][2][ip]; // 2 additions, 3 multiplications

      t0[0][0][ip] = temp[0][0] + temp[0][0]; // 1 additions
      t0[0][1][ip] = temp[0][1] + temp[1][0]; // 1 additions
      t0[0][2][ip] = temp[0][2] + temp[2][0]; // 1 additions
      t0[1][0][ip] = temp[1][0] + temp[0][1]; // 1 additions
      t0[1][1][ip] = temp[1][1] + temp[1][1]; // 1 additions
      t0[1][2][ip] = temp[1][2] + temp[2][1]; // 1 additions
      t0[2][0][ip] = temp[2][0] + temp[0][2]; // 1 additions
      t0[2][1][ip] = temp[2][1] + temp[1][2]; // 1 additions
      t0[2][2][ip] = temp[2][2] + temp[2][2]; // 1 additions
    }

    T t1[3][3][Np] = {0};
    for (int ip = 0; ip < Np; ip++)
    {
      t1[0][0][ip] = (K[0][0] * t0[0][0][ip] + K[1][0] * t0[1][0][ip] + K[2][0] * t0[2][0][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[0][1][ip] = (K[0][0] * t0[0][1][ip] + K[1][0] * t0[1][1][ip] + K[2][0] * t0[2][1][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[0][2][ip] = (K[0][0] * t0[0][2][ip] + K[1][0] * t0[1][2][ip] + K[2][0] * t0[2][2][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[1][0][ip] = (K[0][1] * t0[0][0][ip] + K[1][1] * t0[1][0][ip] + K[2][1] * t0[2][0][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[1][1][ip] = (K[0][1] * t0[0][1][ip] + K[1][1] * t0[1][1][ip] + K[2][1] * t0[2][1][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[1][2][ip] = (K[0][1] * t0[0][2][ip] + K[1][1] * t0[1][2][ip] + K[2][1] * t0[2][2][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[2][0][ip] = (K[0][2] * t0[0][0][ip] + K[1][2] * t0[1][0][ip] + K[2][2] * t0[2][0][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[2][1][ip] = (K[0][2] * t0[0][1][ip] + K[1][2] * t0[1][1][ip] + K[2][2] * t0[2][1][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
      t1[2][2][ip] = (K[0][2] * t0[0][2][ip] + K[1][2] * t0[1][2][ip] + K[2][2] * t0[2][2][ip]) * weights[offset + ip] * w0[ip]; // 2 additions, 3 multiplications
    }

    T fw1[3][3][Nc] = {0};
    // Np * (9 + 18) additions + Np * 9  multiplications
    for (int ip = 0; ip < Np; ip++)
    {
      for (int ic = 0; ic < Nc; ++ic)
      {
        fw1[0][0][ic] += t1[0][0][ip] * D0[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[0][1][ic] += t1[0][1][ip] * D0[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[0][2][ic] += t1[0][2][ip] * D0[offset + ip][ic]; // 1 additions, 1 multiplications
      }
      for (int ic = 0; ic < Nc; ++ic)
      {
        fw1[1][0][ic] += t1[1][0][ip] * D1[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[1][1][ic] += t1[1][1][ip] * D1[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[1][2][ic] += t1[1][2][ip] * D1[offset + ip][ic]; // 1 additions, 1 multiplications
      }
      for (int ic = 0; ic < Nc; ++ic)
      {
        fw1[2][0][ic] += t1[2][0][ip] * D2[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[2][1][ic] += t1[2][1][ip] * D2[offset + ip][ic]; // 1 additions, 1 multiplications
        fw1[2][2][ic] += t1[2][2][ip] * D2[offset + ip][ic]; // 1 additions, 1 multiplications
      }

      for (int ic = 0; ic < Nc; ++ic)
      {
        A[ic] += fw1[0][0][ic] + fw1[0][0][ic] + fw1[1][0][ic] + fw1[0][1][ic] + fw1[2][0][ic] + fw1[0][2][ic];          // 6 additions
        A[ic + Nc] += fw1[0][1][ic] + fw1[1][0][ic] + fw1[1][1][ic] + fw1[1][1][ic] + fw1[2][1][ic] + fw1[1][2][ic];     // 6 additions
        A[ic + 2 * Nc] += fw1[0][2][ic] + fw1[2][0][ic] + fw1[1][2][ic] + fw1[2][1][ic] + fw1[2][2][ic] + fw1[2][2][ic]; // 6 additions
      }
    }
  }
};
