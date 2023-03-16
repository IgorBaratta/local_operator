
#include <cmath>
#include <cstdint>
#define restrict __restrict__

template <typename T, typename S, int degree>
struct Operator
{
  inline void apply(T *restrict A, const T *restrict w, const T *restrict coordinate_dofs){};
};

template <typename T, typename S>
struct Operator<T, S, 1>
{
  constexpr static std::size_t num_dofs = 4;
  inline void apply(T *restrict A, const T *restrict w, const T *restrict coordinate_dofs)
  {
    // Quadrature rules
    static const S weights_fad[4] = {0.04166666666666666, 0.04166666666666666, 0.04166666666666666, 0.04166666666666666};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const S FE8_C0_Qfad[1][1][4][4] =
        {{{{0.1381966011250091, 0.5854101966249688, 0.138196601125011, 0.138196601125011},
           {0.1381966011250091, 0.138196601125011, 0.585410196624969, 0.138196601125011},
           {0.1381966011250091, 0.1381966011250109, 0.138196601125011, 0.585410196624969},
           {0.585410196624967, 0.1381966011250109, 0.138196601125011, 0.138196601125011}}}};
    static const S FE9_C0_D100_Qfad[1][1][1][4] = {{{{-1.0, 1.0, 0.0, 0.0}}}};
    static const S FE9_C1_D010_Qfad[1][1][1][4] = {{{{-1.0, 0.0, 1.0, 0.0}}}};
    static const S FE9_C2_D001_Qfad[1][1][1][4] = {{{{-1.0, 0.0, 0.0, 1.0}}}};
    // Quadrature loop independent computations for quadrature rule fad
    T J_c0 = {0.0};
    T J_c4 = {0.0};
    T J_c8 = {0.0};
    T J_c5 = {0.0};
    T J_c7 = {0.0};
    T J_c1 = {0.0};
    T J_c6 = {0.0};
    T J_c3 = {0.0};
    T J_c2 = {0.0};
    for (int ic = 0; ic < 4; ++ic)
    {
      J_c0 += coordinate_dofs[ic * 3] * FE9_C0_D100_Qfad[0][0][0][ic];
      J_c4 += coordinate_dofs[ic * 3 + 1] * FE9_C1_D010_Qfad[0][0][0][ic];
      J_c8 += coordinate_dofs[ic * 3 + 2] * FE9_C2_D001_Qfad[0][0][0][ic];
      J_c5 += coordinate_dofs[ic * 3 + 1] * FE9_C2_D001_Qfad[0][0][0][ic];
      J_c7 += coordinate_dofs[ic * 3 + 2] * FE9_C1_D010_Qfad[0][0][0][ic];
      J_c1 += coordinate_dofs[ic * 3] * FE9_C1_D010_Qfad[0][0][0][ic];
      J_c6 += coordinate_dofs[ic * 3 + 2] * FE9_C0_D100_Qfad[0][0][0][ic];
      J_c3 += coordinate_dofs[ic * 3 + 1] * FE9_C0_D100_Qfad[0][0][0][ic];
      J_c2 += coordinate_dofs[ic * 3] * FE9_C2_D001_Qfad[0][0][0][ic];
    }
    T sp_fad[15];
    sp_fad[0] = J_c4 * J_c8;
    sp_fad[1] = J_c5 * J_c7;
    sp_fad[2] = sp_fad[0] + -1 * sp_fad[1];
    sp_fad[3] = J_c0 * sp_fad[2];
    sp_fad[4] = J_c5 * J_c6;
    sp_fad[5] = J_c3 * J_c8;
    sp_fad[6] = sp_fad[4] + -1 * sp_fad[5];
    sp_fad[7] = J_c1 * sp_fad[6];
    sp_fad[8] = sp_fad[3] + sp_fad[7];
    sp_fad[9] = J_c3 * J_c7;
    sp_fad[10] = J_c4 * J_c6;
    sp_fad[11] = sp_fad[9] + -1 * sp_fad[10];
    sp_fad[12] = J_c2 * sp_fad[11];
    sp_fad[13] = sp_fad[8] + sp_fad[12];
    sp_fad[14] = sp_fad[13];
    for (int iq = 0; iq < 4; ++iq)
    {
      // Quadrature loop body setup for quadrature rule fad
      // Varying computations for quadrature rule fad
      T w0 = {0.0};
      for (int ic = 0; ic < 4; ++ic)
        w0 += w[ic] * FE8_C0_Qfad[0][0][iq][ic];
      T sv_fad[1];
      sv_fad[0] = sp_fad[14] * w0;
      const T fw0 = sv_fad[0] * weights_fad[iq];
      for (int i = 0; i < 4; ++i)
        A[i] += fw0 * FE8_C0_Qfad[0][0][iq][i];
    };
  }
};

template <typename T, typename S>
struct Operator<T, S, 2>
{
  constexpr static std::size_t num_dofs = 10;
  inline void apply(T *restrict A, const T *restrict w, const T *restrict coordinate_dofs)
  {
    // Quadrature rules
    static const S weights_d8a[14] = {0.003174603174603167, 0.003174603174603167, 0.003174603174603167, 0.003174603174603167, 0.003174603174603167, 0.003174603174603167, 0.01476497079049678, 0.01476497079049678, 0.01476497079049678, 0.01476497079049678, 0.02213979111426512, 0.02213979111426512, 0.02213979111426512, 0.02213979111426512};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static const S FE11_C0_Qd8a[1][1][14][10] =
        {{{{0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0},
           {0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0},
           {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0},
           {-0.08031550417191768, 0.2771604624527406, -0.08031550417191773, -0.08031550417191774, 0.04042252210657356, 0.2808394945810975, 0.2808394945810974, 0.04042252210657353, 0.04042252210657356, 0.2808394945810972},
           {0.2771604624527405, -0.08031550417191775, -0.08031550417191773, -0.08031550417191774, 0.04042252210657366, 0.04042252210657363, 0.04042252210657363, 0.2808394945810975, 0.2808394945810975, 0.2808394945810974},
           {-0.08031550417191777, -0.08031550417191774, -0.08031550417191773, 0.2771604624527406, 0.2808394945810974, 0.2808394945810973, 0.04042252210657359, 0.2808394945810974, 0.04042252210657358, 0.04042252210657351},
           {-0.08031550417191773, -0.08031550417191773, 0.2771604624527406, -0.08031550417191774, 0.2808394945810974, 0.04042252210657354, 0.2808394945810974, 0.04042252210657354, 0.2808394945810974, 0.04042252210657356},
           {-0.116712266316459, -0.05041039684813053, -0.1167122663164589, -0.116712266316459, 0.3953212143534667, 0.07152785091236928, 0.07152785091236939, 0.3953212143534666, 0.3953212143534667, 0.0715278509123693},
           {-0.05041039684813056, -0.116712266316459, -0.116712266316459, -0.1167122663164589, 0.3953212143534665, 0.3953212143534666, 0.3953212143534665, 0.07152785091236931, 0.07152785091236925, 0.07152785091236927},
           {-0.116712266316459, -0.116712266316459, -0.1167122663164589, -0.05041039684813052, 0.07152785091236931, 0.07152785091236931, 0.3953212143534667, 0.0715278509123693, 0.3953212143534666, 0.3953212143534665},
           {-0.116712266316459, -0.116712266316459, -0.05041039684813054, -0.116712266316459, 0.07152785091236939, 0.3953212143534666, 0.07152785091236936, 0.3953212143534666, 0.07152785091236938, 0.3953212143534666}}}};
    static const S FE9_C0_D100_Qd8a[1][1][1][4] = {{{{-1.0, 1.0, 0.0, 0.0}}}};
    static const S FE9_C1_D010_Qd8a[1][1][1][4] = {{{{-1.0, 0.0, 1.0, 0.0}}}};
    static const S FE9_C2_D001_Qd8a[1][1][1][4] = {{{{-1.0, 0.0, 0.0, 1.0}}}};
    // Quadrature loop independent computations for quadrature rule d8a
    T J_c0 = {0.0};
    T J_c4 = {0.0};
    T J_c8 = {0.0};
    T J_c5 = {0.0};
    T J_c7 = {0.0};
    T J_c1 = {0.0};
    T J_c6 = {0.0};
    T J_c3 = {0.0};
    T J_c2 = {0.0};
    for (int ic = 0; ic < 4; ++ic)
    {
      J_c0 += coordinate_dofs[ic * 3] * FE9_C0_D100_Qd8a[0][0][0][ic];
      J_c4 += coordinate_dofs[ic * 3 + 1] * FE9_C1_D010_Qd8a[0][0][0][ic];
      J_c8 += coordinate_dofs[ic * 3 + 2] * FE9_C2_D001_Qd8a[0][0][0][ic];
      J_c5 += coordinate_dofs[ic * 3 + 1] * FE9_C2_D001_Qd8a[0][0][0][ic];
      J_c7 += coordinate_dofs[ic * 3 + 2] * FE9_C1_D010_Qd8a[0][0][0][ic];
      J_c1 += coordinate_dofs[ic * 3] * FE9_C1_D010_Qd8a[0][0][0][ic];
      J_c6 += coordinate_dofs[ic * 3 + 2] * FE9_C0_D100_Qd8a[0][0][0][ic];
      J_c3 += coordinate_dofs[ic * 3 + 1] * FE9_C0_D100_Qd8a[0][0][0][ic];
      J_c2 += coordinate_dofs[ic * 3] * FE9_C2_D001_Qd8a[0][0][0][ic];
    }
    T sp_d8a[18];
    sp_d8a[0] = J_c4 * J_c8;
    sp_d8a[1] = J_c5 * J_c7;
    sp_d8a[2] = sp_d8a[0] + -1 * sp_d8a[1];
    sp_d8a[3] = J_c0 * sp_d8a[2];
    sp_d8a[4] = J_c5 * J_c6;
    sp_d8a[5] = J_c3 * J_c8;
    sp_d8a[6] = sp_d8a[4] + -1 * sp_d8a[5];
    sp_d8a[7] = J_c1 * sp_d8a[6];
    sp_d8a[8] = sp_d8a[3] + sp_d8a[7];
    sp_d8a[9] = J_c3 * J_c7;
    sp_d8a[10] = J_c4 * J_c6;
    sp_d8a[11] = sp_d8a[9] + -1 * sp_d8a[10];
    sp_d8a[12] = J_c2 * sp_d8a[11];
    sp_d8a[13] = sp_d8a[8] + sp_d8a[12];
    sp_d8a[14] = sp_d8a[13];
    for (int iq = 0; iq < 14; ++iq)
    {
      // Quadrature loop body setup for quadrature rule d8a
      // Varying computations for quadrature rule d8a
      T w0 = {0.0};
      for (int ic = 0; ic < 10; ++ic)
        w0 += w[ic] * FE11_C0_Qd8a[0][0][iq][ic];
      T sv_d8a[1];
      sv_d8a[0] = sp_d8a[14] * w0;
      const T fw0 = sv_d8a[0] * weights_d8a[iq];
      for (int i = 0; i < 10; ++i)
        A[i] += fw0 * FE11_C0_Qd8a[0][0][iq][i];
    }
  }
};