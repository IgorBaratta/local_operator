
#include <cmath>
#include <cstdint>
#include <iostream>
#define restrict __restrict__

template <typename T, typename S, int Nq>
void compute_jacobian(T *restrict jac, const T *restrict coordinate_dofs, const S FE_TF2[Nq][2])
{
    static constexpr S FE_TF3[2] = {-1.0, 1.0};
    T J_c4[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 8];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int id = 0; id < 2 * Nq; ++id)
                temp1[id] += FE_TF3[ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                temp1transp[Nq * ic2 + iq0] = temp1[2 * iq0 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c4[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c8[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 16];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2 * Nq; ++id)
                    temp1[2 * Nq * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[2 * Nq * iq1 + 2 * iq0 + ic2];
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int id = 0; id < Nq * Nq; ++id)
                J_c8[id] += FE_TF3[ic2] * temp1transp[Nq * Nq * ic2 + id];
    }
    T J_c5[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 8];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2 * Nq; ++id)
                    temp1[2 * Nq * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[2 * Nq * iq1 + 2 * iq0 + ic2];
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int id = 0; id < Nq * Nq; ++id)
                J_c5[id] += FE_TF3[ic2] * temp1transp[Nq * Nq * ic2 + id];
    }
    T J_c7[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 16];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int id = 0; id < 2 * Nq; ++id)
                temp1[id] += FE_TF3[ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                temp1transp[Nq * ic2 + iq0] = temp1[2 * iq0 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c7[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c0[Nq * Nq] = {0};
    {
        T temp0[4] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[4 * ic0 + 2 * ic1 + ic2];
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int id = 0; id < 4; ++id)
                temp0[id] += FE_TF3[ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                temp0transp[2 * ic1 + ic2] = temp0[2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2; ++id)
                    temp1[2 * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                temp1transp[Nq * ic2 + iq1] = temp1[2 * iq1 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c0[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c1[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[4 * ic0 + 2 * ic1 + ic2];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int id = 0; id < 2 * Nq; ++id)
                temp1[id] += FE_TF3[ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                temp1transp[Nq * ic2 + iq0] = temp1[2 * iq0 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c1[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c6[Nq * Nq] = {0};
    {
        T temp0[4] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 16];
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int id = 0; id < 4; ++id)
                temp0[id] += FE_TF3[ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                temp0transp[2 * ic1 + ic2] = temp0[2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2; ++id)
                    temp1[2 * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                temp1transp[Nq * ic2 + iq1] = temp1[2 * iq1 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c6[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c3[Nq * Nq] = {0};
    {
        T temp0[4] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[(4 * ic0 + 2 * ic1 + ic2) + 8];
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int id = 0; id < 4; ++id)
                temp0[id] += FE_TF3[ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq] = {0};
        T temp0transp[4] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                temp0transp[2 * ic1 + ic2] = temp0[2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2; ++id)
                    temp1[2 * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * ic1 + id];
        T temp1transp[2 * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                temp1transp[Nq * ic2 + iq1] = temp1[2 * iq1 + ic2];
        for (int iq2 = 0; iq2 < Nq; ++iq2)
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int id = 0; id < Nq; ++id)
                    J_c3[Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * ic2 + id];
    }
    T J_c2[Nq * Nq] = {0};
    {
        T temp0[4 * Nq] = {0};
        T coordinate_dofstransp[8] = {0};
        for (int ic0 = 0; ic0 < 2; ++ic0)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    coordinate_dofstransp[4 * ic0 + 2 * ic1 + ic2] = coordinate_dofs[4 * ic0 + 2 * ic1 + ic2];
        for (int iq0 = 0; iq0 < Nq; ++iq0)
            for (int ic0 = 0; ic0 < 2; ++ic0)
                for (int id = 0; id < 4; ++id)
                    temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * coordinate_dofstransp[4 * ic0 + id];
        T temp1[2 * Nq * Nq] = {0};
        T temp0transp[4 * Nq] = {0};
        for (int ic1 = 0; ic1 < 2; ++ic1)
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
        for (int iq1 = 0; iq1 < Nq; ++iq1)
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int id = 0; id < 2 * Nq; ++id)
                    temp1[2 * Nq * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * Nq * ic1 + id];
        T temp1transp[2 * Nq * Nq] = {0};
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[2 * Nq * iq1 + 2 * iq0 + ic2];
        for (int ic2 = 0; ic2 < 2; ++ic2)
            for (int id = 0; id < Nq * Nq; ++id)
                J_c2[id] += FE_TF3[ic2] * temp1transp[Nq * Nq * ic2 + id];
    }
    constexpr int cubNq = Nq * Nq * Nq;
    for (int iq0 = 0; iq0 < Nq; ++iq0)
    {
        for (int iq1 = 0; iq1 < Nq; ++iq1)
        {
            for (int iq2 = 0; iq2 < Nq; ++iq2)
            {
                int giq = iq0 * Nq * Nq + iq1 * Nq + iq2;
                jac[giq] = J_c0[Nq * iq1 + iq2];
                jac[giq + cubNq] = J_c1[Nq * iq0 + iq2];
                jac[giq + cubNq * 2] = J_c2[Nq * iq0 + iq1];
                jac[giq + cubNq * 3] = J_c3[Nq * iq1 + iq2];
                jac[giq + cubNq * 4] = J_c4[Nq * iq0 + iq2];
                jac[giq + cubNq * 5] = J_c5[Nq * iq0 + iq1];
                jac[giq + cubNq * 6] = J_c6[Nq * iq1 + iq2];
                jac[giq + cubNq * 7] = J_c7[Nq * iq0 + iq2];
                jac[giq + cubNq * 8] = J_c8[Nq * iq0 + iq1];
            }
        }
    }
}

template <int Np, typename T, int cubNq>
static inline void transform_precompute(const T *restrict w1_d, const T *restrict w0, const T *restrict G, T *restrict fw)
{
    for (int ip = 0; ip < Np; ip++)
    {
        // Compute fw = G * w1_d * w0 * detJ
        // FLOPS : 7 * 3 * Np
        for (int ip = 0; ip < Np; ip++)
        {
            fw[ip] = (G[ip] * w1_d[ip] + G[cubNq + ip] * w1_d[cubNq + ip] + G[2 * cubNq + ip] * w1_d[2 * cubNq + ip]) * w0[ip];
            fw[cubNq + ip] = (G[cubNq + ip] * w1_d[ip] + G[3 * cubNq + ip] * w1_d[cubNq + ip] + G[4 * cubNq + ip] * w1_d[2 * cubNq + ip]) * w0[ip];
            fw[2 * cubNq + ip] = (G[2 * cubNq + ip] * w1_d[ip] + G[4 * cubNq + ip] * w1_d[cubNq + ip] + G[5 * cubNq + ip] * w1_d[2 * cubNq + ip]) * w0[ip];
        }
    }
}

template <int Np, typename T, int cubNq>
static inline void transform(const T *restrict jac, const T *restrict w1_d, const T *restrict w0, T *restrict fw)
{

    // Compute determinant of the jacobian
    // FLOPS : 4 * 3 * Np
    T detJ[Np] = {0};
    for (int ip = 0; ip < Np; ip++)
    {
        detJ[ip] = jac[ip] * (jac[4 * cubNq + ip] * jac[8 * cubNq + ip] - jac[7 * cubNq + ip] * jac[5 * cubNq + ip]);
        detJ[ip] -= jac[cubNq + ip] * (jac[3 * cubNq + ip] * jac[8 * cubNq + ip] - jac[6 * cubNq + ip] * jac[5 * cubNq + ip]);
        detJ[ip] += jac[2 * cubNq + ip] * (jac[3 * cubNq + ip] * jac[7 * cubNq + ip] - jac[6 * cubNq + ip] * jac[4 * cubNq + ip]);
    }

    // Compute the inverse of the jacobian
    // FLOPS : 9 * 5 * Np
    T invJ[3][3][Np] = {{{0.0}}};
    for (int ip = 0; ip < Np; ip++)
    {
        invJ[0][0][ip] = (jac[4 * cubNq + ip] * jac[8 * cubNq + ip] - jac[5 * cubNq + ip] * jac[7 * cubNq + ip]) / detJ[ip];
        invJ[0][1][ip] = (jac[2 * cubNq + ip] * jac[7 * cubNq + ip] - jac[cubNq + ip] * jac[8 * cubNq + ip]) / detJ[ip];
        invJ[0][2][ip] = (jac[cubNq + ip] * jac[5 * cubNq + ip] - jac[2 * cubNq + ip] * jac[4 * cubNq + ip]) / detJ[ip];
        invJ[1][0][ip] = (jac[5 * cubNq + ip] * jac[6 * cubNq + ip] - jac[3 * cubNq + ip] * jac[8 * cubNq + ip]) / detJ[ip];
        invJ[1][1][ip] = (jac[ip] * jac[8 * cubNq + ip] - jac[2 * cubNq + ip] * jac[6 * cubNq + ip]) / detJ[ip];
        invJ[1][2][ip] = (jac[2 * cubNq + ip] * jac[3 * cubNq + ip] - jac[ip] * jac[5 * cubNq + ip]) / detJ[ip];
        invJ[2][0][ip] = (jac[3 * cubNq + ip] * jac[7 * cubNq + ip] - jac[4 * cubNq + ip] * jac[6 * cubNq + ip]) / detJ[ip];
        invJ[2][1][ip] = (jac[cubNq + ip] * jac[6 * cubNq + ip] - jac[ip] * jac[7 * cubNq + ip]) / detJ[ip];
        invJ[2][2][ip] = (jac[ip] * jac[4 * cubNq + ip] - jac[cubNq + ip] * jac[3 * cubNq + ip]) / detJ[ip];
    }

    // Compute t = K * w_d
    // FLOPS : 3 * 5 * Np
    T t[3][Np] = {{0}};
    for (int ip = 0; ip < Np; ip++)
    {
        t[0][ip] = invJ[0][0][ip] * w1_d[ip] + invJ[0][1][ip] * w1_d[cubNq + ip] + invJ[0][2][ip] * w1_d[2 * cubNq + ip];
        t[1][ip] = invJ[1][0][ip] * w1_d[ip] + invJ[1][1][ip] * w1_d[cubNq + ip] + invJ[1][2][ip] * w1_d[2 * cubNq + ip];
        t[2][ip] = invJ[2][0][ip] * w1_d[ip] + invJ[2][1][ip] * w1_d[cubNq + ip] + invJ[2][2][ip] * w1_d[2 * cubNq + ip];
    }

    // Compute fw = (K * t) * w0 * detJ
    // FLOPS : 7 * 3 * Np
    for (int ip = 0; ip < Np; ip++)
    {
        fw[ip] = (invJ[0][0][ip] * t[0][ip] + invJ[1][0][ip] * t[1][ip] + invJ[2][0][ip] * t[2][ip]) * w0[ip] * detJ[ip];
        fw[cubNq + ip] = (invJ[0][1][ip] * t[0][ip] + invJ[1][1][ip] * t[1][ip] + invJ[2][1][ip] * t[2][ip]) * w0[ip] * detJ[ip];
        fw[2 * cubNq + ip] = (invJ[0][2][ip] * t[0][ip] + invJ[1][2][ip] * t[1][ip] + invJ[2][2][ip] * t[2][ip]) * w0[ip] * detJ[ip];
    }
}

template <typename T, typename S, int P, int block_size, bool precompute = false>
struct Operator
{
    constexpr static std::size_t num_dofs = (P + 1) * (P + 1) * (P + 1);
    constexpr static int Nq = (P + 3) * (P + 3) * (P + 3);
    inline void apply(T *restrict A, const T *restrict w, const T *restrict geom)
    {
        constexpr int Nq = P + 3;
        constexpr int Nd = P + 1;
        constexpr int cubNq = Nq * Nq * Nq;
        // Quadrature loop independent computations for quadrature rule 037
        // Quadrature loop body setup for quadrature rule 037
        // Varying computations for quadrature rule 037
        T w1_d[3 * Nq * Nq * Nq] = {{0}};
        {
#if BATCH_SIZE == 1
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int id = 0; id < Nd * Nd; ++id)
                        temp0[Nd * Nd * iq0 + id] += FE_TF0[iq0][ic0] * w[8 + Nd * Nd * ic0 + id];
            T temp1[Nq * Nq * Nd] = {0};
            T temp0transp[Nq * Nd * Nd] = {0};
            for (int ic1 = 0; ic1 < Nd; ++ic1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int ic2 = 0; ic2 < Nd; ++ic2)
                        temp0transp[Nd * Nq * ic1 + Nd * iq0 + ic2] = temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int id = 0; id < Nd * Nq; ++id)
                        temp1[Nd * Nq * iq1 + id] += FE_TF1[iq1][ic1] * temp0transp[Nd * Nq * ic1 + id];
            T temp1transp[Nq * Nq * Nd] = {0};
            for (int ic2 = 0; ic2 < Nd; ++ic2)
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int id = 0; id < Nq * Nq; ++id)
                        w1_d[Nq * Nq * iq2 + id] += FE_TF1[iq2][ic2] * temp1transp[Nq * Nq * ic2 + id];
#else
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int ic1 = 0; ic1 < Nd; ++ic1)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp0[Nd * Nd * iq0 + Nd * ic1 + ic2] += FE_TF0[iq0][ic0] * w[8 + (Nd * Nd * ic0 + Nd * ic1 + ic2)];
            T temp1[Nq * Nq * Nd] = {0};
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp1[Nd * Nq * iq1 + Nd * iq0 + ic2] += FE_TF1[iq1][ic1] * temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq0 = 0; iq0 < Nq; ++iq0)
                            w1_d[Nq * Nq * iq0 + Nq * iq1 + iq2] += FE_TF1[iq2][ic2] * temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
#endif
        }
        {
#if BATCH_SIZE == 1
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int id = 0; id < Nd * Nd; ++id)
                        temp0[Nd * Nd * iq0 + id] += FE_TF1[iq0][ic0] * w[8 + Nd * Nd * ic0 + id];
            T temp1[Nq * Nq * Nd] = {0};
            T temp0transp[Nq * Nd * Nd] = {0};
            for (int ic1 = 0; ic1 < Nd; ++ic1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int ic2 = 0; ic2 < Nd; ++ic2)
                        temp0transp[Nd * Nq * ic1 + Nd * iq0 + ic2] = temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int id = 0; id < Nd * Nq; ++id)
                        temp1[Nd * Nq * iq1 + id] += FE_TF0[iq1][ic1] * temp0transp[Nd * Nq * ic1 + id];
            T temp1transp[Nq * Nq * Nd] = {0};
            for (int ic2 = 0; ic2 < Nd; ++ic2)
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int id = 0; id < Nq * Nq; ++id)
                        w1_d[cubNq + Nq * Nq * iq2 + id] += FE_TF1[iq2][ic2] * temp1transp[Nq * Nq * ic2 + id];
#else
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int ic1 = 0; ic1 < Nd; ++ic1)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp0[Nd * Nd * iq0 + Nd * ic1 + ic2] += FE_TF1[iq0][ic0] * w[8 + (Nd * Nd * ic0 + Nd * ic1 + ic2)];
            T temp1[Nq * Nq * Nd] = {0};
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp1[Nd * Nq * iq1 + Nd * iq0 + ic2] += FE_TF0[iq1][ic1] * temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq0 = 0; iq0 < Nq; ++iq0)
                             w1_d[cubNq +Nq * Nq * iq0 + Nq * iq1 + iq2] += FE_TF1[iq2][ic2] * temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
#endif
        }
        {
#if BATCH_SIZE == 1
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int id = 0; id < Nd * Nd; ++id)
                        temp0[Nd * Nd * iq0 + id] += FE_TF1[iq0][ic0] * w[8 + Nd * Nd * ic0 + id];
            T temp1[Nq * Nq * Nd] = {0};
            T temp0transp[Nq * Nd * Nd] = {0};
            for (int ic1 = 0; ic1 < Nd; ++ic1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int ic2 = 0; ic2 < Nd; ++ic2)
                        temp0transp[Nd * Nq * ic1 + Nd * iq0 + ic2] = temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int id = 0; id < Nd * Nq; ++id)
                        temp1[Nd * Nq * iq1 + id] += FE_TF1[iq1][ic1] * temp0transp[Nd * Nq * ic1 + id];
            T temp1transp[Nq * Nq * Nd] = {0};
            for (int ic2 = 0; ic2 < Nd; ++ic2)
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int id = 0; id < Nq * Nq; ++id)
                        w1_d[2 * cubNq + Nq * Nq * iq2 + id] += FE_TF0[iq2][ic2] * temp1transp[Nq * Nq * ic2 + id];
#else
            T temp0[Nq * Nd * Nd] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < Nd; ++ic0)
                    for (int ic1 = 0; ic1 < Nd; ++ic1)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp0[Nd * Nd * iq0 + Nd * ic1 + ic2] += FE_TF1[iq0][ic0] * w[8 + (25 * ic0 + 5 * ic1 + ic2)];
            T temp1[Nq * Nq * Nd] = {0};
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < Nd; ++ic1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        for (int ic2 = 0; ic2 < Nd; ++ic2)
                            temp1[Nd * Nq * iq1 + Nd * iq0 + ic2] += FE_TF1[iq1][ic1] * temp0[Nd * Nd * iq0 + Nd * ic1 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < Nd; ++ic2)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq0 = 0; iq0 < Nq; ++iq0)
                            w1_d[2 * cubNq +Nq * Nq * iq0 + Nq * iq1 + iq2] += FE_TF0[iq2][ic2] * temp1[Nd * Nq * iq1 + Nd * iq0 + ic2];
#endif
        }
        T w0[Nq * Nq * Nq] = {0};
        {
#if BATCH_SIZE == 1
            T temp0[4 * Nq] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < 2; ++ic0)
                    for (int id = 0; id < 4; ++id)
                        temp0[4 * iq0 + id] += FE_TF2[iq0][ic0] * w[4 * ic0 + id];
            T temp1[2 * Nq * Nq] = {0};
            T temp0transp[4 * Nq] = {0};
            for (int ic1 = 0; ic1 < 2; ++ic1)
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int ic2 = 0; ic2 < 2; ++ic2)
                        temp0transp[2 * Nq * ic1 + 2 * iq0 + ic2] = temp0[4 * iq0 + 2 * ic1 + ic2];
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < 2; ++ic1)
                    for (int id = 0; id < 2 * Nq; ++id)
                        temp1[2 * Nq * iq1 + id] += FE_TF2[iq1][ic1] * temp0transp[2 * Nq * ic1 + id];
            T temp1transp[2 * Nq * Nq] = {0};
            for (int ic2 = 0; ic2 < 2; ++ic2)
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        temp1transp[Nq * Nq * ic2 + Nq * iq1 + iq0] = temp1[2 * Nq * iq1 + 2 * iq0 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    for (int id = 0; id < Nq * Nq; ++id)
                        w0[Nq * Nq * iq2 + id] += FE_TF2[iq2][ic2] * temp1transp[Nq * Nq * ic2 + id];
#else
            T temp0[4 * Nq] = {0};
            for (int iq0 = 0; iq0 < Nq; ++iq0)
                for (int ic0 = 0; ic0 < 2; ++ic0)
                    for (int ic1 = 0; ic1 < 2; ++ic1)
                        for (int ic2 = 0; ic2 < 2; ++ic2)
                            temp0[4 * iq0 + 2 * ic1 + ic2] += FE_TF2[iq0][ic0] * w[4 * ic0 + 2 * ic1 + ic2];
            T temp1[2 * Nq * Nq] = {0};
            for (int iq1 = 0; iq1 < Nq; ++iq1)
                for (int ic1 = 0; ic1 < 2; ++ic1)
                    for (int iq0 = 0; iq0 < Nq; ++iq0)
                        for (int ic2 = 0; ic2 < 2; ++ic2)
                            temp1[2 * Nq * iq1 + 2 * iq0 + ic2] += FE_TF2[iq1][ic1] * temp0[4 * iq0 + 2 * ic1 + ic2];
            for (int iq2 = 0; iq2 < Nq; ++iq2)
                for (int ic2 = 0; ic2 < 2; ++ic2)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq0 = 0; iq0 < Nq; ++iq0)
                            w0[Nq * Nq * iq0 + Nq * iq1 + iq2] += FE_TF2[iq2][ic2] * temp1[2 * Nq * iq1 + 2 * iq0 + ic2];
#endif
        }

        // Transform K^TK * detJ * grad(u) * weights
        // -------------------------------------------------------------------------------------------
        constexpr int num_blocks = cubNq / block_size;
        constexpr int remainder = cubNq % block_size;
        T fw[3 * cubNq] = {0};

        // Geometric factor G is available, use transform_precompute
        if constexpr (precompute == true)
        {
            for (int i = 0; i < num_blocks; i++)
            {
                int offset = block_size * i;
                transform_precompute<block_size, T, cubNq>(w1_d + offset, w0 + offset, geom + offset, fw + offset);
            }
            if constexpr (remainder > 0)
            {
                int offset = block_size * num_blocks;
                transform_precompute<remainder, T, cubNq>(w1_d + offset, w0 + offset, geom + offset, fw + offset);
            }
        }
        // Compute jacobian on the fly
        else
        {
            T jac[9 * cubNq] = {0};
            compute_jacobian<T, S, Nq>(jac, geom, FE_TF2);
            for (int i = 0; i < num_blocks; i++)
            {
                int offset = block_size * i;
                transform<block_size, T, cubNq>(jac + offset, w1_d + offset, w0 + offset, fw + offset);
            }
            if constexpr (remainder > 0)
            {
                int offset = block_size * num_blocks;
                transform<remainder, T, cubNq>(jac + offset, w1_d + offset, w0 + offset, fw + offset);
            }
        }
        // -------------------------------------------------------------------------------------------
        {
            {
                T temp0[Nq * Nq * Nd] = {0};
                T fw0transp[Nq * Nq * Nq] = {0};
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            fw0transp[Nq * Nq * iq0 + Nq * iq1 + iq2] = fw[Nq * Nq * iq0 + Nq * iq1 + iq2] * (weights[iq0] * weights[iq1] * weights[iq2]);
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int id = 0; id < Nq * Nq; ++id)
                            temp0[Nq * Nq * i0 + id] += FE_TF0[iq0][i0] * fw0transp[Nq * Nq * iq0 + id];
                T temp1[Nq * Nd * Nd] = {0};
                T temp0transp[Nq * Nq * Nd] = {0};
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            temp0transp[Nd * Nq * iq1 + Nq * i0 + iq2] = temp0[Nq * Nq * i0 + Nq * iq1 + iq2];
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int id = 0; id < Nd * Nq; ++id)
                            temp1[Nd * Nq * i1 + id] += FE_TF1[iq1][i1] * temp0transp[Nd * Nq * iq1 + id];
                T temp1transp[Nq * Nd * Nd] = {0};
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int i0 = 0; i0 < Nd; ++i0)
                            temp1transp[Nd * Nd * iq2 + Nd * i1 + i0] = temp1[Nd * Nq * i1 + Nq * i0 + iq2];
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i2 = 0; i2 < Nd; ++i2)
                        for (int id = 0; id < Nd * Nd; ++id)
                            A[Nd * Nd * i2 + id] += FE_TF1[iq2][i2] * temp1transp[Nd * Nd * iq2 + id];
            }
            {
                T temp0[Nq * Nq * Nd] = {0};
                T fw1transp[Nq * Nq * Nq] = {0};
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            fw1transp[Nq * Nq * iq0 + Nq * iq1 + iq2] = fw[cubNq + Nq * Nq * iq0 + Nq * iq1 + iq2] * (weights[iq0] * weights[iq1] * weights[iq2]);
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int id = 0; id < Nq * Nq; ++id)
                            temp0[Nq * Nq * i0 + id] += FE_TF1[iq0][i0] * fw1transp[Nq * Nq * iq0 + id];
                T temp1[Nq * Nd * Nd] = {0};
                T temp0transp[Nq * Nq * Nd] = {0};
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            temp0transp[Nd * Nq * iq1 + Nq * i0 + iq2] = temp0[Nq * Nq * i0 + Nq * iq1 + iq2];
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int id = 0; id < Nd * Nq; ++id)
                            temp1[Nd * Nq * i1 + id] += FE_TF0[iq1][i1] * temp0transp[Nd * Nq * iq1 + id];
                T temp1transp[Nq * Nd * Nd] = {0};
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int i0 = 0; i0 < Nd; ++i0)
                            temp1transp[Nd * Nd * iq2 + Nd * i1 + i0] = temp1[Nd * Nq * i1 + Nq * i0 + iq2];
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i2 = 0; i2 < Nd; ++i2)
                        for (int id = 0; id < Nd * Nd; ++id)
                            A[Nd * Nd * i2 + id] += FE_TF1[iq2][i2] * temp1transp[Nd * Nd * iq2 + id];
            }
            {
                T temp0[Nq * Nq * Nd] = {0};
                T fw2transp[Nq * Nq * Nq] = {0};
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int iq1 = 0; iq1 < Nq; ++iq1)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            fw2transp[Nq * Nq * iq0 + Nq * iq1 + iq2] = fw[2 * cubNq + Nq * Nq * iq0 + Nq * iq1 + iq2] * (weights[iq0] * weights[iq1] * weights[iq2]);
                for (int iq0 = 0; iq0 < Nq; ++iq0)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int id = 0; id < Nq * Nq; ++id)
                            temp0[Nq * Nq * i0 + id] += FE_TF1[iq0][i0] * fw2transp[Nq * Nq * iq0 + id];
                T temp1[Nq * Nd * Nd] = {0};
                T temp0transp[Nq * Nq * Nd] = {0};
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i0 = 0; i0 < Nd; ++i0)
                        for (int iq2 = 0; iq2 < Nq; ++iq2)
                            temp0transp[Nd * Nq * iq1 + Nq * i0 + iq2] = temp0[Nq * Nq * i0 + Nq * iq1 + iq2];
                for (int iq1 = 0; iq1 < Nq; ++iq1)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int id = 0; id < Nd * Nq; ++id)
                            temp1[Nd * Nq * i1 + id] += FE_TF1[iq1][i1] * temp0transp[Nd * Nq * iq1 + id];
                T temp1transp[Nq * Nd * Nd] = {0};
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i1 = 0; i1 < Nd; ++i1)
                        for (int i0 = 0; i0 < Nd; ++i0)
                            temp1transp[Nd * Nd * iq2 + Nd * i1 + i0] = temp1[Nd * Nq * i1 + Nq * i0 + iq2];
                for (int iq2 = 0; iq2 < Nq; ++iq2)
                    for (int i2 = 0; i2 < Nd; ++i2)
                        for (int id = 0; id < Nd * Nd; ++id)
                            A[Nd * Nd * i2 + id] += FE_TF0[iq2][i2] * temp1transp[Nd * Nd * iq2 + id];
            }
        }
    }
#if DEGREE == 1
    // Quadrature rules
    static constexpr S weights[4] = {0.1739274225687268, 0.3260725774312731, 0.326072577431273, 0.1739274225687268};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[4][2] =
        {{-1.0, 1.0},
         {-1.0, 1.0},
         {-1.0, 1.0},
         {-1.0, 1.0}};
    static constexpr S FE_TF1[4][2] =
        {{0.9305681557970263, 0.06943184420297366},
         {0.6699905217924281, 0.3300094782075719},
         {0.3300094782075719, 0.6699905217924281},
         {0.06943184420297371, 0.9305681557970262}};
    static constexpr S FE_TF2[4][2] =
        {{0.9305681557970263, 0.06943184420297366},
         {0.6699905217924281, 0.3300094782075719},
         {0.3300094782075719, 0.6699905217924281},
         {0.06943184420297371, 0.9305681557970262}};
#elif DEGREE == 2
    static constexpr S weights[5] = {0.1184634425280946, 0.2393143352496833, 0.2844444444444444, 0.2393143352496833, 0.1184634425280946};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[5][3] =
        {{-2.812359691877328, -0.8123596918773281, 3.624719383754655},
         {-2.076938620211366, -0.07693862021136644, 2.153877240422732},
         {-1.0, 1.0, 0.0},
         {0.07693862021136599, 2.076938620211366, -2.153877240422732},
         {0.8123596918773279, 2.812359691877328, -3.624719383754655}};
    static constexpr S FE_TF1[5][3] =
        {{0.8636708795620425, -0.0425089663766216, 0.1788380868145792},
         {0.4142092540156867, -0.1242600560899964, 0.7100508020743097},
         {-5.551115123125783e-17, 0.0, 1.0},
         {-0.1242600560899964, 0.4142092540156866, 0.7100508020743098},
         {-0.04250896637662177, 0.8636708795620422, 0.1788380868145795}};
    static constexpr S FE_TF2[5][2] =
        {{0.953089922969332, 0.04691007703066796},
         {0.7692346550528416, 0.2307653449471584},
         {0.5, 0.5},
         {0.2307653449471585, 0.7692346550528415},
         {0.04691007703066802, 0.9530899229693319}};
#elif DEGREE == 3
    static constexpr S weights[6] = {0.08566224618958498, 0.1803807865240693, 0.2339569672863456, 0.2339569672863455, 0.1803807865240693, 0.08566224618958498};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[6][4] =
        {{-5.341796516451364, 0.679448945435603, 6.82706288900378, -2.164715317988019},
         {-3.042515413982264, -0.263531518349059, 2.523953938172386, 0.7820929941589372},
         {-0.5600696500842702, -0.6330262803317145, -1.721088004613099, 2.914183935029083},
         {0.633026280331715, 0.5600696500842705, -2.914183935029084, 1.721088004613098},
         {0.2635315183490586, 3.042515413982264, -0.7820929941589359, -2.523953938172387},
         {-0.6794489454356034, 5.341796516451362, 2.16471531798802, -6.827062889003779}};
    static constexpr S FE_TF1[6][4] =
        {{0.8086169815356241, 0.02825726211220378, 0.2516267987615972, -0.0885010424094253},
         {0.2462720621498502, 0.05022525378901735, 0.8718189146469666, -0.1683162305858341},
         {-0.110748722338305, -0.06807738270293576, 0.9039054209873038, 0.274920684053937},
         {-0.06807738270293578, -0.1107487223383049, 0.2749206840539369, 0.9039054209873036},
         {0.0502252537890174, 0.2462720621498504, -0.1683162305858343, 0.8718189146469665},
         {0.02825726211220386, 0.8086169815356242, -0.08850104240942519, 0.2516267987615973}};
    static constexpr S FE_TF2[6][2] =
        {{0.966234757101576, 0.03376524289842397},
         {0.8306046932331324, 0.1693953067668676},
         {0.6193095930415985, 0.3806904069584015},
         {0.3806904069584015, 0.6193095930415985},
         {0.1693953067668676, 0.8306046932331324},
         {0.03376524289842397, 0.966234757101576}};
#elif DEGREE == 4
    // Quadrature rules
    static constexpr S weights[7] = {0.06474248308443456, 0.1398526957446383, 0.1909150252525594, 0.2089795918367347, 0.1909150252525595, 0.1398526957446383, 0.06474248308443487};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[7][5] =
        {{-8.540302315485206, -0.5818411081330925, 10.76418032131331, -3.304517382958648, 1.662480485263638},
         {-3.878725323521341, 0.3948939182513165, 2.340496486173071, 2.275849191900984, -1.13251427280403},
         {0.02611062162665612, 0.2555686340403505, -3.574843617907995, 4.163460420411717, -0.8702960581707269},
         {0.7499999999999998, -0.7499999999999989, -2.673169155390906, -1.490474506377031e-15, 2.673169155390906},
         {-0.2555686340403516, -0.02611062162665467, 0.8702960581707275, -4.163460420411718, 3.574843617907995},
         {-0.3948939182513161, 3.878725323521341, 1.13251427280403, -2.275849191900981, -2.340496486173072},
         {0.5818411081330948, 8.540302315485192, -1.662480485263634, 3.304517382958645, -10.76418032131329}};
    static constexpr S FE_TF1[7][5] =
        {{0.764393793728545, -0.0199586670952077, 0.3082665410350864, -0.1093000994560194, 0.05659843178759583},
         {0.1370626239500496, -0.02034210710893031, 0.9514732616641501, -0.1273991423835425, 0.05920536387827328},
         {-0.1317289831653047, 0.0556728555595951, 0.7339940232440429, 0.514267827100455, -0.1722057227387883},
         {6.938893903907228e-17, -1.387778780781446e-17, -8.326672684688674e-17, 1.0, 1.942890293094024e-16},
         {0.05567285555959507, -0.1317289831653044, -0.1722057227387883, 0.5142678271004548, 0.7339940232440429},
         {-0.02034210710893049, 0.13706262395005, 0.05920536387827337, -0.127399142383543, 0.95147326166415},
         {-0.0199586670952077, 0.7643937937285441, 0.05659843178759594, -0.1093000994560198, 0.3082665410350875}};
    static constexpr S FE_TF2[7][2] =
        {{0.9745539561713794, 0.02544604382862065},
         {0.8707655927996972, 0.1292344072003028},
         {0.7029225756886985, 0.2970774243113014},
         {0.5, 0.5},
         {0.2970774243113015, 0.7029225756886985},
         {0.1292344072003028, 0.8707655927996972},
         {0.02544604382862081, 0.9745539561713792}};
#elif DEGREE == 5
    // Quadrature rules
    static constexpr S weights[8] = {0.05061426814518795, 0.1111905172266872, 0.1568533229389436, 0.181341891689181, 0.181341891689181, 0.1568533229389437, 0.1111905172266872, 0.05061426814518814};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[8][6] =
        {{-12.38195163760049, 0.5076902511042709, 15.40317094369366, -4.493397485457833, 2.344421303574084, -1.379933375313704},
         {-4.600639828588789, -0.4408464845138189, 1.661591997613039, 4.261197300693675, -2.086369524150141, 1.205066538946033},
         {0.6391752981242027, -0.008468397770068137, -5.370045170243436, 4.97155281710012, -0.2812752753932823, 0.04906072818246332},
         {0.6140057438275459, 0.5404193351472224, -1.918641053180408, -2.060550618022463, 4.466466249279335, -1.641699657051235},
         {-0.5404193351472226, -0.6140057438275457, 1.641699657051236, -4.466466249279337, 2.060550618022463, 1.918641053180406},
         {0.008468397770067693, -0.6391752981242039, -0.04906072818246221, 0.2812752753932815, -4.971552817100115, 5.370045170243435},
         {0.4408464845138187, 4.600639828588788, -1.205066538946034, 2.086369524150142, -4.261197300693674, -1.661591997613041},
         {-0.5076902511042705, 12.38195163760048, 1.379933375313702, -2.344421303574089, 4.49339748545782, -15.40317094369365}};
    static constexpr S FE_TF1[8][6] =
        {{0.7286932192353579, 0.01476134369067711, 0.3531455605723022, -0.123663513318777, 0.06702420996760614, -0.03996082014716595},
         {0.06441507565198347, 0.007290025391268107, 0.9872345947467757, -0.07388270651774483, 0.03492584267165069, -0.01998283194393266},
         {-0.1205882945854509, -0.0375050947147014, 0.5691527543107485, 0.6868992194178709, -0.2035888509000691, 0.1056302664716023},
         {0.04088031273045188, 0.02820725875426095, -0.1367506889625262, 0.9460354747981565, 0.2054841216266144, -0.08385647894695758},
         {0.02820725875426093, 0.04088031273045184, -0.08385647894695764, 0.2054841216266144, 0.9460354747981562, -0.136750688962526},
         {-0.03750509471470145, -0.1205882945854509, 0.1056302664716023, -0.2035888509000693, 0.6868992194178714, 0.5691527543107481},
         {0.007290025391268024, 0.06441507565198337, -0.01998283194393277, 0.03492584267165091, -0.07388270651774463, 0.9872345947467756},
         {0.0147613436906771, 0.7286932192353567, -0.03996082014716615, 0.06702420996760625, -0.1236635133187775, 0.3531455605723036}};
    static constexpr S FE_TF2[8][2] =
        {{0.9801449282487682, 0.01985507175123186},
         {0.8983332387068134, 0.1016667612931866},
         {0.7627662049581646, 0.2372337950418354},
         {0.5917173212478249, 0.4082826787521751},
         {0.4082826787521752, 0.5917173212478248},
         {0.2372337950418355, 0.7627662049581645},
         {0.1016667612931866, 0.8983332387068134},
         {0.01985507175123191, 0.9801449282487681}};
#elif DEGREE == 6
    // Quadrature rules
    static constexpr S weights[9] = {0.04063719418078724, 0.09032408034742866, 0.1303053482014677, 0.1561735385200014, 0.1651196775006299, 0.1561735385200015, 0.1303053482014677, 0.09032408034742866, 0.04063719418078724};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[9][7] =
        {{-16.851751124832, -0.449707565532008, 20.72043768445436, -5.739062983810239, 3.002720780829048, -1.872147702689279, 1.189510911580122},
         {-5.226151890706513, 0.4493111974832293, 0.5355703304835957, 6.679017437266775, -3.145580740472883, 1.90115742167227, -1.193323755726469},
         {1.234143950069098, -0.1421811075118372, -7.045102048728367, 5.362891745969385, 0.7893623299089154, -0.5728347939418792, 0.3737199242346838},
         {0.3471790150208963, -0.3124800835384357, -0.9282374495455421, -4.190783329034813, 5.926501172573628, -1.724988454182371, 0.8828091287066384},
         {-0.6249999999999989, 0.6249999999999993, 1.815088942537641, -4.013938481177506, -1.802653698545193e-15, 4.013938481177506, -1.81508894253764},
         {0.3124800835384342, -0.347179015020894, -0.8828091287066335, 1.724988454182359, -5.92650117257363, 4.19078332903483, 0.9282374495455328},
         {0.142181107511836, -1.2341439500691, -0.3737199242346849, 0.572834793941877, -0.7893623299089136, -5.362891745969382, 7.045102048728366},
         {-0.4493111974832298, 5.226151890706515, 1.193323755726468, -1.901157421672267, 3.145580740472885, -6.679017437266779, -0.5355703304835963},
         {0.4497075655320071, 16.85175112483201, -1.189510911580115, 1.872147702689272, -3.002720780829046, 5.739062983810237, -20.72043768445437}};
    static constexpr S FE_TF1[9][7] =
        {{0.6995098678471076, -0.01131626694169769, 0.3893122348117268, -0.1343122465329731, 0.07361500956703126, -0.04666890777843566, 0.02986030902724049},
         {0.01467769306253675, -0.001310808443745745, 0.9992284227915847, -0.01973608507271843, 0.009211848917572513, -0.005553569626622429, 0.003482498371392639},
         {-0.09993398920103827, 0.02394819586252237, 0.429591744672089, 0.804997956212528, -0.2015735745589442, 0.1075016215674345, -0.0645319545545915},
         {0.05932582155455584, -0.03027307319788081, -0.1910361535956981, 0.834827459597666, 0.3956334686460113, -0.1522025571158339, 0.08372503411117982},
         {1.179611963664229e-16, -1.249000902703301e-16, -3.747002708109903e-16, 9.159339953157541e-16, 1.0, -8.049116928532385e-16, 4.302114220422482e-16},
         {-0.03027307319788063, 0.05932582155455578, 0.08372503411117961, -0.1522025571158334, 0.3956334686460117, 0.8348274595976648, -0.1910361535956975},
         {0.02394819586252244, -0.09993398920103848, -0.06453195455459161, 0.1075016215674344, -0.2015735745589447, 0.8049979562125287, 0.4295917446720894},
         {-0.001310808443745801, 0.01467769306253644, 0.003482498371392306, -0.005553569626622401, 0.009211848917572235, -0.01973608507271798, 0.9992284227915851},
         {-0.01131626694169775, 0.6995098678471076, 0.02986030902724035, -0.04666890777843569, 0.07361500956703104, -0.134312246532973, 0.3893122348117271}};
    static constexpr S FE_TF2[9][2] =
        {{0.984080119753813, 0.01591988024618696},
         {0.9180155536633179, 0.08198444633668206},
         {0.8066857163502952, 0.1933142836497048},
         {0.6621267117019045, 0.3378732882980955},
         {0.5, 0.5},
         {0.3378732882980957, 0.6621267117019043},
         {0.1933142836497048, 0.8066857163502952},
         {0.08198444633668206, 0.9180155536633179},
         {0.01591988024618696, 0.984080119753813}};
#elif DEGREE == 7
    // Quadrature rules
    static constexpr S weights[10] = {0.03333567215434398, 0.07472567457529036, 0.1095431812579911, 0.1346333596549982, 0.1477621123573765, 0.1477621123573765, 0.1346333596549983, 0.1095431812579911, 0.07472567457529028, 0.03333567215434398};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[10][8] =
        {{-21.94050741531485, 0.4032566087993779, 26.70072901270923, -7.037657692209953, 3.663353995167439, -2.320484551147652, 1.580493864293892, -1.049183822297493},
         {-5.770027641063378, -0.4408455073440927, -1.001048737074382, 9.489250414310305, -4.289025858900387, 2.607599632497366, -1.746066208236185, 1.150163905810748},
         {1.796012755457883, 0.2319213705943186, -8.589004509171437, 5.370937463941692, 2.2602810842648, -1.387788835681038, 0.923675133193147, -0.6060344625993661},
         {0.01805716514491984, 0.1321954352873847, 0.1516059092724431, -6.269265585335411, 7.029396953468661, -1.332320843350063, 0.632267812516847, -0.3619368470047817},
         {-0.5725670177505662, -0.4852807370033803, 1.613445364674916, -3.141621121401739, -2.269937311358189, 5.915568516390706, -2.395385347301486, 1.33577765374974},
         {0.4852807370033787, 0.5725670177505658, -1.335777653749735, 2.39538534730148, -5.915568516390707, 2.269937311358199, 3.141621121401734, -1.613445364674917},
         {-0.1321954352873815, -0.01805716514491207, 0.3619368470047732, -0.6322678125168357, 1.332320843350042, -7.029396953468651, 6.269265585335428, -0.1516059092724666},
         {-0.2319213705943182, -1.796012755457883, 0.6060344625993654, -0.9236751331931456, 1.387788835681035, -2.260281084264804, -5.370937463941692, 8.589004509171442},
         {0.4408455073440903, 5.770027641063389, -1.150163905810742, 1.746066208236179, -2.607599632497358, 4.289025858900376, -9.489250414310307, 1.001048737074371},
         {-0.403256608799377, 21.94050741531485, 1.049183822297485, -1.58049386429388, 2.320484551147665, -3.663353995167453, 7.037657692209965, -26.70072901270925}};
    static constexpr S FE_TF1[10][8] =
        {{0.6753202332849628, 0.008927195383580618, 0.4189375202774577, -0.1424879209484241, 0.07831996363276503, -0.05061166829344471, 0.03478510115163191, -0.02319042448852887},
         {-0.02033697832508582, -0.001471372735060462, 0.9983084589880243, 0.03102491432481169, -0.0142212550760154, 0.008680299708264033, -0.00582185723486342, 0.003837790349924944},
         {-0.07821578661000622, -0.01493098126178387, 0.3166739008401074, 0.8835542489150553, -0.1812652107581693, 0.09588591118632797, -0.06096719929291208, 0.03926511698138073},
         {0.0642365592956123, 0.02539196826727205, -0.2016797050153993, 0.7105630890353563, 0.5519468065202021, -0.1924540310533169, 0.1097317570744844, -0.06773644412421119},
         {-0.02016277465916416, -0.0149372775833406, 0.05766362452581389, -0.1197692540051789, 0.9651559489434688, 0.1628245029085569, -0.07161591351391355, 0.04084114338375749},
         {-0.01493727758334067, -0.02016277465916435, 0.04084114338375766, -0.07161591351391386, 0.1628245029085581, 0.9651559489434687, -0.1197692540051798, 0.05766362452581446},
         {0.02539196826727205, 0.06423655929561239, -0.06773644412421105, 0.1097317570744843, -0.192454031053317, 0.5519468065202046, 0.7105630890353543, -0.2016797050153997},
         {-0.0149309812617837, -0.07821578661000582, 0.03926511698138041, -0.06096719929291169, 0.09588591118632708, -0.1812652107581682, 0.8835542489150556, 0.316673900840106},
         {-0.001471372735060407, -0.02033697832508528, 0.003837790349924902, -0.005821857234863337, 0.008680299708263312, -0.01422125507601452, 0.03102491432481083, 0.9983084589880243},
         {0.008927195383580569, 0.6753202332849628, -0.02319042448852868, 0.03478510115163183, -0.05061166829344463, 0.07831996363276489, -0.1424879209484245, 0.4189375202774577}};
    static constexpr S FE_TF2[10][2] =
        {{0.9869532642585859, 0.01304673574141413},
         {0.9325316833444923, 0.06746831665550773},
         {0.8397047841495122, 0.1602952158504877},
         {0.7166976970646236, 0.2833023029353764},
         {0.5744371694908156, 0.4255628305091844},
         {0.4255628305091844, 0.5744371694908156},
         {0.2833023029353766, 0.7166976970646234},
         {0.1602952158504877, 0.8397047841495122},
         {0.06746831665550762, 0.9325316833444923},
         {0.01304673574141413, 0.9869532642585859}};
#elif DEGREE == 8
    // Quadrature rules
    static constexpr S weights[11] = {0.02783428355808692, 0.06279018473245228, 0.09314510546386703, 0.1165968822959953, 0.1314022722551234, 0.1364625433889503, 0.1314022722551234, 0.1165968822959953, 0.09314510546386703, 0.06279018473245228, 0.02783428355808692};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[11][9] =
        {{-27.64230198355163, -0.3652806081099778, 33.3340437333771, -8.385005714893705, 4.332084633075659, -2.756596731266882, 1.920503476416588, -1.377639840260324, 0.9401930352131735},
         {-6.24348898266625, 0.4248983426926981, -2.921608983764842, 12.6640486554475, -5.500091437043726, 3.331041246232262, -2.272797510044438, 1.613720790888854, -1.095722121742056},
         {2.321666704350748, -0.2841370420184911, -10.00875366655493, 5.029621967324916, 4.079641027362844, -2.337443722178846, 1.553722960158015, -1.088836551217006, 0.7345183227727501},
         {-0.3365067938754358, -0.0007893978986616146, 1.246752245206748, -8.246266304567737, 7.783326529969367, -0.5404746770857016, 0.1180811550117271, -0.03095006793087096, 0.006827311170562123},
         {-0.4336151960506915, 0.3293594097505661, 1.187093674209232, -2.049617626062589, -4.579456384292646, 7.53739787152795, -2.542121771297293, 1.430428737228688, -0.8794687150132176},
         {0.5468750000000004, -0.5468750000000008, -1.483564795832508, 2.538826172716299, -5.318620435147836, -2.456047588106909e-16, 5.318620435147832, -2.538826172716297, 1.48356479583251},
         {-0.3293594097505661, 0.4336151960506927, 0.8794687150132191, -1.430428737228692, 2.542121771297303, -7.537397871527951, 4.579456384292644, 2.049617626062593, -1.187093674209238},
         {0.000789397898658839, 0.3365067938754446, -0.006827311170556127, 0.03095006793086075, -0.1180811550117111, 0.5404746770856763, -7.783326529969338, 8.246266304567737, -1.246752245206773},
         {0.2841370420184908, -2.321666704350749, -0.7345183227727485, 1.088836551217012, -1.553722960158021, 2.33744372217885, -4.07964102736285, -5.029621967324919, 10.00875366655493},
         {-0.424898342692697, 6.243488982666269, 1.095722121742053, -1.613720790888854, 2.272797510044438, -3.331041246232259, 5.500091437043707, -12.66404865544748, 2.921608983764821},
         {0.3652806081099751, 27.6423019835516, -0.940193035213186, 1.377639840260322, -1.920503476416595, 2.756596731266871, -4.332084633075624, 8.385005714893669, -33.33404373337704}};
    static constexpr S FE_TF1[11][9] =
        {{0.6550006942915133, -0.007208592379485027, 0.443571669492679, -0.1489197736246798, 0.08187124424843992, -0.05331243279483577, 0.03754431313825726, -0.02708152612254702, 0.01853440375065801},
         {-0.04561330210722985, 0.002729876452755686, 0.9904378374334056, 0.077164741437142, -0.0347217455997732, 0.02123816331526446, -0.01455173558448117, 0.01035324506283436, -0.007037080409918003},
         {-0.05834800558352621, 0.009100409810665525, 0.2265938863115029, 0.9345516899095445, -0.1514943153167289, 0.07886310158332813, -0.05085996468024251, 0.03517221788418264, -0.02357901991872619},
         {0.06191326584222909, -0.01960002966535246, -0.1909179563084406, 0.5920934000867238, 0.6741163330778478, -0.2097664945321857, 0.1191862909108875, -0.07824585058808896, 0.05122104117637962},
         {-0.03218253714354972, 0.0185168612859713, 0.09104799364136065, -0.1812958102869983, 0.8871892731048826, 0.3189542974593472, -0.1312206253851247, 0.07806241955209971, -0.04907187222798869},
         {1.387778780781446e-17, -4.163336342344337e-17, -5.551115123125783e-17, 9.71445146547012e-17, -3.747002708109903e-16, 1.0, 1.665334536937735e-16, -1.804112415015879e-16, 1.52655665885959e-16},
         {0.01851686128597131, -0.03218253714354977, -0.04907187222798857, 0.0780624195520998, -0.1312206253851249, 0.3189542974593471, 0.8871892731048824, -0.1812958102869979, 0.09104799364136069},
         {-0.01960002966535247, 0.06191326584222924, 0.05122104117637975, -0.07824585058808932, 0.1191862909108884, -0.2097664945321865, 0.6741163330778506, 0.5920934000867217, -0.1909179563084413},
         {0.009100409810665448, -0.05834800558352592, -0.02357901991872595, 0.03517221788418254, -0.05085996468024199, 0.07886310158332768, -0.1514943153167277, 0.9345516899095441, 0.2265938863115017},
         {0.00272987645275561, -0.04561330210722939, -0.007037080409917712, 0.01035324506283455, -0.01455173558448075, 0.02123816331526424, -0.03472174559977309, 0.07716474143714078, 0.9904378374334056},
         {-0.007208592379485061, 0.6550006942915124, 0.01853440375065786, -0.02708152612254695, 0.03754431313825737, -0.05331243279483583, 0.08187124424844022, -0.1489197736246806, 0.4435716694926806}};
    static constexpr S FE_TF2[11][2] =
        {{0.9891143290730285, 0.01088567092697146},
         {0.9435312998840477, 0.05646870011595229},
         {0.8650760027870247, 0.1349239972129753},
         {0.759548064603406, 0.240451935396594},
         {0.6347715779761725, 0.3652284220238275},
         {0.5, 0.5},
         {0.3652284220238275, 0.6347715779761725},
         {0.2404519353965944, 0.7595480646034056},
         {0.1349239972129753, 0.8650760027870247},
         {0.05646870011595229, 0.9435312998840477},
         {0.01088567092697151, 0.9891143290730284}};
#elif DEGREE == 9
    // Quadrature rules
    static constexpr S weights[12] = {0.02358766819325592, 0.05346966299765923, 0.08003916427167307, 0.101583713361533, 0.1167462682691774, 0.1245735229067014, 0.1245735229067014, 0.1167462682691774, 0.1015837133615331, 0.08003916427167307, 0.05346966299765918, 0.02358766819325592};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[12][10] =
        {{-33.95316647156921, 0.3336961594984826, 40.6136305874231, -9.777856994581068, 5.01030856708002, -3.189661480152751, 2.242493058697686, -1.652980387979738, 1.226078925089393, -0.8525419635058906},
         {-6.654894290848707, -0.4060116340532713, -5.206734506791688, 16.18381989216396, -6.766821758158816, 4.070220552128887, -2.793576478373895, 2.034284048623785, -1.498966108933714, 1.038680284243463},
         {2.812519793890786, 0.3131549663401144, -11.31658321967248, 4.368584145823384, 6.209880791266303, -3.388549054654375, 2.23620098074559, -1.597222531700741, 1.164798340377058, -0.8027842124156417},
         {-0.6970566944109891, -0.09217338512363782, 2.321103493323892, -10.10632939527811, 8.210100252655799, 0.5983454355498914, -0.5813298813502246, 0.4492618922521725, -0.3372000401769037, 0.2352783225581148},
         {-0.2425348845845069, -0.1934948852223001, 0.6313576758207937, -0.8549891349268768, -6.842249723953703, 8.854196834696038, -2.302845959676309, 1.227573086842827, -0.784735155734784, 0.5077221467388234},
         {0.5256415498404903, 0.4453897733073634, -1.410575169768825, 2.332764066019608, -4.372247769077337, -2.411758788747668, 7.309822855203278, -3.106277362044861, 1.865388918196838, -1.178148072928886},
         {-0.4453897733073613, -0.5256415498404904, 1.178148072928884, -1.865388918196839, 3.106277362044853, -7.309822855203274, 2.411758788747677, 4.372247769077325, -2.332764066019597, 1.410575169768821},
         {0.1934948852222991, 0.2425348845845073, -0.5077221467388227, 0.7847351557347858, -1.227573086842825, 2.302845959676301, -8.854196834696047, 6.842249723953711, 0.8549891349268808, -0.6313576758207928},
         {0.09217338512363993, 0.6970566944109982, -0.2352783225581208, 0.3372000401769117, -0.4492618922521805, 0.5813298813502397, -0.5983454355499207, -8.210100252655756, 10.10632939527809, -2.321103493323909},
         {-0.3131549663401145, -2.812519793890797, 0.802784212415643, -1.164798340377063, 1.597222531700743, -2.236200980745592, 3.388549054654387, -6.209880791266316, -4.368584145823379, 11.31658321967248},
         {0.4060116340532701, 6.654894290848723, -1.038680284243459, 1.498966108933717, -2.034284048623785, 2.793576478373891, -4.070220552128887, 6.766821758158809, -16.18381989216392, 5.206734506791658},
         {-0.3336961594984738, 33.95316647156923, 0.8525419635058871, -1.226078925089379, 1.652980387979738, -2.242493058697644, 3.189661480152715, -5.010308567079949, 9.777856994580997, -40.61363058742305}};
    static constexpr S FE_TF1[12][10] =
        {{0.637722700703627, 0.005934313552769153, 0.4643335745781222, -0.1540790292071483, 0.08464294857588464, -0.05530648764425068, 0.03936528942527023, -0.02920837985957846, 0.02174481166755889, -0.01514974179225542},
         {-0.06426597363274085, -0.00323614412574335, 0.9789546446489455, 0.1185565074622587, -0.05241376541297674, 0.03201919547088078, -0.02212183622595462, 0.01616328609716033, -0.01193173104662563, 0.008275816764795482},
         {-0.04110575553289672, -0.005343979960733072, 0.1548185656178299, 0.9665902146731671, -0.1174340149287504, 0.06005750675066619, -0.03882902711074389, 0.02747824762374872, -0.01994383917020894, 0.01371208203792044},
         {0.05609355993900286, 0.01458359681780183, -0.1706629246147103, 0.486220252550122, 0.7671290862485313, -0.2105786159257141, 0.1180876327865069, -0.07877879800280899, 0.05553230186000614, -0.03762609165873794},
         {-0.0380754652335718, -0.01759727700325901, 0.1068573208272038, -0.2064250160507811, 0.7925852915866173, 0.4562245437344592, -0.1733401391702521, 0.1031720774186616, -0.06919520396473032, 0.0457938678556526},
         {0.01163316534779615, 0.009043727570243451, -0.03137880094299823, 0.05276398375949921, -0.104598261814395, 0.9756052004327661, 0.1344818691325706, -0.06116291296890827, 0.03746824786787005, -0.02385621838444384},
         {0.009043727570243493, 0.01163316534779629, -0.02385621838444406, 0.03746824786787033, -0.0611629129689087, 0.134481869132572, 0.9756052004327658, -0.104598261814396, 0.05276398375949955, -0.03137880094299857},
         {-0.01759727700325893, -0.03807546523357176, 0.04579386785565246, -0.06919520396473025, 0.1031720774186609, -0.173340139170251, 0.4562245437344586, 0.7925852915866165, -0.2064250160507798, 0.1068573208272034},
         {0.01458359681780179, 0.05609355993900299, -0.03762609165873798, 0.05553230186000636, -0.0787787980028091, 0.118087632786507, -0.210578615925715, 0.7671290862485339, 0.4862202525501198, -0.1706629246147099},
         {-0.005343979960733051, -0.04110575553289673, 0.01371208203792031, -0.01994383917020896, 0.02747824762374877, -0.03882902711074394, 0.06005750675066622, -0.1174340149287501, 0.966590214673167, 0.1548185656178297},
         {-0.003236144125743398, -0.06426597363274145, 0.008275816764795468, -0.01193173104662575, 0.01616328609716033, -0.02212183622595487, 0.03201919547088092, -0.0524137654129769, 0.1185565074622591, 0.9789546446489457},
         {0.005934313552769167, 0.6377227007036249, -0.01514974179225528, 0.02174481166755911, -0.02920837985957815, 0.03936528942527043, -0.05530648764425071, 0.08464294857588492, -0.1540790292071481, 0.4643335745781249}};
    static constexpr S FE_TF2[12][2] =
        {{0.9907803171233597, 0.009219682876640323},
         {0.9520586281852375, 0.04794137181476255},
         {0.8849513370971525, 0.1150486629028476},
         {0.7936589771433087, 0.2063410228566913},
         {0.6839157494990901, 0.3160842505009099},
         {0.5626167042557345, 0.4373832957442655},
         {0.4373832957442656, 0.5626167042557344},
         {0.3160842505009099, 0.6839157494990901},
         {0.2063410228566915, 0.7936589771433085},
         {0.1150486629028476, 0.8849513370971525},
         {0.04794137181476249, 0.9520586281852375},
         {0.009219682876640378, 0.9907803171233596}};
#elif DEGREE == 10
    // Quadrature rules
    static constexpr S weights[13] = {0.02024200238265796, 0.04606074991886425, 0.06943675510989362, 0.08907299038097291, 0.1039080237684443, 0.1131415901314486, 0.1162757766154369, 0.1131415901314486, 0.1039080237684442, 0.08907299038097297, 0.06943675510989362, 0.04606074991886425, 0.02024200238265796};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[13][11] =
        {{-40.87034750252714, -0.3070406276450841, 48.53479668857029, -11.21382121993334, 5.698244391365392, -3.623061707443906, 2.557025210303191, -1.907638204281561, 1.458948692083808, -1.107369852920982, 0.7802641324293607},
         {-7.010473414106156, 0.3864193569709464, -7.842228738429828, 20.03457283253744, -8.080725201685709, 4.822916204560638, -3.313906010636505, 2.438607743369224, -1.850666537167987, 1.398409914220251, -0.9829261496323227},
         {3.271495883215244, -0.327699198281609, -12.52513914771286, 3.412163123966593, 8.623086962181407, -4.518936259523545, 2.955961554274436, -2.123530888160652, 1.590385054193773, -1.19268193566295, 0.8348948515101586},
         {-1.053537809868541, 0.1569859096812389, 3.358281707010701, -11.84970106662371, 8.334917342188811, 2.043294288195204, -1.421477572989298, 1.024075203746492, -0.7651735621623259, 0.5725394066182783, -0.4002038457968502},
         {-0.02152303869999261, 0.08418452950265037, 0.00522646671295246, 0.3749442824406155, -9.017405475976499, 9.86820971830036, -1.717992880326335, 0.8041061457085501, -0.4917301205836067, 0.3311227710295741, -0.2191423981082716},
         {0.4461541993026195, -0.3308445231471063, -1.185380467342007, 1.899958267138328, -3.225988008827385, -4.84776269846148, 9.054798007066255, -3.320789080970549, 1.951926963986365, -1.303379864415584, 0.8613072056705414},
         {-0.4921874999999947, 0.492187499999996, 1.293879927664193, -2.013010817153546, 3.216255807451021, -6.610353706756648, 7.466201532446678e-15, 6.610353706756642, -3.216255807451022, 2.013010817153552, -1.293879927664196},
         {0.3308445231471062, -0.44615419930262, -0.8613072056705412, 1.303379864415582, -1.951926963986366, 3.320789080970556, -9.054798007066259, 4.847762698461479, 3.225988008827381, -1.899958267138331, 1.185380467342009},
         {-0.08418452950264976, 0.02152303869999073, 0.2191423981082716, -0.3311227710295739, 0.4917301205836033, -0.8041061457085481, 1.717992880326328, -9.868209718300353, 9.017405475976503, -0.3749442824406228, -0.005226466712950018},
         {-0.1569859096812418, 1.053537809868551, 0.4002038457968549, -0.5725394066182865, 0.7651735621623446, -1.024075203746513, 1.421477572989331, -2.04329428819526, -8.334917342188751, 11.84970106662371, -3.358281707010726},
         {0.3276991982816092, -3.271495883215243, -0.8348948515101562, 1.192681935662945, -1.590385054193773, 2.123530888160649, -2.955961554274435, 4.518936259523539, -8.623086962181416, -3.412163123966582, 12.52513914771286},
         {-0.3864193569709449, 7.010473414106151, 0.9829261496323189, -1.398409914220249, 1.85066653716797, -2.438607743369229, 3.313906010636495, -4.822916204560629, 8.080725201685702, -20.03457283253746, 7.842228738429842},
         {0.3070406276450743, 40.87034750252715, -0.7802641324293766, 1.107369852920993, -1.458948692083798, 1.907638204281557, -2.557025210303177, 3.623061707443917, -5.698244391365389, 11.21382121993334, -48.53479668857027}};
    static constexpr S FE_TF1[13][11] =
        {{0.6228692188668682, -0.004965211414775506, 0.482042914837309, -0.158285690637018, 0.08685782361826266, -0.05684518937586132, 0.04067644727067556, -0.03057456482805972, 0.02348551448028499, -0.01787198751152136, 0.01261072469383538},
         {-0.07829693713180771, 0.003364517282770534, 0.9657556299002976, 0.1555080931788738, -0.06762999182957712, 0.0412128265745299, -0.02857111508974977, 0.02112140513283439, -0.01607079153540832, 0.01216184049499146, -0.008555476977754287},
         {-0.02646552595282459, 0.002914855712509861, 0.09736806830567774, 0.9856213545979238, -0.08206862260458277, 0.04123870444728014, -0.02662096111237641, 0.01900908685368607, -0.01419090918399457, 0.01062305525286523, -0.007429106316164329},
         {0.04884872017475304, -0.01063767410406603, -0.1470815046270966, 0.3943784362169553, 0.8368126195707992, -0.2002283198397243, 0.1105198674272538, -0.07397593161356696, 0.0534372293326759, -0.03928605865596593, 0.02721261611798271},
         {-0.03985910527810998, 0.01517617915750686, 0.1111713065938055, -0.2099185476147686, 0.6955026673457365, 0.5717051620128245, -0.1991694659502353, 0.1173249097044252, -0.0800952195488152, 0.05720360465427266, -0.03904149107664204},
         {0.01961231844691406, -0.01226575201285632, -0.05267218151689103, 0.08740392333658234, -0.166515535740482, 0.9180787227969579, 0.2661142417824078, -0.1139272153706984, 0.07005909100936485, -0.04771108343963381, 0.03182347070833409},
         {-7.28583859910259e-17, 2.42861286636753e-17, 4.163336342344337e-17, -6.938893903907228e-17, 1.249000902703301e-16, -2.914335439641036e-16, 1.0, 3.33066907387547e-16, -1.387778780781446e-16, 8.326672684688674e-17, -8.326672684688674e-17},
         {-0.01226575201285636, 0.01961231844691427, 0.03182347070833412, -0.04771108343963374, 0.07005909100936523, -0.1139272153706989, 0.2661142417824088, 0.9180787227969579, -0.1665155357404826, 0.08740392333658295, -0.05267218151689127},
         {0.01517617915750685, -0.03985910527811003, -0.03904149107664204, 0.05720360465427256, -0.08009521954881538, 0.1173249097044255, -0.1991694659502352, 0.5717051620128244, 0.6955026673457365, -0.2099185476147692, 0.1111713065938058},
         {-0.01063767410406596, 0.04884872017475269, 0.02721261611798248, -0.03928605865596554, 0.05343722933267551, -0.07397593161356622, 0.1105198674272529, -0.2002283198397234, 0.8368126195708019, 0.3943784362169512, -0.1470815046270955},
         {0.002914855712509823, -0.0264655259528244, -0.00742910631616437, 0.01062305525286515, -0.01419090918399411, 0.01900908685368603, -0.02662096111237613, 0.04123870444727976, -0.08206862260458235, 0.985621354597924, 0.09736806830567729},
         {0.003364517282770575, -0.07829693713180774, -0.008555476977754287, 0.01216184049499126, -0.01607079153540843, 0.02112140513283417, -0.02857111508974974, 0.04121282657452971, -0.06762999182957714, 0.1555080931788738, 0.9657556299002973},
         {-0.004965211414775544, 0.6228692188668682, 0.01261072469383509, -0.01787198751152141, 0.02348551448028482, -0.03057456482805997, 0.04067644727067544, -0.0568451893758611, 0.08685782361826261, -0.1582856906370183, 0.4820429148373094}};
    static constexpr S FE_TF2[13][2] =
        {{0.9920915273592941, 0.007908472640705932},
         {0.9587991996114891, 0.04120080038851098},
         {0.9007890453666549, 0.09921095463334501},
         {0.8211746697201702, 0.1788253302798298},
         {0.7242463755182234, 0.2757536244817766},
         {0.6152291579775674, 0.3847708420224326},
         {0.5, 0.5},
         {0.3847708420224326, 0.6152291579775674},
         {0.2757536244817765, 0.7242463755182235},
         {0.1788253302798302, 0.8211746697201698},
         {0.09921095463334501, 0.9007890453666549},
         {0.04120080038851098, 0.9587991996114891},
         {0.007908472640705932, 0.9920915273592941}};
#elif DEGREE == 11
    // Quadrature rules
    static constexpr S weights[14] = {0.01755973016587587, 0.04007904357988007, 0.06075928534395163, 0.07860158357909677, 0.09276919873896884, 0.1025992318606478, 0.1076319267315789, 0.1076319267315788, 0.1025992318606478, 0.09276919873896884, 0.07860158357909694, 0.06075928534395163, 0.04007904357988007, 0.01755973016587587};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[14][12] =
        {{-48.39187833294439, 0.28426009592682, 57.09419151010134, -12.69115064550239, 6.395779796483591, -4.058225124923716, 2.868392294893766, -2.152663646350554, 1.669928157182049, -1.310370815046095, 1.01124957759151, -0.7195128674119164},
         {-7.314909974914882, -0.3672183888470929, -10.81755786881188, 24.20614957281839, -9.435737693943253, 5.587177345561948, -3.835549808024327, 2.83556084099725, -2.180811587932809, 1.70236166308221, -1.309614751247428, 0.9301506512618651},
         {3.701800308814341, 0.3331421254296744, -13.64582359675212, 2.179844327842478, 11.29830460084855, -5.713018526891684, 3.703302612129417, -2.662675098513691, 2.016408314305675, -1.559624827639408, 1.19320653080288, -0.8448667703761157},
         {-1.400975456231268, -0.2016463908165976, 4.352046997463638, -13.48327974051922, 8.181833930928084, 3.76188383067378, -2.372009336235666, 1.665516180627439, -1.243203170068166, 0.9531165621676365, -0.7252832077837623, 0.5119997997941081},
         {0.2154005509404647, 0.0008490829669551836, -0.6551894712631643, 1.60221355024539, -11.08850948329033, 10.59371549365313, -0.8240388141476993, 0.2066547973686936, -0.06731494633028801, 0.0218747601628988, -0.005077140358126231, -0.0005783799479188723},
         {0.3270464169829964, 0.225105742698213, -0.8583912356410901, 1.323346575177579, -1.973912433693747, -7.242553223614861, 10.53287060985332, -3.203010494037918, 1.80238361647081, -1.206497657421755, 0.8536362333926419, -0.5800241501661858},
         {-0.484639883246127, -0.4142684543379381, 1.267799825784922, -1.945427021160182, 3.010643905191678, -5.612780908418826, -2.514067314777457, 8.671342778275219, -3.796091940778579, 2.350088795446613, -1.606450963334347, 1.073851181355024},
         {0.4142684543379382, 0.4846398832461263, -1.073851181355024, 1.606450963334357, -2.350088795446622, 3.796091940778583, -8.671342778275218, 2.514067314777442, 5.612780908418835, -3.01064390519167, 1.945427021160177, -1.267799825784927},
         {-0.2251057426982145, -0.327046416982998, 0.5800241501661878, -0.8536362333926526, 1.206497657421766, -1.802383616470818, 3.20301049403793, -10.53287060985332, 7.242553223614845, 1.973912433693763, -1.323346575177587, 0.8583912356410992},
         {-0.0008490829669536293, -0.2154005509404608, 0.00057837994791754, 0.00507714035812934, -0.02187476016290679, 0.06731494633029511, -0.2066547973687043, 0.8240388141477142, -10.59371549365315, 11.08850948329032, -1.60221355024536, 0.655189471263159},
         {0.2016463908166014, 1.400975456231279, -0.5119997997941139, 0.725283207783779, -0.9531165621676556, 1.24320317006819, -1.665516180627469, 2.372009336235704, -3.761883830673848, -8.181833930927958, 13.48327974051919, -4.352046997463694},
         {-0.333142125429677, -3.701800308814333, 0.8448667703761144, -1.19320653080289, 1.559624827639421, -2.016408314305682, 2.662675098513705, -3.703302612129423, 5.713018526891702, -11.29830460084854, -2.179844327842526, 13.64582359675213},
         {0.3672183888470918, 7.314909974914911, -0.9301506512618647, 1.309614751247436, -1.702361663082213, 2.180811587932802, -2.835560840997242, 3.835549808024315, -5.587177345561928, 9.435737693943164, -24.20614957281837, 10.8175578688119},
         {-0.2842600959268164, 48.39187833294448, 0.7195128674119218, -1.011249577591492, 1.310370815046134, -1.66992815718206, 2.15266364635055, -2.868392294893816, 4.058225124923773, -6.395779796483655, 12.69115064550256, -57.09419151010156}};
    static constexpr S FE_TF1[14][12] =
        {{0.6099745919326945, 0.004212151433949363, 0.4973099551884855, -0.1617643749746479, 0.08866033386165635, -0.05806804274570543, 0.04167205929925986, -0.03153622796272114, 0.02458596332078423, -0.01935154121960078, 0.01496229957961064, -0.01065716771376498},
         {-0.08902830207693754, -0.003303881738161936, 0.9519274159942592, 0.1884730884693804, -0.08073200088280609, 0.04905686576164434, -0.03405258696013914, 0.02532012908863877, -0.01953843181544324, 0.01528269969370298, -0.01177131011177988, 0.008366314577642323},
         {-0.01413728294082632, -0.001336964832782062, 0.05105725093128442, 0.9957209627895485, -0.04708010048410489, 0.02326185510868581, -0.01497029873177877, 0.0107280591734249, -0.008109680191216509, 0.006266001876756355, -0.004790882288693329, 0.00339107958970164},
         {0.04127896322736009, 0.007650257193876878, -0.1232621935679718, 0.3157933806318745, 0.8883701109326113, -0.1826754800965878, 0.09919756250273637, -0.06631345367437995, 0.0483675943469781, -0.0366182009251151, 0.0276657408367485, -0.01945428140813128},
         {-0.03906657180575649, -0.01249799768912401, 0.1084248556109784, -0.2011280618242428, 0.6029418791244694, 0.6664563280804426, -0.2116484277683016, 0.1229674091546514, -0.08430532586993567, 0.06177251514416796, -0.04582077288075857, 0.03190417072340921},
         {0.02448863169498032, 0.01264032300488876, -0.06554407881775405, 0.1076442202323588, -0.199020213135027, 0.8443659996174212, 0.3865598655604829, -0.1548727778679921, 0.09448735837574639, -0.06559463034057966, 0.04729486525334127, -0.03244956357786648},
         {-0.007412905724121573, -0.00596712697690914, 0.01943577541581563, -0.03001827218599386, 0.04718187984803903, -0.0923785024507848, 0.9819593715093887, 0.1144029473789543, -0.05304740776484917, 0.03342018029209184, -0.02302268717410441, 0.01544674783247388},
         {-0.005967126976909144, -0.007412905724121583, 0.01544674783247369, -0.02302268717410418, 0.03342018029209184, -0.05304740776484895, 0.1144029473789533, 0.9819593715093889, -0.09237850245078422, 0.04718187984803852, -0.03001827218599344, 0.01943577541581544},
         {0.01264032300488867, 0.02448863169498017, -0.03244956357786646, 0.04729486525334146, -0.06559463034057966, 0.09448735837574623, -0.154872777867992, 0.3865598655604812, 0.8443659996174224, -0.1990202131350259, 0.1076442202323582, -0.06554407881775424},
         {-0.01249799768912409, -0.03906657180575647, 0.0319041707234091, -0.04582077288075891, 0.06177251514416836, -0.08430532586993569, 0.1229674091546512, -0.2116484277683019, 0.6664563280804418, 0.6029418791244695, -0.2011280618242418, 0.1084248556109787},
         {0.007650257193876938, 0.04127896322736013, -0.01945428140813134, 0.0276657408367488, -0.03661820092511556, 0.0483675943469785, -0.06631345367438066, 0.09919756250273665, -0.182675480096589, 0.8883701109326143, 0.3157933806318732, -0.1232621935679724},
         {-0.001336964832781964, -0.01413728294082534, 0.003391079589701411, -0.004790882288692927, 0.006266001876755883, -0.008109680191215635, 0.010728059173424, -0.01497029873177769, 0.02326185510868413, -0.04708010048410084, 0.9957209627895481, 0.05105725093128087},
         {-0.003303881738161932, -0.08902830207693596, 0.008366314577642031, -0.01177131011177976, 0.01528269969370276, -0.01953843181544278, 0.02532012908863841, -0.03405258696013831, 0.04905686576164336, -0.08073200088280383, 0.1884730884693774, 0.9519274159942587},
         {0.004212151433949381, 0.6099745919326969, -0.01065716771376506, 0.01496229957961091, -0.01935154121960059, 0.02458596332078428, -0.03153622796272142, 0.04167205929925979, -0.05806804274570525, 0.08866033386165578, -0.1617643749746487, 0.4973099551884844}};
    static constexpr S FE_TF2[14][2] =
        {{0.9931419043484062, 0.006858095651593843},
         {0.9642174418317868, 0.03578255816821324},
         {0.9136006575348825, 0.08639934246511749},
         {0.8436464524058427, 0.1563535475941573},
         {0.757624318179077, 0.242375681820923},
         {0.6595561844639448, 0.3404438155360551},
         {0.5540274743536718, 0.4459725256463282},
         {0.4459725256463282, 0.5540274743536718},
         {0.3404438155360551, 0.6595561844639449},
         {0.242375681820923, 0.757624318179077},
         {0.1563535475941575, 0.8436464524058425},
         {0.08639934246511749, 0.9136006575348825},
         {0.03578255816821313, 0.9642174418317868},
         {0.006858095651593787, 0.9931419043484062}};
#elif DEGREE == 12
    // Quadrature rules
    static constexpr S weights[15] = {0.01537662099805849, 0.03518302374405406, 0.05357961023358595, 0.06978533896307716, 0.08313460290849697, 0.09308050000778108, 0.0992157426635558, 0.1012891209627806, 0.0992157426635558, 0.09308050000778108, 0.08313460290849697, 0.06978533896307727, 0.05357961023358595, 0.03518302374405393, 0.01537662099805849};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[15][13] =
        {{-56.51631860898831, -0.2645775861457755, 66.28936531124901, -14.2085562218269, 7.10272659090143, -4.495809224743322, 3.178600007125166, -2.39278939341645, 1.870098078027137, -1.491053760228102, 1.192126416324829, -0.9314787674929228, 0.6676671592142185},
         {-7.571766303566354, 0.3489267310897819, -14.12479485591715, 28.69105139657755, -10.82747893427988, 6.361496481507607, -4.35908265054771, 3.228844091383086, -2.500038792218148, 1.981838394297988, -1.578717328703357, 1.230712179246562, -0.8809904088699745},
         {4.106444370865473, -0.3328443290286365, -14.68838007340076, 0.6870913097786069, 14.21966316384789, -6.959549538916682, 4.471628552983729, -3.211159274231061, 2.443692622119178, -1.916994168227776, 1.517064512201818, -1.177827418131551, 0.8411702701397732},
         {-1.737143816631121, 0.2319590660110628, 5.301290871021097, -15.0162826274648, 7.772116855780085, 5.728303266714294, -3.410880285993951, 2.35603018161182, -1.754709608008372, 1.358843890011153, -1.066740502934251, 0.8240886083557819, -0.5868758984727962},
         {0.4594785146444369, -0.06591471013026884, -1.32795665290836, 2.806064750460993, -13.05206342554232, 11.04928297256127, 0.3470730829550628, -0.5309851784456217, 0.4530368952002823, -0.368949206672252, 0.296659229999428, -0.2320831566031554, 0.1663568844805134},
         {0.1817167012139599, -0.1345353724313182, -0.465160971285226, 0.6583130205838998, -0.674743636765651, -9.561585018482234, 11.74343587821872, -2.779080992368082, 1.457377567704301, -0.9528869455309693, 0.6809714321492117, -0.498211812693806, 0.3443901496871967},
         {-0.4362081861266402, 0.3259619260664609, 1.136016946304157, -1.721491069895759, 2.588477375714954, -4.432781284580366, -5.043891709432063, 10.51191850854494, -4.069733476405315, 2.453975271310762, -1.695931226074918, 1.220633843733094, -0.8369469191593061},
         {0.4511718749999988, -0.4511718750000012, -1.165868981244306, 1.72973025901819, -2.48850010973042, 3.876958399506667, -7.895847941765695, -2.835019399347238e-15, 7.895847941765696, -3.876958399506657, 2.488500109730417, -1.729730259018192, 1.165868981244306},
         {-0.3259619260664604, 0.4362081861266434, 0.8369469191593057, -1.220633843733094, 1.695931226074924, -2.45397527131077, 4.069733476405314, -10.51191850854493, 5.043891709432053, 4.43278128458037, -2.588477375714958, 1.721491069895766, -1.136016946304158},
         {0.1345353724313196, -0.1817167012139634, -0.3443901496871994, 0.4982118126938109, -0.68097143214922, 0.9528869455309783, -1.45737756770431, 2.779080992368095, -11.74343587821873, 9.561585018482221, 0.6747436367656765, -0.6583130205839158, 0.4651609712852323},
         {0.06591471013026795, -0.4594785146444418, -0.1663568844805161, 0.2320831566031611, -0.2966592299994306, 0.3689492066722577, -0.4530368952002894, 0.5309851784456296, -0.3470730829550717, -11.04928297256123, 13.05206342554232, -2.806064750461029, 1.327956652908373},
         {-0.2319590660110635, 1.737143816631129, 0.5868758984728035, -0.8240886083557921, 1.066740502934259, -1.358843890011168, 1.754709608008387, -2.356030181611841, 3.410880285993988, -5.728303266714367, -7.772116855780032, 15.01628262746479, -5.301290871021092},
         {0.3328443290286349, -4.106444370865497, -0.841170270139773, 1.177827418131546, -1.517064512201819, 1.916994168227781, -2.443692622119157, 3.211159274231059, -4.47162855298372, 6.959549538916649, -14.21966316384789, -0.687091309778581, 14.68838007340077},
         {-0.3489267310897824, 7.571766303566387, 0.8809904088699723, -1.23071217924656, 1.578717328703366, -1.981838394297996, 2.500038792218164, -3.228844091383092, 4.359082650547721, -6.361496481507605, 10.82747893427993, -28.69105139657756, 14.12479485591706},
         {0.2645775861457809, 56.51631860898826, -0.6676671592141901, 0.9314787674928908, -1.19212641632484, 1.491053760228066, -1.870098078027119, 2.392789393416408, -3.178600007125148, 4.495809224743219, -7.102726590901266, 14.20855622182667, -66.28936531124874}};
    static constexpr S FE_TF1[15][13] =
        {{0.5986824254606026, -0.003616044009224089, 0.5105963662124857, -0.1646768154687279, 0.09014931938246759, -0.05906104435577272, 0.04245543894213302, -0.03225394002182433, 0.02534760440043568, -0.02028055213011537, 0.0162511877797582, -0.01271609341599191, 0.009122147223773733},
         {-0.09735537783320683, 0.003152251306903584, 0.9380915097898307, 0.2179262754802424, -0.09205135217831618, 0.05577292256230335, -0.03872470527620099, 0.02888241045398283, -0.02245319726356629, 0.01784360217542505, -0.01423667980841432, 0.01110949437317989, -0.007957153782163454},
         {-0.003773178855648693, 0.0003098916079206059, 0.01342412340269761, 0.9996815317859097, -0.0133793746827137, 0.006506233795745872, -0.004172035838045832, 0.0029932753284692, -0.002276746264940233, 0.001785498347431497, -0.001412739863484905, 0.00109670542088959, -0.0007831841842308382},
         {0.03394158813488776, -0.005424265645930914, -0.1006686872240984, 0.248919295379035, 0.9260332355472924, -0.1607190031541478, 0.08588681771229753, -0.05723772122243409, 0.04190949692624224, -0.03214944487387349, 0.025096704827452, -0.0193222299128394, 0.0137342135061174},
         {-0.0367524594043664, 0.0100369873410083, 0.1015936544049714, -0.1857940719150531, 0.5179609498556228, 0.7430135782688393, -0.2138439645489539, 0.1224178936838795, -0.08386429040032609, 0.06208619494660283, -0.04747010078193093, 0.0360984991938588, -0.02548287064415207},
         {0.02700117465208697, -0.01173375137827151, -0.07206986049835554, 0.1173829425044841, -0.2120058309524833, 0.7648074989274414, 0.4926532456739119, -0.1839804878726293, 0.1109165054897634, -0.07739841755942062, 0.05730935007058956, -0.04278933428024449, 0.02990696522312807},
         {-0.01294222526584445, 0.008606707301371664, 0.03386158982510016, -0.05199628274241051, 0.08067752771293495, -0.1522840046849972, 0.9378342723407903, 0.2277838717727104, -0.1001221084653868, 0.06271093931546662, -0.0440938877598737, 0.03202314866947548, -0.0220595480193364},
         {5.898059818321144e-17, 6.938893903907228e-18, 6.245004513516506e-17, -3.469446951953614e-17, -5.551115123125783e-17, 2.775557561562891e-17, 0.0, 1.0, -1.387778780781446e-16, -6.938893903907228e-17, 1.110223024625157e-16, 4.163336342344337e-17, -2.081668171172169e-17},
         {0.008606707301371553, -0.01294222526584448, -0.02205954801933632, 0.03202314866947517, -0.04409388775987357, 0.06271093931546654, -0.1001221084653865, 0.2277838717727094, 0.9378342723407904, -0.1522840046849964, 0.08067752771293443, -0.0519962827424103, 0.03386158982509999},
         {-0.01173375137827145, 0.02700117465208707, 0.02990696522312802, -0.04278933428024452, 0.05730935007058951, -0.07739841755942078, 0.110916505489763, -0.1839804878726293, 0.4926532456739107, 0.7648074989274414, -0.2120058309524828, 0.1173829425044843, -0.07206986049835551},
         {0.01003698734100834, -0.03675245940436676, -0.02548287064415235, 0.03609849919385912, -0.04747010078193137, 0.06208619494660339, -0.08386429040032688, 0.1224178936838802, -0.2138439645489554, 0.7430135782688405, 0.5179609498556239, -0.1857940719150552, 0.101593654404972},
         {-0.005424265645930824, 0.03394158813488711, 0.01373421350611705, -0.01932222991283871, 0.02509670482745133, -0.03214944487387274, 0.04190949692624082, -0.05723772122243291, 0.08588681771229542, -0.1607190031541436, 0.9260332355472937, 0.2489192953790286, -0.1006686872240956},
         {0.0003098916079206579, -0.003773178855649442, -0.0007831841842309423, 0.001096705420889632, -0.001412739863485363, 0.001785498347431899, -0.002276746264940566, 0.00299327532846981, -0.00417203583804672, 0.006506233795747191, -0.01337937468271611, 0.9996815317859097, 0.01342412340270008},
         {0.003152251306903529, -0.09735537783320866, -0.007957153782163634, 0.01110949437317994, -0.01423667980841445, 0.01784360217542529, -0.02245319726356551, 0.02888241045398317, -0.0387247052762011, 0.0557729225623036, -0.09205135217831764, 0.2179262754802433, 0.9380915097898326},
         {-0.003616044009224075, 0.5986824254605991, 0.009122147223773809, -0.01271609341599191, 0.01625118777975777, -0.02028055213011554, 0.02534760440043639, -0.03225394002182419, 0.04245543894213288, -0.05906104435577278, 0.09014931938246802, -0.1646768154687287, 0.5105963662124894}};
    static constexpr S FE_TF2[15][2] =
        {{0.9939962590102428, 0.0060037409897572},
         {0.968636696200353, 0.03136330379964702},
         {0.9241032917052137, 0.07589670829478634},
         {0.862208865680085, 0.1377911343199149},
         {0.7854860863042694, 0.2145139136957306},
         {0.6970756735387817, 0.3029243264612183},
         {0.6005970469987172, 0.3994029530012828},
         {0.5, 0.5},
         {0.3994029530012827, 0.6005970469987173},
         {0.3029243264612183, 0.6970756735387817},
         {0.2145139136957306, 0.7854860863042694},
         {0.1377911343199152, 0.8622088656800848},
         {0.07589670829478629, 0.9241032917052137},
         {0.03136330379964702, 0.968636696200353},
         {0.006003740989757256, 0.9939962590102427}};
#elif DEGREE == 13
    // Quadrature rules
    static constexpr S weights[16] = {0.0135762297058771, 0.03112676196932399, 0.04757925584124639, 0.06231448562776695, 0.07479799440828835, 0.08457825969750127, 0.09130170752246179, 0.09472530522753428, 0.09472530522753428, 0.09130170752246179, 0.0845782596975013, 0.07479799440828835, 0.06231448562776706, 0.04757925584124639, 0.03112676196932391, 0.0135762297058771};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[16][14] =
        {{-65.24259042017246, 0.2474086336500605, 76.11848888304911, -15.76507458598777, 7.81889987942111, -4.936135366498856, 3.488679876045751, -2.63033387039971, 2.064382083202275, -1.660284896268507, 1.350707980243037, -1.095289062974004, 0.8640064202326236, -0.6228655535426597},
         {-7.783785951133226, -0.3317618097750386, -17.75788917896018, 33.48365314152839, -12.25274253879006, 7.144756320170913, -4.884670062751259, 3.620229640247117, -2.813093260902834, 2.248384238085237, -1.821748949927357, 1.473315829704748, -1.160215387164717, 0.835567969668265},
         {4.488088472635351, 0.3289410261614163, -15.66094366899875, -1.053830860523352, 17.3750292183484, -8.25027338416516, 5.256324923218533, -3.766537812724022, 2.872135718212155, -2.26929931142976, 1.825165074053807, -1.468979016497268, 1.153242408077327, -0.829062786368695},
         {-2.061297953918778, -0.2520444824457473, 6.207409399046653, -16.45822948887939, 7.123903143964943, 7.922066321297179, -4.52164224070632, 3.083625759335451, -2.289769873997558, 1.780778202497574, -1.418022698268205, 1.133952432795363, -0.886588127743682, 0.6358596070225122},
         {0.7052934873772758, 0.1152240663232728, -1.99989888216647, 3.975812146971755, -14.91110610088411, 11.25400138894261, 1.768108796735846, -1.382952471585469, 1.04618215797028, -0.8160938280350964, 0.6497366788726611, -0.5191621444457968, 0.4056071427670878, -0.2907524388438403},
         {0.0195681157039424, 0.059618078427573, -0.0302894582315909, -0.05836813453598499, 0.6352525932170912, -11.78853094584071, 12.69574763055755, -2.075812341881165, 0.9491579324379826, -0.5794196553231545, 0.4013247804521805, -0.2946644478595815, 0.2186173798973048, -0.1522015270214443},
         {-0.3578041501682075, -0.2413526829237108, 0.9273464413922475, -1.386296530705356, 2.019836301005494, -3.147308718852417, -7.537752693509572, 12.11529221879407, -4.049616156902181, 2.353772906901153, -1.616240753413044, 1.187158071773337, -0.8826725235582562, 0.6156382701664466},
         {0.4502138084456351, 0.3889653815493711, -1.160417566586427, 1.710001747242835, -2.426845584173903, 3.673034443009455, -6.860780416294057, -2.591309826739841, 10.01178309459895, -4.473068920355241, 2.815602044202505, -1.985516026487646, 1.444782118464764, -0.9964442968763974},
         {-0.3889653815493705, -0.4502138084456291, 0.9964442968763929, -1.444782118464754, 1.985516026487629, -2.815602044202478, 4.47306892035519, -10.01178309459892, 2.591309826739855, 6.860780416294043, -3.673034443009444, 2.426845584173896, -1.710001747242826, 1.160417566586418},
         {0.2413526829237104, 0.3578041501682014, -0.6156382701664399, 0.8826725235582457, -1.187158071773321, 1.616240753413019, -2.353772906901115, 4.049616156902137, -12.11529221879409, 7.537752693509624, 3.147308718852403, -2.019836301005484, 1.386296530705346, -0.9273464413922388},
         {-0.05961807842757527, -0.01956811570394462, 0.1522015270214443, -0.2186173798973043, 0.2946644478595783, -0.4013247804521782, 0.5794196553231559, -0.9491579324379917, 2.075812341881184, -12.69574763055756, 11.7885309458407, -0.635252593217086, 0.05836813453597967, 0.0302894582315949},
         {-0.1152240663232704, -0.7052934873772636, 0.2907524388438387, -0.4056071427670842, 0.5191621444457941, -0.6497366788726562, 0.8160938280350867, -1.046182157970264, 1.38295247158547, -1.768108796735841, -11.25400138894265, 14.91110610088412, -3.975812146971733, 1.999898882166458},
         {0.2520444824457468, 2.061297953918758, -0.635859607022514, 0.8865881277436809, -1.133952432795357, 1.418022698268198, -1.780778202497569, 2.289769873997554, -3.083625759335463, 4.521642240706333, -7.922066321297217, -7.123903143964863, 16.45822948887935, -6.207409399046642},
         {-0.3289410261614197, -4.488088472635319, 0.8290627863686963, -1.153242408077322, 1.468979016497254, -1.825165074053794, 2.269299311429735, -2.872135718212145, 3.766537812724033, -5.256324923218549, 8.250273384165183, -17.37502921834836, 1.053830860523208, 15.66094366899879},
         {0.3317618097750465, 7.78378595113338, -0.835567969668261, 1.160215387164719, -1.473315829704756, 1.821748949927362, -2.248384238085233, 2.813093260902846, -3.620229640247152, 4.884670062751335, -7.144756320171, 12.25274253879014, -33.48365314152844, 17.75788917896003},
         {-0.2474086336500871, 65.24259042017222, 0.6228655535426562, -0.8640064202326379, 1.095289062974018, -1.350707980243005, 1.6602848962685, -2.064382083202304, 2.630333870399781, -3.488679876045879, 4.936135366499041, -7.818899879421252, 15.76507458598801, -76.11848888304907}};
    static constexpr S FE_TF1[16][14] =
        {{0.5887161694568323, 0.003136542685683374, 0.5222567423035873, -0.1671419591325821, 0.09139495612950495, -0.05988103461527772, 0.04308777419583505, -0.03281156071227898, 0.02590748738593346, -0.02091685859686544, 0.01706026384009961, -0.0138577452126192, 0.01094358206049575, -0.007894359788348085},
         {-0.1038984354879807, -0.002961350588853953, 0.9245976394101282, 0.2443124806940474, -0.1018736324678446, 0.06155078506216595, -0.04272611822605302, 0.0319196884227837, -0.02491949446314518, 0.01997574508236218, -0.01621643442050337, 0.01313148328789368, -0.01034930835685249, 0.007456952051852699},
         {0.004954925482891399, 0.0003568697706931717, -0.01741150301705895, 0.9994260753944068, 0.0185771266283329, -0.008899807860241984, 0.00568563332833405, -0.004079267464554323, 0.003112736600590763, -0.00246041768426733, 0.001979400322378602, -0.00159338759221081, 0.001251047512011411, -0.0008994314213053509},
         {0.0270965966750142, 0.003775601829268367, -0.07991624990939346, 0.1920680564099745, 0.9531063821528746, -0.1362862227410225, 0.07171622768286268, -0.04761186409913382, 0.03490440365195888, -0.02695120223115599, 0.02136731470874251, -0.01703962076917048, 0.0132994902656293, -0.009528913626449361},
         {-0.03360778587941771, -0.007937772374114128, 0.09259195728712591, -0.1673707985721918, 0.4415310153938323, 0.8042312693484232, -0.2084197174538296, 0.1175537651995886, -0.08030206788784641, 0.05969928676565467, -0.04627577772630397, 0.03638890269630059, -0.02815558599883622, 0.02007330920161454},
         {0.02781843821865784, 0.01034084587868578, -0.07408250842928064, 0.1198418254670398, -0.2124208627711468, 0.6854027146998194, 0.5839617981052572, -0.2025791565576443, 0.1205160681490063, -0.0840803549842319, 0.06293898756143176, -0.04847526030675595, 0.03704075668463098, -0.0262232917154692},
         {-0.01674838750935199, -0.009388224704570699, 0.04374468471885496, -0.06685568338574976, 0.1026555793614485, -0.1883091134745148, 0.8790598564526161, 0.3342263609777889, -0.1388118525850633, 0.08605826193643962, -0.06087393421043388, 0.04543841477321476, -0.03409646569043604, 0.02390050333975724},
         {0.00506007858145514, 0.004181968492884575, -0.01305728005639408, 0.01930119865916761, -0.02756890472141912, 0.04231894251408772, -0.08255085283899724, 0.9861152223209536, 0.09947434898739851, -0.04671484966299998, 0.02986384445757928, -0.02121024252448664, 0.01549183120386086, -0.01070530541309046},
         {0.004181968492884705, 0.005060078581455364, -0.01070530541309082, 0.01549183120386134, -0.02121024252448722, 0.02986384445758027, -0.04671484966300132, 0.09947434898740165, 0.9861152223209537, -0.08255085283900054, 0.04231894251408927, -0.02756890472142002, 0.01930119865916826, -0.01305728005639451},
         {-0.009388224704570576, -0.01674838750935159, 0.02390050333975705, -0.03409646569043568, 0.0454384147732142, -0.06087393421043293, 0.0860582619364379, -0.1388118525850611, 0.3342263609777885, 0.8790598564526148, -0.1883091134745135, 0.1026555793614479, -0.06685568338574915, 0.04374468471885445},
         {0.01034084587868581, 0.02781843821865758, -0.02622329171546908, 0.03704075668463078, -0.04847526030675567, 0.06293898756143138, -0.08408035498423119, 0.1205160681490055, -0.2025791565576453, 0.5839617981052573, 0.6854027146998211, -0.2124208627711474, 0.1198418254670397, -0.07408250842928053},
         {-0.0079377723741141, -0.03360778587941724, 0.02007330920161445, -0.02815558599883597, 0.03638890269630009, -0.04627577772630338, 0.05969928676565403, -0.08030206788784591, 0.1175537651995883, -0.2084197174538292, 0.8042312693484226, 0.4415310153938319, -0.1673707985721908, 0.09259195728712537},
         {0.003775601829268484, 0.02709659667501407, -0.009528913626449361, 0.01329949026562927, -0.01703962076917027, 0.02136731470874229, -0.02695120223115589, 0.03490440365195891, -0.04761186409913386, 0.07171622768286301, -0.1362862227410234, 0.9531063821528759, 0.1920680564099742, -0.07991624990939306},
         {0.000356869770693078, 0.004954925482891181, -0.0008994314213054133, 0.001251047512011327, -0.00159338759221074, 0.001979400322378422, -0.002460417684267621, 0.003112736600590735, -0.004079267464553921, 0.005685633328333412, -0.008899807860241707, 0.0185771266283324, 0.9994260753944064, -0.01741150301705803},
         {-0.002961350588853908, -0.1038984354879779, 0.007456952051852657, -0.01034930835685237, 0.01313148328789322, -0.0162164344205032, 0.01997574508236143, -0.0249194944631445, 0.03191968842278332, -0.04272611822605221, 0.06155078506216448, -0.1018736324678426, 0.2443124806940422, 0.924597639410129},
         {0.003136542685683332, 0.5887161694568334, -0.007894359788347995, 0.01094358206049566, -0.01385774521261941, 0.01706026384009997, -0.02091685859686583, 0.0259074873859339, -0.0328115607122792, 0.04308777419583533, -0.05988103461527826, 0.09139495612950609, -0.167141959132584, 0.5222567423035871}};
    static constexpr S FE_TF2[16][2] =
        {{0.9947004674958251, 0.005299532504174975},
         {0.9722875115366163, 0.0277124884633837},
         {0.9328156011939159, 0.06718439880608412},
         {0.8777022041775016, 0.1222977958224984},
         {0.8089381222013219, 0.1910618777986781},
         {0.7290083888286136, 0.2709916111713863},
         {0.6408017753896295, 0.3591982246103705},
         {0.5475062549188188, 0.4524937450811813},
         {0.4524937450811812, 0.5475062549188188},
         {0.3591982246103705, 0.6408017753896295},
         {0.2709916111713863, 0.7290083888286137},
         {0.1910618777986781, 0.8089381222013219},
         {0.1222977958224986, 0.8777022041775013},
         {0.06718439880608401, 0.9328156011939159},
         {0.0277124884633837, 0.9722875115366163},
         {0.005299532504174975, 0.9947004674958251}};
#elif DEGREE == 14
    // Quadrature rules
    static constexpr S weights[17] = {0.01207415143427396, 0.02772976468699366, 0.04251807415858961, 0.05594192359670195, 0.06756818423426275, 0.07702288053840516, 0.084002051078225, 0.08828135268349632, 0.08972323517810327, 0.08828135268349632, 0.084002051078225, 0.07702288053840516, 0.06756818423426275, 0.05594192359670211, 0.04251807415858961, 0.02772976468699356, 0.01207415143427396};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[17][15] =
        {{-74.56987170560609, -0.2323056863488286, 86.58017070156946, -17.35997543203418, 8.544140134053826, -5.379367207761597, 3.799201783420415, -2.866542660930257, 2.255269605071788, -1.823032943385584, 1.497515686496676, -1.237095927409065, 1.014217015955218, -0.8060660516511042, 0.5837426885593366},
         {-7.953108271656233, 0.3157841336319767, -21.71215925798467, 38.57967099646395, -13.70914946509735, 7.936141501021857, -5.412338253940373, 4.01065728850653, -3.122397118339001, 2.507327113941719, -2.050654706125805, 1.689050870182431, -1.381981521838163, 1.096915811499975, -0.7937591202668095},
         {4.849017879676917, -0.3228077189307487, -16.57024285274164, -3.033158440997129, 20.75503438513338, -9.579044440354995, 6.054107186225423, -4.327086551235796, 3.301444690051894, -2.618469230870635, 2.124392147074117, -1.74038549639723, 1.418833610608631, -1.123508292496529, 0.8118731252543386},
         {-2.373477189874506, 0.2648145414441194, 7.072922450460601, -17.81806020631921, 6.252395275289849, 10.32688834241101, -5.691833034621116, 3.839765536846143, -2.841858942681712, 2.213647752591426, -1.775461132850203, 1.44354569322049, -1.170908716386444, 0.9241573497862026, -0.6665377193166527},
         {0.9495412598590501, -0.1523167913732081, -2.663258222943347, 5.10663240407005, -16.67172868837535, 11.22577194399963, 3.416235454656085, -2.32906696497637, 1.695693352286815, -1.304117362984703, 1.036773154567562, -0.8378605223083038, 0.6768231645045208, -0.5327522438938104, 0.3836300629113827},
         {-0.1528796764825688, 0.001229010431821775, 0.4293664073923515, -0.8023692832693765, 1.934059388117647, -13.91764315663144, 13.40338249972946, -1.118380394493812, 0.3041411138713954, -0.1115612114596987, 0.04447725995532181, -0.01663170308897177, 0.004275435253424886, 0.0009442390249521093, -0.00240992835050946},
         {-0.2578483527426911, 0.1657113362626074, 0.6637998432932503, -0.9733699638180622, 1.354036682718476, -1.807366537252254, -9.96601298185724, 13.47939076065436, -3.754292598022678, 2.076285907597545, -1.400567778194427, 1.029707574713964, -0.7849613849648022, 0.5962419042865063, -0.4207544126745537},
         {0.4194317199890948, -0.3184464019696864, -1.078563440995271, 1.579606630187228, -2.214247176183476, 3.265128421381387, -5.65818613752888, -5.193345553215223, 11.92827037176191, -4.796590843036384, 2.941386540692482, -2.071351315215259, 1.543149796076603, -1.15701789343893, 0.8107752814944074},
         {-0.4189453125000034, 0.4189453125000054, 1.071371014580488, -1.546477685099753, 2.107599018099593, -2.943450903493874, 4.529627283809281, -9.177829929249931, -2.069520888089573e-15, 9.177829929249938, -4.529627283809281, 2.943450903493874, -2.107599018099591, 1.546477685099742, -1.071371014580486},
         {0.3184464019696851, -0.4194317199890969, -0.8107752814944079, 1.157017893438936, -1.543149796076601, 2.071351315215256, -2.941386540692483, 4.796590843036375, -11.92827037176191, 5.193345553215231, 5.658186137528873, -3.265128421381378, 2.214247176183464, -1.579606630187215, 1.078563440995268},
         {-0.1657113362626081, 0.2578483527426956, 0.4207544126745575, -0.5962419042865141, 0.7849613849648067, -1.029707574713969, 1.400567778194432, -2.076285907597549, 3.754292598022689, -13.47939076065438, 9.966012981857226, 1.807366537252276, -1.354036682718482, 0.9733699638180662, -0.6637998432932563},
         {-0.001229010431820887, 0.1528796764825698, 0.002409928350509905, -0.0009442390249487786, -0.004275435253423776, 0.01663170308897488, -0.04447725995532803, 0.1115612114597018, -0.3041411138714007, 1.118380394493825, -13.40338249972946, 13.91764315663143, -1.934059388117644, 0.8023692832693681, -0.4293664073923484},
         {0.1523167913732071, -0.9495412598590536, -0.3836300629113836, 0.532752243893817, -0.6768231645045248, 0.8378605223083029, -1.036773154567559, 1.304117362984702, -1.695693352286819, 2.329066964976378, -3.416235454656092, -11.22577194399966, 16.67172868837534, -5.106632404070006, 2.663258222943339},
         {-0.2648145414441185, 2.373477189874522, 0.6665377193166574, -0.9241573497862237, 1.170908716386452, -1.443545693220499, 1.775461132850215, -2.21364775259144, 2.841858942681732, -3.839765536846185, 5.691833034621163, -10.3268883424111, -6.252395275289674, 17.81806020631909, -7.072922450460602},
         {0.3228077189307508, -4.849017879676952, -0.8118731252543413, 1.123508292496544, -1.418833610608636, 1.74038549639724, -2.124392147074138, 2.618469230870644, -3.301444690051909, 4.327086551235836, -6.054107186225446, 9.579044440355059, -20.75503438513343, 3.033158440997209, 16.5702428527416},
         {-0.31578413363198, 7.953108271656181, 0.7937591202668046, -1.096915811499966, 1.381981521838163, -1.689050870182423, 2.050654706125816, -2.507327113941717, 3.122397118339005, -4.010657288506509, 5.412338253940341, -7.936141501021842, 13.70914946509727, -38.57967099646382, 21.71215925798469},
         {0.2323056863488162, 74.5698717056063, -0.5837426885593615, 0.8060660516510758, -1.014217015955186, 1.237095927409079, -1.497515686496641, 1.823032943385549, -2.255269605071753, 2.866542660930143, -3.799201783420337, 5.379367207761497, -8.544140134053695, 17.35997543203385, -86.58017070156946}};
    static constexpr S FE_TF1[17][15] =
        {{0.5798584174116606, -0.002745381944473001, 0.5325673072832593, -0.1692489008779454, 0.09244847325490366, -0.06056746263458869, 0.0436081780688085, -0.03325761680298735, 0.02633701255831444, -0.02137937726044856, 0.01761174478332808, -0.01457728371149952, 0.01196683117091846, -0.00951914039534961, 0.006897199096097956},
         {-0.1090961437668428, 0.002758588144613459, 0.9116346388809646, 0.2680285739812601, -0.1104382277237163, 0.06654786471949757, -0.04617203567014827, 0.03452471632941317, -0.0270216738673759, 0.02177197831521485, -0.01784632578214435, 0.01472165636771863, -0.01205764561276292, 0.009576932474507971, -0.006932896790200767},
         {0.01233118531955398, -0.0007854283783485239, -0.0428831507329332, 0.996287964839225, 0.04859899523327148, -0.02295857345326716, 0.01461343509683424, -0.01047872731395368, 0.008009246470861639, -0.006359250191445509, 0.005162939889550175, -0.004231657485022952, 0.0034509007726123, -0.00273316466249468, 0.001975284595556834},
         {0.02084403186354419, -0.002556325156328709, -0.06118447661416334, 0.1436644875486481, 0.9721223642043628, -0.1106812504629042, 0.0573932856535023, -0.03794911370714309, 0.02781898177646176, -0.02155277623142484, 0.01722898221451805, -0.01397790940447111, 0.01132187283945984, -0.008927807945221015, 0.006435653421158545},
         {-0.03007189813450972, 0.0062102051853306, 0.08261790184313714, -0.1478982159075032, 0.3735686219307343, 0.8527918678796005, -0.1975190641004753, 0.1098009309856624, -0.07472539061919622, 0.05563051962985166, -0.04341848920275455, 0.03469259240776938, -0.02782292735929738, 0.0218010658447426, -0.01565772038309241},
         {0.0274706872790269, -0.008849613014431613, -0.07301686176671388, 0.1174476832953934, -0.2049968802715718, 0.6094320347929062, 0.6613822847329978, -0.2124034411603437, 0.1246498408776687, -0.0867465021652718, 0.0652309182874638, -0.05095885092780711, 0.0402905357825601, -0.03129015433332404, 0.02235831859144706},
         {-0.01913026638168025, 0.009185033626930465, 0.04989450771142878, -0.07595580541621172, 0.1156274299152616, -0.2072943225756693, 0.8132525816653133, 0.4309209283175289, -0.1686900520388859, 0.1033428645573957, -0.07312741633666481, 0.05521182972135208, -0.04276450430309477, 0.032798329838023, -0.0232711383010272},
         {0.009049065804643893, -0.00630806150640398, -0.02332367473965373, 0.03437298192392864, -0.04880770014401384, 0.0740225176894746, -0.1396181726376212, 0.9512282193140776, 0.1988348952635831, -0.08905071921823127, 0.05643102397550051, -0.04035049262347534, 0.03031222963705841, -0.0228362606631201, 0.01604414792425294},
         {-3.816391647148976e-17, 2.775557561562891e-17, 0.0, 6.938893903907228e-17, -1.040834085586084e-16, 4.85722573273506e-17, 2.498001805406602e-16, -5.551115123125783e-17, 1.0, 1.387778780781446e-17, -1.387778780781446e-16, 4.163336342344337e-17, 1.595945597898663e-16, 0.0, -1.387778780781446e-17},
         {-0.006308061506403971, 0.00904906580464394, 0.01604414792425296, -0.02283626066312026, 0.03031222963705846, -0.04035049262347547, 0.05643102397550074, -0.08905071921823145, 0.1988348952635837, 0.9512282193140779, -0.1396181726376216, 0.07402251768947485, -0.04880770014401394, 0.03437298192392835, -0.02332367473965366},
         {0.009185033626930359, -0.01913026638168032, -0.02327113830102718, 0.03279832983802309, -0.04276450430309466, 0.05521182972135186, -0.07312741633666457, 0.1033428645573953, -0.1686900520388853, 0.4309209283175276, 0.813252581665314, -0.2072943225756685, 0.1156274299152609, -0.07595580541621101, 0.04989450771142855},
         {-0.008849613014431606, 0.02747068727902717, 0.02235831859144718, -0.03129015433332437, 0.04029053578256039, -0.05095885092780726, 0.06523091828746376, -0.08674650216527184, 0.1246498408776691, -0.2124034411603445, 0.661382284732998, 0.6094320347929072, -0.2049968802715721, 0.117447683295393, -0.073016861766714},
         {0.006210205185330567, -0.03007189813450965, -0.0156577203830923, 0.02180106584474289, -0.02782292735929718, 0.03469259240776924, -0.04341848920275453, 0.05563051962985131, -0.07472539061919584, 0.1098009309856623, -0.1975190641004741, 0.8527918678796, 0.3735686219307328, -0.1478982159075017, 0.08261790184313661},
         {-0.002556325156328686, 0.02084403186354426, 0.006435653421158517, -0.008927807945220981, 0.0113218728394598, -0.01397790940447114, 0.01722898221451819, -0.02155277623142499, 0.0278189817764617, -0.03794911370714308, 0.05739328565350224, -0.1106812504629045, 0.9721223642043643, 0.1436644875486467, -0.06118447661416303},
         {-0.0007854283783484371, 0.01233118531955458, 0.001975284595557049, -0.002733164662495013, 0.003450900772612495, -0.004231657485023091, 0.005162939889550119, -0.006359250191445884, 0.008009246470861903, -0.01047872731395434, 0.01461343509683524, -0.02295857345326845, 0.0485989952332742, 0.9962879648392244, -0.04288315073293505},
         {0.002758588144613563, -0.1090961437668451, -0.006932896790200746, 0.00957693247450829, -0.01205764561276307, 0.01472165636771898, -0.01784632578214494, 0.02177197831521505, -0.02702167386737621, 0.03452471632941423, -0.04617203567014881, 0.06654786471949868, -0.1104382277237173, 0.2680285739812629, 0.9116346388809644},
         {-0.002745381944473055, 0.5798584174116592, 0.006897199096097754, -0.009519140395349429, 0.01196683117091862, -0.01457728371149915, 0.01761174478332804, -0.02137937726044854, 0.02633701255831442, -0.03325761680298694, 0.04360817806880818, -0.06056746263458837, 0.0924484732549023, -0.1692489008779431, 0.5325673072832595}};
    static constexpr S FE_TF2[17][2] =
        {{0.9952877376572087, 0.004712262342791262},
         {0.9753377608843838, 0.0246622391156161},
         {0.940119576863493, 0.05988042313650699},
         {0.8907570019484007, 0.1092429980515993},
         {0.8288355796083454, 0.1711644203916546},
         {0.7563452685432385, 0.2436547314567615},
         {0.6756158817269382, 0.3243841182730618},
         {0.5892420907479239, 0.4107579092520761},
         {0.5, 0.5},
         {0.4107579092520761, 0.5892420907479239},
         {0.3243841182730618, 0.6756158817269382},
         {0.2436547314567615, 0.7563452685432385},
         {0.1711644203916547, 0.8288355796083453},
         {0.1092429980515995, 0.8907570019484006},
         {0.05988042313650704, 0.9401195768634929},
         {0.0246622391156161, 0.9753377608843838},
         {0.004712262342791262, 0.9952877376572087}};
#elif DEGREE == 15
    // Quadrature rules
    static constexpr S weights[18] = {0.01080800676324169, 0.02485727444748486, 0.03821286512744452, 0.05047102205314358, 0.06127760335573925, 0.07032145733532535, 0.07734233756313257, 0.08213824187291641, 0.08457119148157179, 0.08457119148157179, 0.08213824187291641, 0.07734233756313257, 0.07032145733532535, 0.06127760335573925, 0.0504710220531437, 0.03821286512744452, 0.02485727444748486, 0.01080800676324169};
    // Precomputed values of basis functions and precomputations
    // FE* dimensions: [permutation][entities][points][dofs]
    static constexpr S FE_TF0[18][16] =
        {{-84.49752500376539, 0.2189205636851393, 97.67333455808924, -18.99269730917851, 9.278316230645647, -5.825590403219591, 4.110498862284734, -3.102132625981371, 2.444135668779152, -1.981872298378487, 1.637474626283375, -1.367011302097875, 1.142877470950593, -0.9451499041355507, 0.7556916569998139, -0.5492707909608254},
         {-8.08142145304345, -0.3009734036903823, -25.98393598246374, 43.97579665306263, -15.19491127758103, 8.735060363804623, -5.94207279711021, 4.400670022781415, -3.429310127980366, 2.761467171834247, -2.271077002448039, 1.889953397386763, -1.576589403855399, 1.301833201332839, -1.039818122454977, 0.7553287604250896},
         {5.191170252422818, 0.3153404307039362, -17.42182898646096, -5.243073710147115, 24.35237259651862, -10.94123168717278, 6.862609067000532, -4.891585863681448, 3.73133603223116, -2.965455559731448, 2.417965049742915, -2.00047544867473, 1.662062760222864, -1.368602269794658, 1.091136902956511, -0.7917395661362209},
         {-2.674119448418854, -0.2723327597098668, 7.900752112933331, -19.1038023989952, 5.170238564145539, 12.92977849955001, -6.9119049379811, 4.618226330063626, -3.406470680213315, 2.654049519255646, -2.13679806543817, 1.752942926581103, -1.447992252437515, 1.187619525461955, -0.9443777272959875, 0.6841907924987914},
         {1.190255424608693, 0.1800023664830708, -3.313636979814022, 6.19708125364627, -18.34119471138829, 10.98064596540873, 5.272484808190494, -3.353435530823079, 2.389438443812519, -1.8225961810224, 1.447552410861825, -1.176770381549673, 0.9660445371252364, -0.7889913761250205, 0.6256505059796955, -0.4525305553940433},
         {-0.33113796793065, -0.05013204427962242, 0.9023254597626094, -1.557667596325716, 3.2086987616646, -15.94916959884358, 13.88126447006364, 0.07062863518686099, -0.4559279415565338, 0.4316986163675574, -0.370597386708296, 0.3128012401148045, -0.2622053222952942, 0.2168285651247098, -0.1732317361026112, 0.1258238457575236},
         {-0.1426796004201693, -0.1003938841382555, 0.3618228843111805, -0.5072568660318488, 0.6256837541481485, -0.4465453918702073, -12.31295546742206, 14.61009993171836, -3.204156298185301, 1.645595187970087, -1.073584645906425, 0.7785765844053193, -0.5950417481388127, 0.4640692559593959, -0.3573632372223781, 0.2541295408229702},
         {0.3656390959431681, 0.2491521360159747, -0.9380068660022135, 1.365117989957902, -1.889463017063691, 2.712553389946635, -4.352076843884221, -7.764057694138114, 13.63146121488286, -4.856215944984164, 2.88473756231785, -2.013189652244092, 1.50915210285227, -1.164353201306531, 0.8910703272721276, -0.631520599565764},
         {-0.4213249112886785, -0.367821498693786, 1.075869007037176, -1.547118820010539, 2.093712982550063, -2.886957410766961, 4.327759729629552, -8.114147300578777, -2.651670633098477, 11.33786615603536, -5.141456887678253, 3.2704907705395, -2.343888788002134, 1.76365140345654, -1.33020802235479, 0.935244223224205},
         {0.367821498693792, 0.4213249112886705, -0.935244223224209, 1.330208022354791, -1.763651403456526, 2.343888788002112, -3.270490770539479, 5.141456887678236, -11.33786615603536, 2.651670633098497, 8.114147300578772, -4.327759729629566, 2.886957410766974, -2.09371298255007, 1.547118820010541, -1.075869007037169},
         {-0.2491521360159777, -0.3656390959431591, 0.6315205995657645, -0.8910703272721265, 1.164353201306516, -1.509152102852247, 2.013189652244069, -2.884737562317826, 4.856215944984141, -13.63146121488285, 7.764057694138129, 4.352076843884201, -2.712553389946631, 1.889463017063689, -1.365117989957897, 0.9380068660022016},
         {0.1003938841382564, 0.1426796004201656, -0.2541295408229697, 0.3573632372223758, -0.4640692559593899, 0.5950417481388053, -0.7785765844053111, 1.073584645906412, -1.645595187970073, 3.204156298185283, -14.61009993171838, 12.31295546742209, 0.4465453918702129, -0.6256837541481532, 0.5072568660318456, -0.3618228843111754},
         {0.05013204427962338, 0.3311379679306435, -0.1258238457575232, 0.1732317361026089, -0.2168285651247077, 0.2622053222952831, -0.3128012401147975, 0.3705973867083023, -0.4316986163675514, 0.4559279415565234, -0.07062863518685196, -13.88126447006364, 15.94916959884359, -3.208698761664618, 1.557667596325711, -0.902325459762603},
         {-0.180002366483073, -1.190255424608657, 0.4525305553940385, -0.6256505059796839, 0.788991376125012, -0.966044537125211, 1.176770381549651, -1.447552410861819, 1.822596181022372, -2.38943844381249, 3.353435530823065, -5.27248480819046, -10.98064596540885, 18.34119471138834, -6.197081253646203, 3.313636979813967},
         {0.2723327597098713, 2.674119448418787, -0.684190792498788, 0.9443777272959836, -1.187619525461943, 1.447992252437494, -1.752942926581085, 2.136798065438156, -2.654049519255635, 3.406470680213299, -4.618226330063625, 6.91190493798109, -12.92977849955006, -5.170238564145402, 19.10380239899508, -7.900752112933227},
         {-0.3153404307039431, -5.191170252422793, 0.7917395661362305, -1.091136902956524, 1.368602269794647, -1.662062760222863, 2.000475448674729, -2.417965049742913, 2.965455559731468, -3.731336032231178, 4.891585863681494, -6.862609067000586, 10.9412316871729, -24.35237259651858, 5.243073710146824, 17.4218289864611},
         {0.3009734036903866, 8.081421453043486, -0.7553287604251109, 1.039818122454994, -1.301833201332851, 1.576589403855432, -1.889953397386782, 2.271077002448042, -2.761467171834274, 3.429310127980413, -4.400670022781492, 5.942072797110296, -8.735060363804791, 15.19491127758124, -43.97579665306272, 25.98393598246374},
         {-0.2189205636851433, 84.49752500376522, 0.5492707909608997, -0.7556916569998954, 0.945149904135595, -1.142877470950726, 1.367011302097979, -1.637474626283302, 1.981872298378651, -2.444135668779433, 3.102132625981572, -4.110498862285061, 5.825590403220083, -9.27831623064627, 18.99269730917956, -97.67333455808971}};
    static constexpr S FE_TF1[18][16] =
        {{0.5719361722450154, 0.002422308525655256, 0.5417460539150715, -0.1710654613414243, 0.0933481057827369, -0.06114870920578144, 0.04404309107539373, -0.0336223966531938, 0.02667711739329708, -0.02173043597548865, 0.01800980100434651, -0.01506745067709045, 0.01261604928379824, -0.01044434492749175, 0.0083566139058206, -0.006076514350665072},
         {-0.113264556602765, -0.002558298135858997, 0.8992953658759282, 0.2894202521406157, -0.1179434241447955, 0.07089309094459281, -0.0491564603952166, 0.03677271029072582, -0.02882647332022455, 0.02330052429895587, -0.01921122393425534, 0.0160151381570001, -0.01337600628031888, 0.01105425288582433, -0.008834365059463792, 0.006419473279256347},
         {0.01859254613980583, 0.001055051777799288, -0.06408763915193927, 0.9911983119727897, 0.07664936678246356, -0.03573751181735645, 0.02266765263104013, -0.01624013911674187, 0.01242300216724307, -0.009890071052810921, 0.00807315817913494, -0.006684290464859136, 0.005556434487149715, -0.004577009703620468, 0.003649950621608204, -0.002648813451706288},
         {0.01520001904181312, 0.001654199788379335, -0.04443665489166863, 0.1023326661421899, 0.9850075467936659, -0.08477287209043582, 0.04335015063971009, -0.02854797698033722, 0.02091238767643242, -0.01622973311685689, 0.01303501181767197, -0.01067635460168711, 0.00880953880836342, -0.007220160245179736, 0.005738613651847287, -0.004156382433908096},
         {-0.02641717271627805, -0.004814576471193247, 0.07240305598373604, -0.1285511157972771, 0.3134970151045354, 0.891031235515196, -0.1828001590858791, 0.1002063455667021, -0.06791756497345117, 0.05056140304894842, -0.0396024425278735, 0.03191583501343046, -0.026052139329255, 0.02119735184883346, -0.01676805519547617, 0.01211098401530129},
         {0.02635314397161464, 0.007437902054639902, -0.06993263461638302, 0.1119452008099547, -0.1928977681354388, 0.5385245126526723, 0.7263675668478893, -0.215170029368396, 0.1245793077110377, -0.0863907967629534, 0.06506664818268788, -0.05119528411839215, 0.04114721938937242, -0.0331400944897936, 0.02604405205384045, -0.01873894618235262},
         {-0.0203998301555398, -0.008500206497380812, 0.05314147655237219, -0.08063102638535041, 0.1218592472755317, -0.2143943186180332, 0.7452988379770575, 0.5169178036633419, -0.1903060442193441, 0.1151399933482081, -0.08130310861750042, 0.06168546146969968, -0.04848171133812787, 0.03849427520888209, -0.0299809240177792, 0.02146007435396263},
         {0.01200886981000515, 0.007176371729611137, -0.03092256318279928, 0.04545677852951981, -0.06422825086982678, 0.09643920139249527, -0.1770567599909346, 0.903394764244587, 0.2937526576979588, -0.1252146687650399, 0.07850510726872532, -0.05618563851032349, 0.04271282752129792, -0.03323237643552818, 0.02556264580926974, -0.01816896624901806},
         {-0.003633042940042022, -0.003065199361566499, 0.009283128204513104, -0.01337183775934564, 0.01815448522613197, -0.02518457129191434, 0.03824071657184539, -0.07454047979562203, 0.9889829645511662, 0.08795702038942722, -0.04167994897377672, 0.02688240653398868, -0.0193920388771963, 0.01464381910554336, -0.0110677200911923, 0.007790298508039748},
         {-0.003065199361566571, -0.003633042940042032, 0.007790298508039862, -0.0110677200911922, 0.01464381910554357, -0.01939203887719652, 0.02688240653398893, -0.04167994897377725, 0.08795702038942843, 0.9889829645511661, -0.07454047979562317, 0.0382407165718461, -0.02518457129191488, 0.01815448522613241, -0.0133718377593456, 0.009283128204513167},
         {0.007176371729611314, 0.01200886981000497, -0.01816896624901821, 0.02556264580926991, -0.03323237643552807, 0.04271282752129766, -0.05618563851032328, 0.07850510726872545, -0.1252146687650403, 0.2937526576979604, 0.9033947642445868, -0.1770567599909362, 0.09643920139249616, -0.06422825086982722, 0.04545677852952011, -0.03092256318279926},
         {-0.008500206497380884, -0.02039983015553929, 0.0214600743539626, -0.02998092401777903, 0.03849427520888146, -0.04848171133812688, 0.06168546146969882, -0.08130310861749959, 0.1151399933482071, -0.1903060442193425, 0.5169178036633413, 0.7452988379770568, -0.2143943186180326, 0.1218592472755315, -0.08063102638534977, 0.05314147655237146},
         {0.007437902054640054, 0.02635314397161414, -0.01873894618235264, 0.02604405205384034, -0.03314009448979347, 0.04114721938937194, -0.05119528411839182, 0.06506664818268795, -0.08639079676295301, 0.1245793077110377, -0.2151700293683972, 0.7263675668478899, 0.5385245126526744, -0.1928977681354404, 0.1119452008099546, -0.06993263461638261},
         {-0.004814576471193246, -0.02641717271627714, 0.01211098401530116, -0.01676805519547593, 0.02119735184883303, -0.02605213932925455, 0.03191583501342985, -0.03960244252787249, 0.0505614030489477, -0.0679175649734502, 0.100206345566701, -0.1828001590858769, 0.8910312355151937, 0.313497015104535, -0.1285511157972754, 0.07240305598373441},
         {0.001654199788379379, 0.01520001904181307, -0.004156382433908197, 0.005738613651847587, -0.007220160245179471, 0.008809538808363595, -0.01067635460168726, 0.0130350118176718, -0.01622973311685729, 0.02091238767643283, -0.02854797698033758, 0.04335015063971061, -0.08477287209043773, 0.9850075467936662, 0.1023326661421911, -0.04443665489166856},
         {0.001055051777799294, 0.0185925461398043, -0.002648813451705974, 0.003649950621607948, -0.004577009703620166, 0.005556434487149142, -0.006684290464858707, 0.008073158179134577, -0.009890071052810347, 0.01242300216724225, -0.01624013911674084, 0.02266765263103865, -0.0357375118173545, 0.07664936678245904, 0.991198311972791, -0.06408763915193562},
         {-0.002558298135859074, -0.1132645566027624, 0.006419473279256334, -0.008834365059463865, 0.01105425288582396, -0.01337600628031847, 0.01601513815699981, -0.01921122393425558, 0.02330052429895596, -0.02882647332022408, 0.03677271029072577, -0.04915646039521661, 0.07089309094459324, -0.1179434241447954, 0.2894202521406129, 0.8992953658759275},
         {0.002422308525655303, 0.5719361722450236, -0.006076514350665131, 0.008356613905820549, -0.01044434492749203, 0.01261604928379852, -0.01506745067709063, 0.01800980100434715, -0.02173043597548877, 0.02667711739329692, -0.03362239665319432, 0.04404309107539418, -0.06114870920578232, 0.09334810578273844, -0.1710654613414243, 0.5417460539150631}};
    static constexpr S FE_TF2[18][2] =
        {{0.9957825842104655, 0.004217415789534495},
         {0.9779119747856988, 0.02208802521430114},
         {0.9463012332487779, 0.05369876675122209},
         {0.9018524794862617, 0.09814752051373837},
         {0.8458435215301766, 0.1541564784698234},
         {0.7798854155369738, 0.2201145844630262},
         {0.7058755807314213, 0.2941244192685787},
         {0.6259431128457528, 0.3740568871542472},
         {0.5423875065208676, 0.4576124934791324},
         {0.4576124934791324, 0.5423875065208676},
         {0.3740568871542472, 0.6259431128457528},
         {0.2941244192685787, 0.7058755807314213},
         {0.2201145844630262, 0.7798854155369738},
         {0.1541564784698234, 0.8458435215301766},
         {0.09814752051373865, 0.9018524794862613},
         {0.05369876675122204, 0.9463012332487779},
         {0.0220880252143012, 0.9779119747856988},
         {0.00421741578953444, 0.9957825842104655}};
#endif
};
