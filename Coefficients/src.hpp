#include "blis.h"
#include <iostream>

//------------------------------------------------------------
template <typename T, int nq, int np>
inline void compute_coefficients_sf(T __restrict__ *phi, __restrict__ T *u,
                                    T __restrict__ *w) {
  constexpr int qoffset = nq * nq;
  constexpr int doffset = nq * nq;

  T T0[nq][np][np] = {0};
  T T1[nq][nq][np] = {0};

  // First tensor contraction
  for (int q2 = 0; q2 < nq; q2++) {
    for (int i0 = 0; i0 < np; i0++) {
      for (int i1 = 0; i1 < np; i1++) {
        for (int i2 = 0; i2 < np; i2++) {
          T0[q2][i0][i1] += phi[q2 * np + i2] * u[i0 * doffset + i1 * np + i2];
        }
      }
    }
  }

  // Second tensor contraction
  for (int q1 = 0; q1 < nq; q1++) {
    for (int q2 = 0; q2 < nq; q2++) {
      for (int i0 = 0; i0 < np; i0++) {
        for (int i1 = 0; i1 < np; i1++) {
          T1[q1][q2][i0] += phi[q1 * np + i1] * T0[q2][i0][i1];
        }
      }
    }
  }

  // Third tensor contraction
  for (int q0 = 0; q0 < nq; q0++) {
    for (int q1 = 0; q1 < nq; q1++) {
      for (int q2 = 0; q2 < nq; q2++) {
        for (int i0 = 0; i0 < np; i0++) {
          w[q0 * qoffset + q1 * nq + q2] += phi[q0 * np + i0] * T1[q1][q2][i0];
        }
      }
    }
  }
}

//------------------------------------------------------------
template <typename T, int nq, int np>
inline void compute_coefficients_sfgemm(T __restrict__ *phi, __restrict__ T *u,
                                        T __restrict__ *w) {

  // First tensor contraction as a GEMM
  // [q0, i1, i2] <- [q0, i0] [i0, i1, i2]
  T T0[nq * np * np] = {0};
  for (int q0 = 0; q0 < nq; q0++) {
    for (int i0 = 0; i0 < np; i0++) {
      for (int d = 0; d < np * np; d++) {
        T0[q0 * np * np + d] = phi[q0 * np + i0] * u[i0 * np * np + d];
      }
    }
  }

  // Tensor transposition
  // [i1, q0, i2] <- [q0, i1, i2]
  T W0[nq * np * np] = {0};
  for (int i1 = 0; i1 < np; i1++) {
    for (int q0 = 0; q0 < nq; q0++) {
      for (int i2 = 0; i2 < np; i2++) {
        W0[i1 * np * nq + q0 * np + i2] = T0[q0 * np * np + i1 * np + i2];
      }
    }
  }

  // Second tensor contraction as a GEMM
  // [q1, q0, i2] <- [q1, i1][i1, q0, i2]
  T T1[nq * nq * np] = {0};
  for (int q1 = 0; q1 < nq; q1++) {
    for (int i1 = 0; i1 < np; i1++) {
      for (int d = 0; d < nq * np; d++) {
        T1[q1 * nq * np + d] = phi[q1 * np + i1] * W0[i1 * nq * np + d];
      }
    }
  }

  // Tensor transposition
  // [i2, q1, q0] <- [q1, q0, i2]
  T W1[nq * nq * np] = {0};
  for (int i2 = 0; i2 < np; i2++) {
    for (int q1 = 0; q1 < nq; q1++) {
      for (int q0 = 0; q0 < nq; q0++) {
        W1[i2 * nq * nq + q1 * nq + q0] = T1[q1 * nq * np + q0 * np + i2];
      }
    }
  }

  // Third tensor contraction as a GEMM
  // [q2, q1, q0] <-[q2, i2][i2, q1, q0]
  T T2[nq * nq * nq] = {0};
  for (int q2 = 0; q2 < nq; q2++) {
    for (int i2 = 0; i2 < np; i2++) {
      for (int d = 0; d < nq * nq; d++) {
        T2[q2 * nq * nq + d] = phi[q2 * np + i2] * W1[i2 * nq * nq + d];
      }
    }
  }

  // Tensor transposition
  // [q0, q1, qw] <- [q2, q1, q0]
  for (int q0 = 0; q0 < np; q0++) {
    for (int q1 = 0; q1 < nq; q1++) {
      for (int q2 = 0; q2 < nq; q2++) {
        w[q0 * nq * nq + q1 * nq + q2] = T2[q2 * nq * nq + q1 * nq + q0];
      }
    }
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_auto(T __restrict__ *phi, __restrict__ T *u,
                                      T __restrict__ *w) {
  T acc = 0;

  for (int q = 0; q < nq; q++) {
    for (int i = 0; i < nd; i++)
#pragma omp simd
      for (int j = 0; j < 4; j++)
        w[q * 4 + j] += phi[q * nd + i] * u[i * 4 + j];
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_auto(T __restrict__ *dphi0,
                                      T __restrict__ *dphi1,
                                      T __restrict__ *dphi2, T __restrict__ *u,
                                      T __restrict__ *w) {
  for (int q = 0; q < nq; q++) {
    T _w0 = 0;
    T _w1 = 0;
    T _w2 = 0;
    for (int i = 0; i < nd; i++) {
      _w0 += dphi0[q * nd + i] * u[i];
      _w1 += dphi1[q * nd + i] * u[i];
      _w2 += dphi2[q * nd + i] * u[i];
    }
    w[q] = _w0;
    w[q + nq] = _w1;
    w[q + 2 * nq] = _w2;
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_unfused(T __restrict__ *dphi0,
                                         T __restrict__ *dphi1,
                                         T __restrict__ *dphi2,
                                         T __restrict__ *u, T __restrict__ *w) {
  for (int q = 0; q < nq; q++) {
    T _w0 = 0;
    for (int i = 0; i < nd; i++) {
      _w0 += dphi0[q * nd + i] * u[i];
    }
    w[q] = _w0;

    T _w1 = 0;
    for (int i = 0; i < nd; i++) {

      _w1 += dphi1[q * nd + i] * u[i];
    }
    w[q + nq] = _w1;

    T _w2 = 0;
    for (int i = 0; i < nd; i++) {
      _w2 += dphi2[q * nd + i] * u[i];
    }
    w[q + 2 * nq] = _w2;
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_openmp(T __restrict__ *phi, T __restrict__ *u,
                                        T __restrict__ *w) {
  for (int q = 0; q < nq; q++) {
    T w0 = 0;
#pragma omp simd reduction(+ : w0)
    for (int i = 0; i < nd; i++)
      w0 += phi[q * nd + i] * u[i];
    w[q] = w0;
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_openmp(T __restrict__ *dphi0,
                                        T __restrict__ *dphi1,
                                        T __restrict__ *dphi2,
                                        T __restrict__ *u, T __restrict__ *w) {
  for (int q = 0; q < nq; q++) {
    T _w0 = 0;
    T _w1 = 0;
    T _w2 = 0;
#pragma omp simd reduction(+ : _w0, _w1, _w2)
    for (int i = 0; i < nd; i++) {
      _w0 += dphi0[q * nd + i] * u[i];
      _w1 += dphi1[q * nd + i] * u[i];
      _w2 += dphi2[q * nd + i] * u[i];
    }
    w[q] = _w0;
    w[q + nq] = _w1;
    w[q + 2 * nq] = _w2;
  }
}

// ------------------------------------------------------------
template <typename T, int nq, int nd>
inline void compute_coefficients_blis(T __restrict__ *dphi0,
                                      T __restrict__ *dphi1,
                                      T __restrict__ *dphi2, T __restrict__ *u,
                                      T __restrict__ *w) {
  rntm_t rntm;
  bli_rntm_init(&rntm);
  bli_rntm_set_num_threads(1, &rntm);

  obj_t _dphi0, _dphi1, _dphi2, _u, _w0, _w1, _w2;

  num_t dt = BLIS_DOUBLE;
  dim_t m = nq;
  dim_t n = nd;
  inc_t rs = nd;
  inc_t cs = 1;

  // Create tabulated basis objects
  bli_obj_create_with_attached_buffer(dt, m, n, dphi0, rs, cs, &_dphi0);
  bli_obj_create_with_attached_buffer(dt, m, n, dphi1, rs, cs, &_dphi1);
  bli_obj_create_with_attached_buffer(dt, m, n, dphi2, rs, cs, &_dphi2);

  // create blis objects for vectors
  bli_obj_create_with_attached_buffer(dt, n, 1, u, 1, 1, &_u);
  bli_obj_create_with_attached_buffer(dt, m, 1, w, 1, 1, &_w0);
  bli_obj_create_with_attached_buffer(dt, m, 1, w + nq, 1, 1, &_w1);
  bli_obj_create_with_attached_buffer(dt, m, 1, w + 2 * nq, 1, 1, &_w2);

  // Set the scalars to use.
  obj_t *alpha = &BLIS_ONE;
  obj_t *beta = &BLIS_ZERO;

  // Compute
  bli_gemv_ex(alpha, &_dphi0, &_u, beta, &_w0, NULL, &rntm);
  bli_gemv_ex(alpha, &_dphi1, &_u, beta, &_w1, NULL, &rntm);
  bli_gemv_ex(alpha, &_dphi2, &_u, beta, &_w2, NULL, &rntm);
}

// ------------------------------------------------------------
// Compute coefficient mass action
template <typename T, int nq, int nd>
inline void compute_coefficients_blis(T __restrict__ *phi, T __restrict__ *u,
                                      T __restrict__ *w) {
  obj_t _phi, _u, _w;
  num_t dt = BLIS_DOUBLE;
  dim_t m = nq;
  dim_t n = nd;
  inc_t rs = nd;
  inc_t cs = 1;

  bli_obj_create_with_attached_buffer(dt, m, n, phi, rs, cs, &_phi);
  bli_obj_create_with_attached_buffer(dt, n, 4, u, 4, 1, &_u);
  bli_obj_create_with_attached_buffer(dt, m, 4, w, 4, 1, &_w);

  // Set the scalars to use.
  obj_t *alpha = &BLIS_ONE;
  obj_t *beta = &BLIS_ZERO;

  bli_gemm(alpha, &_phi, &_u, beta, &_w);
}