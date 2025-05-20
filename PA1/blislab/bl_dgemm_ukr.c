#include "bl_config.h"
#include "bl_dgemm_kernel.h"
#include <arm_sve.h>

void bl_dgemm_ukr(int k,
                  int m,
                  int n,
                  double *a,
                  double *b,
                  double *c,
                  unsigned long long ldc,
                  aux_t *data) {
    svbool_t n0pred, n1pred;

    if (n > 4) {
        n0pred = svwhilelt_b64(0, 4);
        n1pred = svwhilelt_b64(0, n - 4);
    } else {
        n0pred = svwhilelt_b64(0, n);
    }
    
    svfloat64_t c0x0, c0x1, c1x0, c1x1, c2x0, c2x1, c3x0, c3x1, 
                c4x0, c4x1, c5x0, c5x1, c6x0, c6x1, c7x0, c7x1, 
                c8x0, c8x1, c9x0, c9x1, c10x0, c10x1, c11x0, c11x1, 
                c12x0, c12x1, c13x0, c13x1;
                
    int p;

    // Load C vectors
    if (m > 0) {
        c0x0 = svld1_f64(n0pred, c + 0 * ldc);
        if (n > 4) { 
            c0x1 = svld1_f64(n1pred, c + 0 * ldc + 4);
        }
    }
    if (m > 1) {
        c1x0 = svld1_f64(n0pred, c + 1 * ldc);
        if (n > 4) {
            c1x1 = svld1_f64(n1pred, c + 1 * ldc + 4);
        }
    }
    if (m > 2) {
        c2x0 = svld1_f64(n0pred, c + 2 * ldc);
        if (n > 4) {
            c2x1 = svld1_f64(n1pred, c + 2 * ldc + 4);
        }
    }
    if (m > 3) {
        c3x0 = svld1_f64(n0pred, c + 3 * ldc);
        if (n > 4) {
            c3x1 = svld1_f64(n1pred, c + 3 * ldc + 4);
        }
    }
    if (m > 4) {
        c4x0 = svld1_f64(n0pred, c + 4 * ldc);
        if (n > 4) {
            c4x1 = svld1_f64(n1pred, c + 4 * ldc + 4);
        }
    }
    if (m > 5) {
        c5x0 = svld1_f64(n0pred, c + 5 * ldc);
        if (n > 4) {
            c5x1 = svld1_f64(n1pred, c + 5 * ldc + 4);
        }
    }
    if (m > 6) {
        c6x0 = svld1_f64(n0pred, c + 6 * ldc);
        if (n > 4) {
            c6x1 = svld1_f64(n1pred, c + 6 * ldc + 4);
        }
    }
    if (m > 7) {
        c7x0 = svld1_f64(n0pred, c + 7 * ldc);
        if (n > 4) {
            c7x1 = svld1_f64(n1pred, c + 7 * ldc + 4);
        }
    }
    if (m > 8) {
        c8x0 = svld1_f64(n0pred, c + 8 * ldc);
        if (n > 4) {
            c8x1 = svld1_f64(n1pred, c + 8 * ldc + 4);
        }
    }
    if (m > 9) {
        c9x0 = svld1_f64(n0pred, c + 9 * ldc);
        if (n > 4) {
            c9x1 = svld1_f64(n1pred, c + 9 * ldc + 4);
        }
    }
    if (m > 10) {
        c10x0 = svld1_f64(n0pred, c + 10 * ldc);
        if (n > 4) {
            c10x1 = svld1_f64(n1pred, c + 10 * ldc + 4);
        }
    }
    if (m > 11) {
        c11x0 = svld1_f64(n0pred, c + 11 * ldc);
        if (n > 4) {
            c11x1 = svld1_f64(n1pred, c + 11 * ldc + 4);
        }
    }
    if (m > 12) {
        c12x0 = svld1_f64(n0pred, c + 12 * ldc);
        if (n > 4) {
            c12x1 = svld1_f64(n1pred, c + 12 * ldc + 4);
        }
    }
    if (m > 13) {
        c13x0 = svld1_f64(n0pred, c + 13 * ldc);
        if (n > 4) {
            c13x1 = svld1_f64(n1pred, c + 13 * ldc + 4);
        }
    }



    for (p = 0; p < k; p++) {
        svfloat64_t bx0 = svld1_f64(svptrue_b64(), b + p * DGEMM_NR);
        svfloat64_t bx1;

        if (n > 4) {
            bx1 = svld1_f64(svptrue_b64(), b + p * DGEMM_NR + 4);
        }

        if (m > 0) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 0]);
            c0x0 = svmla_f64_m(n0pred, c0x0, ax, bx0);
            if (n > 4) {
                c0x1 = svmla_f64_m(n1pred, c0x1, ax, bx1);
            }
        }

      if (m > 1) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 1]);
            c1x0 = svmla_f64_m(n0pred, c1x0, ax, bx0);
            if (n > 4) {
                c1x1 = svmla_f64_m(n1pred, c1x1, ax, bx1);
            }
        }
        if (m > 2) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 2]);
            c2x0 = svmla_f64_m(n0pred, c2x0, ax, bx0);
            if (n > 4) {
                c2x1 = svmla_f64_m(n1pred, c2x1, ax, bx1);
            }
        }
        if (m > 3) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 3]);
            c3x0 = svmla_f64_m(n0pred, c3x0, ax, bx0);
            if (n > 4) {
                c3x1 = svmla_f64_m(n1pred, c3x1, ax, bx1);
            }
        }
        if (m > 4) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 4]);
            c4x0 = svmla_f64_m(n0pred, c4x0, ax, bx0);
            if (n > 4) {
                c4x1 = svmla_f64_m(n1pred, c4x1, ax, bx1);
            }
        }
        if (m > 5) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 5]);
            c5x0 = svmla_f64_m(n0pred, c5x0, ax, bx0);
            if (n > 4) {
                c5x1 = svmla_f64_m(n1pred, c5x1, ax, bx1);
            }
        }
        if (m > 6) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 6]);
            c6x0 = svmla_f64_m(n0pred, c6x0, ax, bx0);
            if (n > 4) {
                c6x1 = svmla_f64_m(n1pred, c6x1, ax, bx1);
            }
        }
        if (m > 7) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 7]);
            c7x0 = svmla_f64_m(n0pred, c7x0, ax, bx0);
            if (n > 4) {
                c7x1 = svmla_f64_m(n1pred, c7x1, ax, bx1);
            }
        }
        if (m > 8) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 8]);
            c8x0 = svmla_f64_m(n0pred, c8x0, ax, bx0);
            if (n > 4) {
                c8x1 = svmla_f64_m(n1pred, c8x1, ax, bx1);
            }
        }
        if (m > 9) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 9]);
            c9x0 = svmla_f64_m(n0pred, c9x0, ax, bx0);
            if (n > 4) {
                c9x1 = svmla_f64_m(n1pred, c9x1, ax, bx1);
            }
        }
        if (m > 10) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 10]);
            c10x0 = svmla_f64_m(n0pred, c10x0, ax, bx0);
            if (n > 4) {
                c10x1 = svmla_f64_m(n1pred, c10x1, ax, bx1);
            }
        }
        if (m > 11) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 11]);
            c11x0 = svmla_f64_m(n0pred, c11x0, ax, bx0);
            if (n > 4) {
                c11x1 = svmla_f64_m(n1pred, c11x1, ax, bx1);
            }
        }
        if (m > 12) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 12]);
            c12x0 = svmla_f64_m(n0pred, c12x0, ax, bx0);
            if (n > 4) {
                c12x1 = svmla_f64_m(n1pred, c12x1, ax, bx1);
            }
        }
        if (m > 13) {
            svfloat64_t ax = svdup_f64(a[p * DGEMM_MR + 13]);
            c13x0 = svmla_f64_m(n0pred, c13x0, ax, bx0);
            if (n > 4) {
                c13x1 = svmla_f64_m(n1pred, c13x1, ax, bx1);
            }
        }
    }

    if (m > 0) {
        svst1_f64(n0pred, c + 0 * ldc, c0x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 0 * ldc + 4, c0x1);
        }
    }
    if (m > 1) {
        svst1_f64(n0pred, c + 1 * ldc, c1x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 1 * ldc + 4, c1x1);
        }
    }
    if (m > 2) {
        svst1_f64(n0pred, c + 2 * ldc, c2x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 2 * ldc + 4, c2x1);
        }
    }
    if (m > 3) {
        svst1_f64(n0pred, c + 3 * ldc, c3x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 3 * ldc + 4, c3x1);
        }
    }
    if (m > 4) {
        svst1_f64(n0pred, c + 4 * ldc, c4x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 4 * ldc + 4, c4x1);
        }
    }
    if (m > 5) {
        svst1_f64(n0pred, c + 5 * ldc, c5x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 5 * ldc + 4, c5x1);
        }
    }
    if (m > 6) {
        svst1_f64(n0pred, c + 6 * ldc, c6x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 6 * ldc + 4, c6x1);
        }
    }
    if (m > 7) {
        svst1_f64(n0pred, c + 7 * ldc, c7x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 7 * ldc + 4, c7x1);
        }
    }
    if (m > 8) {
        svst1_f64(n0pred, c + 8 * ldc, c8x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 8 * ldc + 4, c8x1);
        }
    }
    if (m > 9) {
        svst1_f64(n0pred, c + 9 * ldc, c9x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 9 * ldc + 4, c9x1);
        }
    }
    if (m > 10) {
        svst1_f64(n0pred, c + 10 * ldc, c10x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 10 * ldc + 4, c10x1);
        }
    }
    if (m > 11) {
        svst1_f64(n0pred, c + 11 * ldc, c11x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 11 * ldc + 4, c11x1);
        }
    }
    if (m > 12) {
        svst1_f64(n0pred, c + 12 * ldc, c12x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 12 * ldc + 4, c12x1);
        }
    }
    if (m > 13) {
        svst1_f64(n0pred, c + 13 * ldc, c13x0);
        if (n > 4) {
            svst1_f64(n1pred, c + 13 * ldc + 4, c13x1);
        }
    }
}




