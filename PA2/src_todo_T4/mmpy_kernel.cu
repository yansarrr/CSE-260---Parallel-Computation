// ;-*- mode: c;-*-
// Matrix multiply device code
#include <assert.h>
#include <math.h>
#include "../src/utils.h"
#include "../src/types.h"
#include "mytypes.h"
using namespace std;

#include <stdio.h>

#if defined(NAIVE)
__global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {

    int I =  blockIdx.y*blockDim.y + threadIdx.y;
    int J =  blockIdx.x*blockDim.x + threadIdx.x;

    if((I < N) && (J < N)){
        _FTYPE_ _c = 0;
        for (unsigned int k = 0; k < N; k++) {
            _FTYPE_ a = A[I * N + k];
            _FTYPE_ b = B[k * N + J];
            _c += a * b;
        }
        C[I * N + J] = _c;
    }
}

#elif defined(SMEM)
//You should be changing the kernel here for the non naive implementation.
__global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {

    extern __shared__ _FTYPE_ sMem[];
    _FTYPE_ (*As)[TILEDIM_K] = (_FTYPE_ (*)[TILEDIM_K])sMem;
    _FTYPE_ (*Bs)[TILEDIM_N] = (_FTYPE_ (*)[TILEDIM_N])&sMem[TILEDIM_M * TILEDIM_K];

    int ty = threadIdx.y;
    int tx = threadIdx.x;
    int by = blockIdx.y;
    int bx = blockIdx.x;

    int I =  by*TILEDIM_M + ty;
    int J =  bx*TILEDIM_N + tx;

    float Cij = 0.0;

    for (int kk = 0; kk < N; kk+=TILEDIM_K)
    {
        if (I < N && kk+tx < N){
            As[ty][tx] = A[I*N + kk + tx];
        }else{
            As[ty][tx] = 0;
        }

        if (J < N && kk+ty < N){
            Bs[ty][tx] = B[(kk + ty)*N + J];
        }else{
            Bs[ty][tx] = 0;
        }

        __syncthreads();

        for (int k = 0; k < TILEDIM_K; k++)
        {
            /* code */
            Cij += As[ty][k] * Bs[k][tx];
        }
        __syncthreads();
    }

    if (I < N && J < N)
    {
        C[I*N + J] = Cij;
    }
    
}

#else
//You should be changing the kernel here for the non naive implementation.
__global__ void matMul(int N, _FTYPE_ *C, _FTYPE_ *A, _FTYPE_ *B) {

    extern __shared__ _FTYPE_ sMem[];
    _FTYPE_ (*As) = (_FTYPE_ (*))sMem;
    _FTYPE_ (*Bs) = (_FTYPE_ (*))&sMem[TILEDIM_M * TILEDIM_K];

    const uint ty = threadIdx.y;
    const uint tx = threadIdx.x;
    const uint by = blockIdx.y;
    const uint bx = blockIdx.x;

    _FTYPE_ Cijs[TILESCALE_M * TILESCALE_N] = {0.0};
    _FTYPE_ A_buf[TILESCALE_M] = {0.0};
    _FTYPE_ B_buf[TILESCALE_N] = {0.0};

    // assign (by, bx) threads to load values into As (T_M, T_K)
    const uint Ay = (ty * blockDim.x + tx) / TILEDIM_K;
    const uint Ax = (ty * blockDim.x + tx) % TILEDIM_K;
    // assign (by, bx) threads to load values into Bs (T_K, T_N)
    const uint By = (ty * blockDim.x + tx) / TILEDIM_N;
    const uint Bx = (ty * blockDim.x + tx) % TILEDIM_N;

    const uint As_offset = (blockDim.x * blockDim.y) / TILEDIM_K;
    const uint Bs_offset = (blockDim.x * blockDim.y) / TILEDIM_N;

    #pragma unroll
    for (uint kk = 0; kk < N; kk+=TILEDIM_K){
        #pragma unroll
        for (uint i = 0; i < TILEDIM_M; i+=As_offset){
            if (((by * TILEDIM_M + Ay+i) < N) && ((kk + Ax) < N)){
                As[(Ay + i) * TILEDIM_K + Ax] = A[(by * TILEDIM_M + Ay + i) * N + kk + Ax];
            }else{
                As[(Ay + i) * TILEDIM_K + Ax] = 0;
            }
        }
        
        #pragma unroll
        for (uint i = 0; i < TILEDIM_K; i+=Bs_offset){
            if (((kk + By + i) < N) && ((bx*TILEDIM_N + Bx) < N)){
                Bs[(By + i) * TILEDIM_N + Bx] = B[(kk + By + i) * N + bx*TILEDIM_N + Bx];
            }else{
                Bs[(By + i) * TILEDIM_N + Bx] = 0;
            }
        }

        __syncthreads();

        #pragma unroll
        for (uint k = 0; k < TILEDIM_K; ++k)
        {
            #pragma unroll
            for (uint i = 0; i < TILESCALE_M; ++i){
                // seperate
                //A_buf[i] = As[(ty + (i * blockDim.y)) * TILEDIM_K + k];
                // adjacent
                A_buf[i] = As[((ty * TILESCALE_M) + i) * TILEDIM_K + k];
            }
            #pragma unroll
            for (uint i = 0; i < TILESCALE_N; ++i){
                // seperate
                //B_buf[i] = Bs[(k * TILEDIM_N) + tx + (i * blockDim.x)];
                // adjacent
                B_buf[i] = Bs[(k * TILEDIM_N) + (tx * TILESCALE_N) + i];
            }
            
            #pragma unroll
            for (uint i = 0; i < TILESCALE_M; ++i){
                #pragma unroll
                for (uint j = 0; j < TILESCALE_N; ++j){
                    Cijs[i * TILESCALE_M + j] += A_buf[i] * B_buf[j];
                }
            }
        }
        __syncthreads();
    }

    #pragma unroll
    for (uint i = 0; i < TILESCALE_M; ++i){
        #pragma unroll
        for (uint j = 0; j < TILESCALE_N; ++j){
            // seperate
            //if ((by*TILEDIM_M + ty + (i * blockDim.y)) < N && (bx*TILEDIM_N + tx + (j * blockDim.x)) < N){
            //    C[(by*TILEDIM_M + ty + (i * blockDim.y)) * N + (bx*TILEDIM_N + tx + (j * blockDim.x))] = Cijs[i * TILESCALE_M + j];
            //}
            // adjacent
            if ((by*TILEDIM_M + (ty * TILESCALE_M) + i) < N && (bx*TILEDIM_N + (tx * TILESCALE_N) + j) < N){
                C[(by*TILEDIM_M + (ty * TILESCALE_M) + i) * N + (bx*TILEDIM_N + (tx * TILESCALE_N) + j)] = Cijs[i * TILESCALE_M + j];
            }
        }
    }
}
#endif