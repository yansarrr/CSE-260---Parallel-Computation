

## PA1: SIMD-Optimized Matrix Multiplication

**Goal:**  
Optimize double-precision general matrix-matrix multiplication (`DGEMM`) using SIMD instructions (ARM SVE) and register blocking.

**Techniques Implemented:**
- Custom microkernel using SVE intrinsics (broadcast-based).
- Matrix packing for cache-efficient memory access.
- Parameter tuning (`MR`, `NR`, `MC`, `KC`, `NC`) for cache fit and performance.
- Manual vectorization exploiting register reuse and instruction-level parallelism.

**Best Performance Results:**
- **Peak GFLOPS:** ~24.74 (N = 513)
- **Performance trend:** Optimal performance occurs when matrix dimensions are multiples of 8, due to register block alignment with 256-bit SIMD vectors.
- **Irregular dips:** Observed for sizes like 513, 1025, 2049 due to vector misalignment and increased register/memory overhead.

---

## PA2: CUDA GPU Matrix Multiplication

**Goal:**  
Implement high-performance square matrix multiplication on NVIDIA GPUs using CUDA and shared memory optimizations.

**Techniques Implemented:**
- Tiling and thread blocking.
- Shared memory staging for reduced global memory access.
- Loop unrolling to enhance instruction throughput.
- 2D tiling (TileScale MxN) with thread responsibility scaling.
- Handling of edge cases for non-divisible matrix sizes.

**Best Performance Results (GFLOPS by Matrix Size):**

| N     | Peak GFLOPS | Optimal Block Size |
|-------|--------------|--------------------|
| 256   | 878.5        | 16×16              |
| 512   | 2972.7       | 32×32              |
| 1024  | 3144.8       | 32×32              |
| 2048  | 3370.0       | 32×32              |

**Speedups Compared to Naive GPU Version:**
- Over **5×** speedup with shared memory tiling and loop unrolling at N=2048.

**Roofline Peak Estimate:** 7680 GFLOPS theoretical; achieved ~54 ops/word at N=2048 using shared memory.

---

## PA3: MPI-Based 2D Wave Equation Simulation

**Goal:**  
Simulate a time-evolving 2D wave equation using domain decomposition and message passing with MPI.

**Techniques Implemented:**
- 5-point stencil update method.
- Ghost cell exchange via non-blocking `MPI_Isend`/`MPI_Irecv`.
- Overlapping communication and computation to reduce idle time.
- Absorbing boundary conditions for physical realism.
- Sine wave stimuli injection and dynamic load balancing.

**Best Performance Results (Strong Scaling):**

| Grid Size | Cores | GFlops/sec |
|-----------|-------|-------------|
| 500×500   | 16    | 14.78       |
| 1000×1000 | 64    | 19.45       |
| 2000×2000 | 128   | 31.75       |
| 4000×4000 | 384   | 88.69       |

**Observations:**
- Communication overhead stayed low due to effective use of non-blocking operations.
- Balanced core geometries (e.g. 16×16) yielded best performance.
- Weak scaling efficiency remained near ideal when problem size scaled with core count.


