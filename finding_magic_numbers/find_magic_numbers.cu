#include <cstdint>
#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>
#include <tuple>
#include <typeinfo>

#ifndef TARGET_ROOT
#define TARGET_ROOT 2
#endif

#ifndef FIND_DOUBLE

#define TARGET_FTYPE float
#define TARGET_ITYPE uint32_t
#define MANTISSA_SIZE ((1 << 24) - 1)
#define INITIAL_MAGIC (((127 - 127 / TARGET_ROOT) - 1) << 23)
#define SKIP_AMOUNT_MUL 1

#else

#define TARGET_FTYPE double
#define TARGET_ITYPE uint64_t
#define MANTISSA_SIZE ((IntT(1) << 53) - 1)
#define INITIAL_MAGIC ((IntT(1023 - 1023 / TARGET_ROOT) - 1) << 52)
#define SKIP_AMOUNT_MUL 179424673

#endif

template <uint32_t n, typename T>
__device__ inline T ct_pow(const T x)
{
    if (n == 0) {
        return 1;
    } else {
        const T part = ct_pow<n / 2, T>(x);
        return ((n & 1) ? (x * part * part) : (part * part));
    }
}

template <uint32_t n, typename T>
__device__ T halley_step(const T x0, const T value)
{
    const T fx = ct_pow<n>(x0) - value;
    const T fpx = n * ct_pow<n - 1>(x0);
    const T fppx = n * (n - 1) * ct_pow<n - 2>(x0);
    const T numer = 2 * fx * fpx;
    const T denom = 2 * fpx * fpx - fx * fppx;
    const T x1 = x0 - (numer / denom);
    return x1;
}

template <uint32_t n, typename T>
__device__ T newton_step(const T x0, const T value)
{
    // x1 = x0 - (f(x0)-y)/f'(x0)
    const T x1 = x0 - ((ct_pow<n>(x0) - value) / (n * ct_pow<n - 1>(x0)));
    return x1;
}

template <typename FloatT, typename IntT>
struct SetInitialErrorLevel {
    __device__ void operator()(const IntT index)
    {
        magic_max_error[index] = 10000;
    }
    FloatT* __restrict__ magic_max_error;
};

template <typename FloatT, typename IntT, uint32_t n>
struct FindErrorForMagics {
    __device__ void operator()(const IntT index)
    {
        IntT magic = INITIAL_MAGIC + magic_offset + index;

        FloatT error = 0.0;

        union {
            FloatT f;
            IntT i;

        } packed;

        for (FloatT root_value = 0.001; root_value < 4.0001; root_value += 0.001) {
            const FloatT powered = ct_pow<n>(root_value);

            packed.f = powered;
            packed.i = magic + (packed.i / n);

            const FloatT approx = halley_step<n>(packed.f, powered);

            const FloatT current_error = abs((approx - root_value) / root_value);
            error = max(current_error, error);
        }

        // Reuse memory
        if (error < magic_max_error[index]) {
            magic_max_error[index] = error;
            magics[index] = magic;
        }
    }

    const IntT magic_offset;
    FloatT* __restrict__ magic_max_error;
    IntT* __restrict__ magics;
};

std::pair<TARGET_FTYPE, TARGET_ITYPE> find_magic_number(void)
{
    using FloatT = TARGET_FTYPE;
    using IntT = TARGET_ITYPE;

    using KernelT = FindErrorForMagics<FloatT, IntT, TARGET_ROOT>;

    const IntT max_per_round = (1 << 14);

    // Allocate memory
    FloatT* local_error = new FloatT[max_per_round];
    IntT* local_magics = new IntT[max_per_round];

    FloatT* gpu_errors;
    IntT* gpu_magics;

    cudaMalloc(&gpu_errors, sizeof(FloatT) * max_per_round);
    cudaMalloc(&gpu_magics, sizeof(IntT) * max_per_round);

    // Set the initial error level to be high
    thrust::for_each(thrust::counting_iterator<IntT>(0),
        thrust::counting_iterator<IntT>(max_per_round),
        SetInitialErrorLevel<FloatT, IntT>{ gpu_errors });

    // Find the best magic
    for (IntT magic_offset = 0; magic_offset < MANTISSA_SIZE;
         magic_offset += (max_per_round * SKIP_AMOUNT_MUL)) {
        thrust::for_each(thrust::counting_iterator<IntT>(0),
            thrust::counting_iterator<IntT>(max_per_round),
            KernelT{ magic_offset, gpu_errors, gpu_magics }

            );
        // Print progress
        std::cerr << std::fixed << ((100.0 * magic_offset) / MANTISSA_SIZE)
                  << "%                    \r";
    }
    std::cerr << "\n";

    // Copy stuff back
    cudaMemcpy(local_error, gpu_errors, sizeof(FloatT) * max_per_round,
        cudaMemcpyDeviceToHost);
    cudaMemcpy(local_magics, gpu_magics, sizeof(FloatT) * max_per_round,
        cudaMemcpyDeviceToHost);

    // Find the best magic
    FloatT max_error = 1000;
    IntT current_magic = 0;

    for (IntT i = 0; i < max_per_round; ++i) {
        if (local_error[i] < max_error) {
            max_error = local_error[i];
            current_magic = local_magics[i];
        }
    }

    cudaFree(&gpu_errors);
    cudaFree(&gpu_magics);
    delete[] local_error;
    delete[] local_magics;

    return std::make_pair(max_error, current_magic);
}

int main(void)
{
    std::cout << "Finding x^" << TARGET_ROOT << " for type ";
#ifdef FIND_DOUBLE
    std::cout << "Double\n";
#else
    std::cout << "Float\n";
#endif
    auto result = find_magic_number();
    std::cout << "Relative error: " << std::fixed << result.first << "\n";
    std::cout << "Magic:          0x" << std::hex << result.second << " + (i / "
              << std::dec << TARGET_ROOT << ")\n";
}