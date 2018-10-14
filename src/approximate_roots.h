#pragma once

#include <cmath>
#include <cstdint>

namespace approximate_roots {

// Perform x^n, where n is an unsigned integer
template <uint32_t n, typename T>
inline T ct_pow(const T x) {
    if (n == 0) {
        return 1;
    } else {
        const T part = ct_pow<n / 2, T>(x);
        return ((n & 1) ? (x * part * part) : (part * part));
    }
}

// Perform a single step of Halleys method, where the
// function is f(x) = x^n - value
template <uint32_t n, typename T>
T halley_step(const T x0, const T value) {
    const T fx = ct_pow<n>(x0) - value;
    const T fpx = n * ct_pow<n - 1>(x0);
    const T fppx = n * (n - 1) * ct_pow<n - 2>(x0);
    const T numer = 2 * fx * fpx;
    const T denom = 2 * fpx * fpx - fx * fppx;
    const T x1 = x0 - (numer / denom);
    return x1;
}

// Perform a fast approximation of a float, followed
// by a single Halley step
template <uint32_t magic, uint32_t n>
float float_approx_with_halley_step(const float x) {
    union {
        float f;
        uint32_t i;
    } packed = {.f = x};
    packed.i = magic - (packed.i / n);
    return halley_step<n>(packed.f, x);
}

// Perform a fast approximation of a double, followed
// by a single Halley step
template <uint64_t magic, uint32_t n>
double double_approx_with_halley_step(const double x) {
    union {
        double f;
        uint64_t i;
    } packed = {.f = x};
    packed.i = magic - (packed.i / n);
    return halley_step<n>(packed.f, x);
}

// Approximations for floats

inline float approx_sqrt(const float x) {
    return float_approx_with_halley_step<0x1fbb67ad, 2>(x);
}

inline float approx_cbrt(const float x) {
    return float_approx_with_halley_step<0x2a511949, 3>(x);
}

inline float approx_forth_root(const float x) {
    return float_approx_with_halley_step<0x2f9b5088, 4>(x);
}

inline float approx_fifth_root(const float x) {
    return float_approx_with_halley_step<0x32c82ec7, 5>(x);
}

inline float approx_sixth_root(const float x) {
    return float_approx_with_halley_step<0x34e5e317, 6>(x);
}

inline float approx_seventh_root(const float x) {
    return float_approx_with_halley_step<0x3668ef49, 7>(x);
}

inline float approx_eighth_root(const float x) {
    return float_approx_with_halley_step<0x378b0a48, 8>(x);
}

inline float approx_ninth_root(const float x) {
    return float_approx_with_halley_step<0x38714eaf, 9>(x);
}

inline float approx_tenth_root(const float x) {
    return float_approx_with_halley_step<0x391cabf0, 10>(x);
}

// Approximations for doubles

inline double approx_sqrt(const double x) {
    return double_approx_with_halley_step<0x1ff76cf48689feb3, 2>(x);
}

inline double approx_cbrt(const double x) {
    return double_approx_with_halley_step<0x2a9f77a7a61a7e7a, 3>(x);
}

inline double approx_forth_root(const double x) {
    return double_approx_with_halley_step<0x2ff36a476a29c002, 4>(x);
}

inline double approx_fifth_root(const double x) {
    return double_approx_with_halley_step<0x3325d3e38d3b8000, 5>(x);
}

inline double approx_sixth_root(const double x) {
    return double_approx_with_halley_step<0x3547679b9fb94000, 6>(x);
}

inline double approx_seventh_root(const double x) {
    return double_approx_with_halley_step<0x36cd1e16377a0000, 7>(x);
}

inline double approx_eighth_root(const double x) {
    return double_approx_with_halley_step<0x37f16397f5297fd4, 8>(x);
}

inline double approx_ninth_root(const double x) {
    return double_approx_with_halley_step<0x38d4b87bc4113ff5, 9>(x);
}

inline double approx_tenth_root(const double x) {
    return double_approx_with_halley_step<0x398a9464ab857fdf, 10>(x);
}

} // namespace fast_root
