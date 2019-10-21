#ifndef USED_MATH_HPP
#define USED_MATH_HPP

#include <math.h>
#include <stdint.h>
#include <cmath>
#include <cstring>

namespace pv {

constexpr double pi = 3.14159265358979323846;
constexpr double pi2 = 1.57079632679489661923;


template <typename scalar_t>
static constexpr scalar_t used_sin(scalar_t x) {
    return std::sin(x);
}
template <typename scalar_t>
static constexpr scalar_t used_cos(scalar_t x) {
    return std::cos(x);
}
template <typename scalar_t>
static constexpr scalar_t used_sqrt(scalar_t x) {
    return std::sqrt(x);
}

//https://www.dsprelated.com/showarticle/1052.php
static float used_atan2(float y, float x) {
    const float n1 = 0.97239411f;
    const float n2 = -0.19194795f;
    float result = 0.0f;
    if (x != 0.0f) {
        const union {
            float flVal;
            uint32_t nVal;
        } tYSign = { y };
        const union {
            float flVal;
            uint32_t nVal;
        } tXSign = { x };
        if (fabsf(x) >= fabsf(y)) {
            union {
                float flVal;
                uint32_t nVal;
            } tOffset = { static_cast<float>(pi) };
            // Add or subtract PI based on y's sign.
            tOffset.nVal |= tYSign.nVal & 0x80000000u;
            // No offset if x is positive, so multiply by 0 or based on x's sign.
            tOffset.nVal *= tXSign.nVal >> 31;
            result = tOffset.flVal;
            const float z = y / x;
            result += (n1 + n2 * z * z) * z;
        } else // Use atan(y/x) = pi/2 - atan(x/y) if |y/x| > 1.
        {
            union {
                float flVal;
                uint32_t nVal;
            } tOffset = { static_cast<float>(pi2) };
            // Add or subtract PI/2 based on y's sign.
            tOffset.nVal |= tYSign.nVal & 0x80000000u;
            result = tOffset.flVal;
            const float z = x / y;
            result -= (n1 + n2 * z * z) * z;
        }
    } else if (y > 0.0f) {
        result = static_cast<float>(pi2);
    } else if (y < 0.0f) {
        result = -static_cast<float>(pi2);
    }
    return result;
}

template <typename scalar_t>
static scalar_t used_fmod(scalar_t x, scalar_t y) {
    unsigned times = static_cast<unsigned> (x / y);
    scalar_t remain = x - times * y;
    return remain;
}

template<typename scalar_t>
static void used_fill(scalar_t *p_dest, scalar_t value, uint32_t size) {
    std::memset(p_dest, value, size);
}
} // end namespace pv

#endif // USED_MATH_HPP
