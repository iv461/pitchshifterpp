#ifndef USED_MATH_HPP
#define USED_MATH_HPP
/**
  * @brief used math functions which can be changed if
  * e.g for different platfoms
*/

#include <math.h>
#include <stdint.h>
#include <cstring>

namespace pv {

constexpr double pi = 3.14159265358979323846;
constexpr double pi2 = 1.57079632679489661923;

template <typename scalar_t>
static scalar_t used_sin(scalar_t x) {
    return std::sin(x);
}
template <typename scalar_t>
static scalar_t used_cos(scalar_t x) {
    return std::cos(x);
}
template <typename scalar_t>
static scalar_t used_sqrt(scalar_t x) {
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

/**
 * @brief cmplx_cartesian_to_polar in place-capable
 * implementation of conversion from cartesian to
 * polar on a complex-valued vector
 */
template<typename scalar_t>
static void cmplx_cartesian_to_polar(scalar_t *out,
                                     scalar_t const *in,
                                     unsigned int length) {
    scalar_t imag_tmp, real_tmp, squared_tmp;
    for(unsigned i = 0; i < length; i++) {
        real_tmp = in[0];
        imag_tmp = in[1];
        in += 2;
        squared_tmp = real_tmp * real_tmp + imag_tmp * imag_tmp;
        out[0] = used_sqrt(squared_tmp);
        out[1] = used_atan2(imag_tmp, real_tmp);
        out += 2;
    }
}

/**
 * @brief cmplx_polar_to_cartesian in place-capable
 * implementation of conversion from polar to
 * cartesian on a complex-valued vector
 */
template <typename scalar_t>
static void cmplx_polar_to_cartesian(scalar_t *complex_out,
                                     scalar_t const *complex_in,
                                     unsigned int length) {
    scalar_t mag_tmp, phase_tmp;
    for(unsigned i = 0; i < length; i++) {
        mag_tmp = complex_in[0];
        phase_tmp = complex_in[1];
        complex_out[0] = mag_tmp * used_cos(phase_tmp);
        complex_out[1] = mag_tmp * used_sin(phase_tmp);
        complex_in += 2;
        complex_out += 2;
    }
}

} // end namespace pv

#endif // USED_MATH_HPP
