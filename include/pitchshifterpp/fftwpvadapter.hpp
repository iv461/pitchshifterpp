#ifndef FFTWPVADAPTER_H
#define FFTWPVADAPTER_H

#include <fftw/fftw3.h>
#include <stdint.h>
/**
 * @brief FFTWPVAdapter adapter class to use FFTW with
 * PhaseVocoder, implicit interface containing @ref rfft() and @ref rifft()
 * methods.
 * The arrays allocated are not used as in the @ref rfft() and @ref rifft()
 * methods the in and out arrays are passed, but are needed if the FFTW
 * plans are created without FFTW_ESTIMATE.
 * The fftw_execute* functions are not able to determine if the momory
 * alignment of the input and output arrays are correct
 * so unless FFTW_NO_SIMD is specified on plan creation
 * the memory must be 16-byte aligned.
 */
class FFTWPVAdapter
{
public:
    FFTWPVAdapter(uint32_t size) :
        m_fft_size(size),
        m_rfft_in(fftw_alloc_real(size)),
        m_rifft_out(fftw_alloc_real(size)),
        m_rifft_in(fftw_alloc_complex(size / 2 + 1)),
        m_rfft_out(fftw_alloc_complex(size / 2 + 1)),
        m_r2c_plan(fftw_plan_dft_r2c_1d(static_cast<int>(size),
                                        m_rfft_in, m_rfft_out, FFTW_ESTIMATE)),
        m_c2r_plan(fftw_plan_dft_c2r_1d(static_cast<int>(size),
                                        m_rifft_in, m_rifft_out,
                                        FFTW_ESTIMATE))
    {
    }

    FFTWPVAdapter(FFTWPVAdapter const &other) = delete;
    FFTWPVAdapter &operator=(FFTWPVAdapter const &other) = delete;

    void rfft(double *in, double *out) {
        fftw_execute_dft_r2c(m_r2c_plan, in,
                             reinterpret_cast<fftw_complex *>(out));
    }

    /**
     * @brief rifft inverse fft, is allowed to destroy input array
     * @param in
     * @param out
     */
    void rifft(double *in, double *out) {
        fftw_execute_dft_c2r(m_c2r_plan, reinterpret_cast<fftw_complex *>(in),
                             out);
        // normalize
        for(size_t i = 0; i < m_fft_size ; i++) {
            out[i] *= (1. / static_cast<double>(m_fft_size));
        }
    }

    ~FFTWPVAdapter() {
        fftw_destroy_plan(m_r2c_plan);
        fftw_destroy_plan(m_c2r_plan);
        fftw_free(m_rfft_out);
        fftw_free(m_rifft_in);
        fftw_free(m_rfft_in);
        fftw_free(m_rifft_out);
    }

private:
    uint32_t m_fft_size;

    double *m_rfft_in;
    double *m_rifft_out;
    fftw_complex *m_rifft_in;
    fftw_complex *m_rfft_out;

    fftw_plan m_r2c_plan;
    fftw_plan m_c2r_plan;
};

#endif // FFTWPVADAPTER_H
