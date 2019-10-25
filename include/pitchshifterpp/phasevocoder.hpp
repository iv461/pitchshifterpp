#ifndef PHASEVOCODER_HPP
#define PHASEVOCODER_HPP

#include <cmath>
#include <stdint.h>

#include "used_math.hpp"

namespace pv {

template <typename scalar_t>
struct Complex {
    scalar_t real, imag;
};

template <typename scalar_t>
struct ComplexPolar {
    scalar_t magnitude, phase;
};

/**
 * @brief simple circular iterator for plain arrays, no checks
 */
template<typename T>
class circular_iterator {
public:
    circular_iterator(T *begin, uint32_t size) :
        m_begin(begin),  m_size(size) {}
    T& operator*() {
        return m_begin[m_counter];
    }
    T& operator[](uint32_t idx) const {
        return m_begin[(m_counter + (idx % m_size)) % m_size];
    }
    circular_iterator& operator+=(uint32_t idx) {
        auto inc_counter = m_counter + (idx % m_size);
        m_is_full = !m_is_full ? inc_counter >= m_size:true;
        m_counter = inc_counter % m_size;
        return *this;
    }
    circular_iterator& operator-=(uint32_t idx) {
        idx %= m_size;
        m_counter = m_counter < idx ? m_size - (idx - m_counter)
                                    : m_counter - idx;
        return *this;
    }
    circular_iterator& operator--(int) {
        m_counter = !m_counter ? m_size - 1: m_counter - 1;
        return *this;
    }
    circular_iterator& operator++(int) {
        auto inc_counter = m_counter + 1;
        m_is_full = !m_is_full ? inc_counter >= m_size:true;
        m_counter = inc_counter % m_size;
        return *this;
    }
    circular_iterator operator+(uint32_t num) const {
        circular_iterator copy(*this);
        copy += num;
        return copy;
    }
    circular_iterator operator-(uint32_t num) const {
        circular_iterator copy(*this);
        copy -= num;
        return copy;
    }
    bool is_full() const {
        return m_is_full;
    }
protected:
    T *m_begin;
    uint32_t m_size;
    uint32_t m_counter = 0;
    bool m_is_full = false;
};

/**
 * @brief wraps x to the inteval [-pi, pi[
 * @param x input to be wrapped
 */
template <typename scalar_t>
static scalar_t wrap_phase(scalar_t x) {
    constexpr auto pi_times_2 = static_cast<scalar_t>(2. * pi);
    return used_fmod(pi_times_2 +
                     used_fmod(static_cast<scalar_t>(x + pi), pi_times_2),
                     pi_times_2) - static_cast<scalar_t>(pi);
}

/**
 * @brief Hann window function
 */
template <typename scalar_t>
static scalar_t hann(scalar_t x) {
    return .5 * (1. - used_cos(static_cast<scalar_t>(pi) * 2 * x));
}

template <typename scalar_t>
static void resample_to_real
(circular_iterator<scalar_t> input, uint32_t input_length,
 scalar_t *output, uint32_t output_length, scalar_t ratio,
 unsigned max_out)
{
    for (uint32_t i = 1; i < output_length; i++) {
        scalar_t index = ratio * i;
        int32_t trunc_idx = static_cast<int32_t>(index);
        scalar_t trunc_idx_f = static_cast<scalar_t>(trunc_idx);
        scalar_t fac = index - trunc_idx_f;
        scalar_t a = input[trunc_idx];
        scalar_t b = input[trunc_idx + 1];
        scalar_t val = a + fac * (b - a);
        if((i - 1) == max_out) {
            return;
        }
        output[i - 1] = val;
    }
    output[(output_length - 1)] = input[(input_length - 1)];
}

static auto semitone_to_scale_factor(int32_t semitone) {
    return std::pow(2., semitone / 12.);
}

/**
 * @brief PhaseVocoder optimized for speed and low memory consumption,
 * currently implemented with 1.5 * window_size latency
 * @tparam window_size window size of FFT, is a template parameter  as
 * all needed memory is allocated at compile time;
 * @tparam max_syn_size this is equivalent to @ref window_size *
 *  maximum scaling factor which is planned to be used
 * @tparam fft_adapter_t non-virtual adapter and implicit interface type
 * for the used FFT
 * @tparam scalar_t scalar type
 */
template <uint32_t window_size, uint32_t max_syn_size,
          typename fft_adapter_t, typename scalar_t = double>
class PhaseVocoder
{
    static_assert (!(window_size % 2), "window size must be even");
public:
    PhaseVocoder(fft_adapter_t &fft_adapter, uint32_t overlap_factor) :
        m_overlap_factor(overlap_factor),
        m_ana_hop_size(window_size / overlap_factor),
        m_ana_hop_inv(1. / static_cast<scalar_t>(m_ana_hop_size)),
        m_fft_adapter(fft_adapter)
    {
        used_fill<scalar_t>(m_pre_ana_phase_frame, 0,
                            sizeof(m_pre_ana_phase_frame));
        used_fill<scalar_t>(m_pre_syn_phase_frame, 0,
                            sizeof(m_pre_syn_phase_frame));
        used_fill<scalar_t>(m_ana_frame, 0, sizeof(m_ana_frame));
        used_fill<scalar_t>(m_syn_frame, 0, sizeof(m_syn_frame));
        init_window();
    }
    void correct_phase(scalar_t syn_hop_size) {
        // Processing
        for (unsigned iBin = 1; iBin < (window_size / 2); iBin++) {

            scalar_t omega_bin = m_omega * iBin;
            scalar_t delta_phi = m_fft_buffer2[iBin].phase
                    - m_pre_ana_phase_frame[iBin];
            m_pre_ana_phase_frame[iBin] = m_fft_buffer2[iBin].phase;
            scalar_t frequency_deviation = wrap_phase
                    (delta_phi - m_ana_hop_size * omega_bin);
            scalar_t true_frequency = frequency_deviation * m_ana_hop_inv
                    + omega_bin;
            // calc output phase, first synth phase is 0
            scalar_t syn_phase = m_pre_syn_phase_frame[iBin]
                    + syn_hop_size * true_frequency;
            // set last
            m_pre_syn_phase_frame[iBin] = syn_phase;
            // save corrected phase
            m_fft_buffer2[iBin].phase = syn_phase;
        }
    }
    /**
     * @brief process
     * @param input_buffer
     * @param output_buffer can be the same as output
     * @param input_size maximum size of samples which is copied to
     * @ref output_buffer, if it is bigger than @ref window_size then
     * @ref window_size number of samples are copied; this
     * parameter exist only to be able to process blocks smaller than
     * @ref window_size, e.g to process a whole wav file which don't
     * necessarily has a number of samples which are a multiple of
     * the @ref window_size
     * @param scaling_factor pitch scaling factor
     */
    void process(scalar_t *input_buffer, scalar_t *output_buffer,
                 uint32_t input_size, scalar_t scaling_factor) {

        uint32_t syn_hop_size = std::round(m_ana_hop_size * scaling_factor);
        if(window_size * scaling_factor > max_syn_size) {
            // not enought memory for scaling with this high scaling factor
            return;
        }
        if(input_size == 0) {
            return;
        }
        if(input_size > window_size) {
            input_size = window_size;
        }

        // copy input to ana frame, which saves the current and the last frame
        for(uint32_t iSample = 0; iSample < input_size; iSample++,
            m_ana_frame_it++) {
            *m_ana_frame_it = input_buffer[iSample];
        }

        // we begin at the last frame + ana_hopsize, assumes wrap around
        auto ana_beginning_it = m_ana_frame_it + m_ana_hop_size;

        for(uint32_t iHop = 0; iHop < m_overlap_factor
            && (iHop * m_ana_hop_size < input_size); iHop++) {
            // get frame and window window it
            auto ana_it = ana_beginning_it + iHop * m_ana_hop_size;
            for(uint32_t iSample = 0; iSample < window_size; iSample++) {
                m_fft_buffer1[iSample] = *ana_it * m_fft_window[iSample];
                ana_it++;
            }

            m_fft_adapter.rfft(m_fft_buffer1,
                               reinterpret_cast<scalar_t*>(m_fft_buffer2));

            // convert to polar, in-place
            cmplx_cartesian_to_polar(reinterpret_cast<scalar_t*>
                                     (m_fft_buffer2),
                                     reinterpret_cast<scalar_t const *>
                                     (m_fft_buffer2), (window_size / 2));

            correct_phase(syn_hop_size);

            // covert to cartesian again
            cmplx_polar_to_cartesian(reinterpret_cast<scalar_t*>
                                     (m_fft_buffer2),
                                     reinterpret_cast<scalar_t const*>
                                     (m_fft_buffer2), (window_size / 2));

            // do inverse FFT
            m_fft_adapter.rifft(reinterpret_cast<scalar_t*>(m_fft_buffer2),
                                m_fft_buffer1);

            // overlap add
            auto syn_it = m_syn_frame_it + iHop * syn_hop_size;

            uint32_t iOverlapSample = 0;
            while(iOverlapSample < (window_size - syn_hop_size)) {
                *syn_it += m_fft_buffer1[iOverlapSample]
                        * m_fft_window[iOverlapSample];
                iOverlapSample++;
                syn_it++;
            }

            // last chunk of syn_hop_size must be overwritten as
            // this is a circular buffer and we don't want to overlap-add
            // the old data
            while(iOverlapSample < window_size) {
                *syn_it = m_fft_buffer1[iOverlapSample]
                        * m_fft_window[iOverlapSample];
                iOverlapSample++;
                syn_it++;
            }
        }

        auto copy_it = m_syn_frame_it;// + syn_hop_size
                //* (m_overlap_factor / 2 - 1);
        // on esample is shared between last and current frame
        // we have to scale by the (integer) input_length / output_lenth ratios
        // to have exact periodicity in a real time application where we don't
        // have access to the next buffer, this leads to small detune.
        // detunes for syn_hop == 114 and ana_hop == 128 under 1 cent
        // for window_size == 512 , ana_hop == 128 and factor 0,8789 it
        // detunes ~1.92 cents, so it's negligible
        resample_to_real((copy_it - 1), (m_overlap_factor * syn_hop_size) + 1,
                         output_buffer, window_size,
                         static_cast<scalar_t>(
                             static_cast<scalar_t>((m_overlap_factor
                                                    * syn_hop_size))
                             / (window_size)), input_size);

        m_syn_frame_it += syn_hop_size * m_overlap_factor;
    }

    /**
     * @brief init_window window for the window function
     */
    void init_window() {
        // gain correction
        scalar_t fac = 1. / used_sqrt(((window_size / m_ana_hop_size) / 2.));
        for (unsigned iSample = 0; iSample < window_size; iSample++) {
            m_fft_window[iSample] = hann((static_cast<scalar_t>(iSample) /
                                          window_size)) * fac;
        }
    }

    uint32_t const m_overlap_factor;
    uint32_t const m_ana_hop_size;

    scalar_t const m_omega = static_cast<scalar_t>(2 * pi) / window_size;
    scalar_t const m_ana_hop_inv;

    // the '__attribute__((aligned(16)))' (require 16 byte alignment)
    // is needed only for the FFTW adapter (you can also
    // disable SIMD in the FFTW adapter)
    scalar_t m_fft_buffer1[window_size] __attribute__((aligned(16)));
    ComplexPolar<scalar_t> m_fft_buffer2[window_size / 2 + 1]
    __attribute__((aligned(16)));

    /**
     * @brief m_pre_ana_phase_frame previous analysis frame phase
     */
    scalar_t m_pre_ana_phase_frame[window_size / 2]
    __attribute__((aligned(16)));
    /**
     * @brief m_pre_syn_phase_frame previous synthesis frame phase
     */
    scalar_t m_pre_syn_phase_frame[window_size / 2]
    __attribute__((aligned(16)));
    /**
     * @brief m_fft_window precomputed values of the window function
     * for all samples of the window_size
     */
    scalar_t m_fft_window[window_size] __attribute__((aligned(16)));
    scalar_t m_ana_frame[window_size * 2] __attribute__((aligned(16)));
    circular_iterator<scalar_t> m_ana_frame_it
    { m_ana_frame, window_size * 2 };
    scalar_t m_syn_frame[max_syn_size * 2] __attribute__((aligned(16)));

    circular_iterator<scalar_t> m_syn_frame_it
    { m_syn_frame, max_syn_size * 2 };

    fft_adapter_t &m_fft_adapter;
};

} // namespace pv

#endif // PHASEVOCODER_HPP
