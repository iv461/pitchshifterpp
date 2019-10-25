// std
#include <vector>
#include <iostream>
#include <chrono>
#include <memory>
#include <fstream>

// C std lib
#include <stdio.h>
#include <stdlib.h>

// libsnd
#include <libsndfile/sndfile.hh>

#include <pitchshifterpp/fftwpvadapter.hpp>
#include <pitchshifterpp/phasevocoder.hpp>

using std::cout;
using std::string;
using std::vector;
using namespace std::chrono;

int main(int argc, char* argv[])
{
    if(argc != 3) {
        cout << "Invalid number of arguments provided\n";
        return -1;
    }
    string input_file_str(argv[1]);
    SndfileHandle in_snd_handle(input_file_str, SFM_READ);
    if(!in_snd_handle) {
        cout << "Could not open input sound file";
        return -1;
    }
    if(in_snd_handle.channels() != 1) {
        cout << "no handling for non-mono files implemented, exiting";
        return -1;
    }
    using scalar_t = double;
    constexpr uint32_t window_size = 512;
    constexpr scalar_t max_scale_factor = 1.5;
    constexpr scalar_t overlap_factor = 4;

    int semitones_to_scale;
    try {
        semitones_to_scale = std::stoi(string(argv[2]));
    } catch(...){
        cout << "could not parse scaling factor\n";
        return -1;
    }
    scalar_t scaling_factor = pv::semitone_to_scale_factor(semitones_to_scale);
    if(scaling_factor > max_scale_factor) {
        cout << "scaling factor too big\n";
        return -1;
    }
    auto input_size = in_snd_handle.frames()  * in_snd_handle.channels();
    vector<scalar_t> input_samples;
    input_samples.resize(static_cast<size_t>(input_size));
    in_snd_handle.readf(&input_samples[0], in_snd_handle.frames());


    FFTWPVAdapter fftw_adapter(window_size);
    using pv_t = pv::PhaseVocoder<window_size,
    static_cast<uint32_t>(window_size * max_scale_factor),
    FFTWPVAdapter, scalar_t>;
    auto phase_vocoder = std::make_unique<pv_t>(fftw_adapter, overlap_factor);

    vector<scalar_t> output_samples(input_samples.size());

    auto start_clock = steady_clock::now();

    auto samples_remaining =  input_samples.size();
    uint32_t iFrames = 0;
    while(true) {
        phase_vocoder->process(&input_samples[iFrames * window_size],
                &output_samples[iFrames * window_size],
                static_cast<uint32_t>(samples_remaining), scaling_factor);
        iFrames++;
        if(samples_remaining > window_size) {
            samples_remaining -= window_size;
        } else {
            break;
        }
    }

    auto end_clock = steady_clock::now();
    auto clock_ms_count = duration_cast<milliseconds>
            (end_clock - start_clock).count();

    cout << "processing of " << input_samples.size() / window_size <<
            " frames took " << clock_ms_count << "ms\n";

    SndfileHandle out_snd_handle("output.wav", SFM_WRITE, in_snd_handle.format(),
                                 in_snd_handle.channels(),
                                 in_snd_handle.samplerate());
    if(!out_snd_handle) {
        cout << "could not open file for writing\n";
        return -1;
    }

    out_snd_handle.writef(&output_samples[0], input_size);

    return 1;
}
