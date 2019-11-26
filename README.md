## Real Time Phase Vocoder

Real-time capable implementation of the phase vocoder algorithm for audio pitch shifting.
Single C++ Header implementation with additional headers for the used math functions
and the adapter for the used fft functions.

Memory relevant parameters like the window size are template parameters, all needed memory is 
reserved at compile time, thus making it suitable for memory-constrained embedded platforms.

The ``` phasevocoder.hpp ``` header is dependency-free but a FFT implementation is needed
for which ``` FFTW ``` is used and for the test application which operates
on Wave files, ```libsdnfile``` is needed.

## Building

### On Linux

Install dependencies, e.g with APT:

``` sudo apt install libfftw3-dev libsndfile1-dev ```

Build:

``` mkdir build && cd build && cmake .. && make ```

Launch:

``` ./pitchshifterpp <wav_file> <semitones_to_scale> ```

With example sound files:

``` ./pitchshifterpp sine1k.wav <semitones_to_scale> ```

The file ``` output.wav ``` should get generated.

### On Windows

For windows the dependencies (x64) are included, build in the same way like for Linux.
Tested with MinGW and MSVC.

## License

Licensed under the GNU General Public License, version 3.

## Third Party Licenses

[FFTW](http://www.fftw.org/) is licensed under the GNU General Public License.

[libsndfile](http://www.mega-nerd.com/libsndfile) is licensed under the GNU Lesser General Public License, either version 2.1 or version 3.


## References

[TRADITIONAL (?) IMPLEMENTATIONS OF A PHASE-VOCODER : THE TRICKS OF THE TRADE](https://pdfs.semanticscholar.org/f1ec/2695adfb65c439d75837b342d6d7b3cc642a.pdf)
