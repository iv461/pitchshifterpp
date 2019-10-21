### Real Time Phase Vocoder

Real-time capable implementation of the phase vocoder algorithm for audio pitch shifting.

### Building

#### On Linux

Install dependencies:
``` sudo apt install libfftw3-dev && sudo apt install libsndfile1-dev ```

Build:
``` mkdir build && cd build && cmake .. && make ```

Launch:

``` ./pitchshifterpp ../sounds/sine1k.wav ```

A file named ``` output.wav ``` should get generated.

#### On Windows

For windows the dependencies are included, build in the same way like for Linux.


### License

Licensed under the GPLv3 License.

### Third Party Licenses

[FFTW](http://www.fftw.org/) is licensed under the GNU General Public License.

[libsndfile](http://www.mega-nerd.com/libsndfile) is licensed under the GNU Lesser General Public License, either version 2.1 or version 3.


### References

[TRADITIONAL (?) IMPLEMENTATIONS OF A PHASE-VOCODER : THE TRICKS OF THE TRADE](https://pdfs.semanticscholar.org/f1ec/2695adfb65c439d75837b342d6d7b3cc642a.pdf)
