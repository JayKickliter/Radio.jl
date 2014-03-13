# Radio.jl

A digital communications package for the Julia language. 

## Status
This package is brand new as of 9 March 2014. It has almost no functionality and should not be used by anyone. That said, if you have experience in DSP, digital communications, or RF test and measurment, please feel free to contribute.

## Proposed Structure
This is a growing list of proposed functionality and package strcture.

* Modulation
	* **PSK.jl**: Phase Shift Keying modulation/demodulation
	* **APSK.jl**: Amplitude Phase Shift Keying modulation/demodulation
	* **QAM.jl**: Quadrature Amplitude Modulation/demodulation 	
* Random
	* **WGN.jl**: White Gaussian Noise
* Math
  * **CZT.jl**: Chirp-z Transform
  * **FFT.jl**: ?. Need a non-GPL FFT. Possibly a native Julia implemantation or an interface to [FFTS](https://github.com/anthonix/ffts)
* Filter
	* **FIR.jl**: Fir filter design and execution
	* **Polyphase.jl**: Polyphase filter and execution
	* **Resample.jl**: Decimation, interpolation, and rational resampling. Maybe cubic interpolation.
* Support
  * **Types.jl**: IQ