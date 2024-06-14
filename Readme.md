# Digital Signal Processing

## General Info
This project was developed with the aim of creating and testing a digital demodulator before implementing it on an embedded system.
This project, programmed in C++, offers the possibility of synthesizing real and complex signals and carrying out operations between them.

In addition, in order to carry out demodulation, an initial implementation of a digital impulse response filter has been put in place. This filter offers the possibility of choosing its template (high pass, low pass, band pass and band stop), its type allowing its coefficients to be calculated (Butterworth, Bessel, Chebyshev, elliptic) and its cutoff frequency.

For the moment only the Butterworth filter is implemented, which is sufficient to carry out amplitude and phase demodulation.
Lately, in order to visualize the signals the project allows complex signals to be exported in csv format and displayed in a python script with matplotlib and numpy.

## Technologies
This project was developed under Linux with the vscode editor.

The code was completed with version g++ 11.4.0.

The Python version used is 3.10.12.