# Amplitude modulator with Chebyshev polynomials for MaxMSP

This repository contains source codes for an external object of MaxMSP that applies self amplitude modulation in higher-order terms to the input signal. This object utilises Chebyshev polynomials and provides pure harmonics for sine waves. For other waveforms, some complex harmonics will be applied due to their cross terms in the spectrum structure.

## Features

- **Higher Term Self Amplitude Modulation:** Applies amplitude modulation up to higher-order power terms.
- **Additive Synthesis:** Provides pure harmonics for sine waves as an input, effectively functioning as additive synthesis.
- **Harmonics design:** Harmonic structure can be specified via list message or `mc.*` signal.

## Compatibility

- The object is optimized for **macOS** and does not work on other operating systems.

## Disclaimer

The developers take no responsibility for any issues or damages that may occur from using this object.

## License
This project is licensed under the [Unlicense](https://unlicense.org/).

## Installation
- Requires [MaxSDK](https://github.com/Cycling74/max-sdk).
- Clone this repository within the MaxSDK directory and use CMake as with other objects.

## Usage
See `example.maxpat` for usage instructions.

Feel free to provide feedback to the project via GitHub.
