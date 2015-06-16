# SynthesisFilters

[![Build Status](https://travis-ci.org/r9y9/SynthesisFilters.jl.svg?branch=master)](https://travis-ci.org/r9y9/SynthesisFilters.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/lox04xqpp3qo9646/branch/master?svg=true)](https://ci.appveyor.com/project/r9y9/synthesisfilters-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/r9y9/SynthesisFilters.jl/badge.svg?branch=master)](https://coveralls.io/r/r9y9/SynthesisFilters.jl?branch=master)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)

SynthesisFilters.jl provides waveform generation filters for speech synthesis, especially from mel-generalized cepstrum. The core is re-coded from [Speech Signal Processing Toolkit (SPTK)](http://sp-tk.sourceforge.net/).

## Features

- All Pole Digital Filter (AllPoleDF)
- Log Magnitude Approximation Digital Filter (LMADF)
- Mel-Log Spectrum Approximation Digital Filter (MLSADF)
- Mel Generalized-Log Spectrum Approximation Digital Filter (MGLSADF)

## Compare vocoders

try [examples/compare_vocoders.jl](examples/compare_vocoders.jl)
