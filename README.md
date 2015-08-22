# SynthesisFilters

[![Build Status](https://travis-ci.org/r9y9/SynthesisFilters.jl.svg?branch=master)](https://travis-ci.org/r9y9/SynthesisFilters.jl)
[![Build status](https://ci.appveyor.com/api/projects/status/lox04xqpp3qo9646/branch/master?svg=true)](https://ci.appveyor.com/project/r9y9/synthesisfilters-jl/branch/master)
[![Coverage Status](https://coveralls.io/repos/r9y9/SynthesisFilters.jl/badge.svg?branch=master)](https://coveralls.io/r/r9y9/SynthesisFilters.jl?branch=master)
[![License](http://img.shields.io/badge/license-MIT-brightgreen.svg?style=flat)](LICENSE.md)

SynthesisFilters.jl provides waveform generation filters for speech synthesis, especially from mel-generalized cepstrum.

Note that this package is built on top of [SPTK.jl](https://github.com/r9y9/SPTK.jl). A part of the core is re-writen in Julia language from [Speech Signal Processing Toolkit (SPTK)](http://sp-tk.sourceforge.net/).

## Features

- **LMADF**: Log magnitude approximation digital filter for synthesis from cepstrum
- **MLSADF**: Mel-log spectrum approximation digital filter for synthesis from mel-cepstrum
- **MGLSADF**: Mel generalized-log spectrum approximation digital filter for synthesis from mel-generalized cepstrum
- **AllPoleDF**: All-pole digital filter for synthesis from LPC
- **AllPoleLatticeDF**: All-pole lattice digital filter for synthesis from PARCOR
- **LSPDF**: LSP digital filter for synthesis from LSP

## Demonstration notebook

- [Introduction notebook](http://nbviewer.ipython.org/github/r9y9/SynthesisFilters.jl/blob/master/examples/Introduction%20to%20SynthesisFilters.jl.ipynb): a brief introduction to SynthesisFilters.jl. Synthesized audio examples (Japanese) are avaiable.

## Supported Platforms

- Linux
- Mac OS X
- Windows

## Installation

```julia
Pkg.clone("https://github.com/r9y9/SynthesisFilters.jl")
```

## Related packages

- [MelGeneralizedCepstrums.jl](https://github.com/r9y9/MelGeneralizedCepstrums.jl): spectral paramter estimation based on mel-generalized cepstrum analysis
- [SPTK.jl](https://github.com/r9y9/SPTK.jl): a thin wrapepr for the SPTK
- [WORLD.jl](https://github.com/r9y9/WORLD.jl): a lightweight julia wrapper for WORLD - a high-quality speech analysis, modification and synthesis system