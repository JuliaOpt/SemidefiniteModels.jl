# SemidefiniteModel

[![Build Status](https://travis-ci.org/blegat/SemidefiniteModel.jl.svg?branch=master)](https://travis-ci.org/blegat/SemidefiniteModel.jl)
[![Coverage Status](https://coveralls.io/repos/blegat/SemidefiniteModel.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/blegat/SemidefiniteModel.jl?branch=master)
[![codecov.io](http://codecov.io/github/blegat/SemidefiniteModel.jl/coverage.svg?branch=master)](http://codecov.io/github/blegat/SemidefiniteModel.jl?branch=master)

This package extends MathProgBase with `SDModel` representing a semidefinite programming problem in the following form
```
max ⟨C, X⟩            min ⟨b, y⟩
    ⟨A_i, X⟩ = b_i        ∑ A_i y_i ⪰ C
          X  ⪰ 0
```
