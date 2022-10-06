#pragma once

namespace poddsp{

    std::vector<std::complex<float>> sequenceCentralizer(const std::vector<std::complex<float>> &,
                                                         int)
    noexcept;

    std::complex<float> dispersion(std::vector<std::complex<float>> const &,
                                   int)
    noexcept;
    /// Коррелятор двух последовательностей семплов комплексных сигналов, так же работает с некомплексными сигналами.

    float complexSequenceCorrelation(const std::vector<std::complex<float>> &,
                                     const std::vector<std::complex<float>> &) noexcept;

}
