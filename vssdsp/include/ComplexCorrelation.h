#pragma once

namespace poddsp{

    std::vector<std::complex<float>> sequenceCentralizer(const std::vector<std::complex<float>> &) noexcept;

    std::complex<float> dispersion(std::vector<std::complex<float>> const &) noexcept;
    float dispersion(s_sig_t const &) noexcept;
    /// Коррелятор двух последовательностей семплов комплексных сигналов, так же работает с некомплексными сигналами.

    float complexSequenceCorrelation(const std::vector<std::complex<float>> &,
                                     const std::vector<std::complex<float>> &) noexcept;

}
