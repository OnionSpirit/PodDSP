#pragma once

namespace poddsp{
    /// Коррелятор двух последовательностей семплов комплексных сигналов.
    std::complex<float> complexSequenceCorrelation(const std::vector<std::complex<float>> &,
                                                   const std::vector<std::complex<float>> &)
                                                   noexcept;

    std::vector<std::complex<float>> sequenceCentralizer(const std::vector<std::complex<float>> &,
                                                         int)
                                                         noexcept;

    std::complex<float> dispersion(std::vector<std::complex<float>> const &,
                                   int)
                                   noexcept;
}
