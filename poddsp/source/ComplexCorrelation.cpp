#include "../include/poddsp.h"


namespace poddsp {

    std::vector<std::complex<float>> sequenceCentralizer(const std::vector<std::complex<float>> &sequence,
                                                         int sequence_size)
    noexcept {

        std::complex<float> average;
        std::vector<std::complex<float>> centralized_sequence;

        for (auto e: sequence) {
            average += e;
        }
        average *= (1.0f / static_cast<float>(sequence_size));

        centralized_sequence.reserve(sequence.size());
        for (auto e: sequence) {
            centralized_sequence.emplace_back(e - average);
        }

        return centralized_sequence;
    }

    std::complex<float> dispersion(const std::vector<std::complex<float>> &sequence,
                                   int sequence_size)
    noexcept {

        std::complex<float> D;

        for (std::complex<float> e: sequence) {
            D += std::norm(e);
        }
        D *= (1.0f / static_cast<float>(sequence_size));

        return D;
    }
}
