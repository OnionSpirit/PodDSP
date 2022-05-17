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

    float complexSequenceCorrelation(const std::vector<std::complex<float>> &original_sequence,
                                     const std::vector<std::complex<float>> &incoming_sequence) noexcept {


        int sequence_size = (int)original_sequence.size();
        std::complex<float> correlation_result;
        std::vector<std::complex<float>> conjugated_incoming_sequence;

        conjugated_incoming_sequence.reserve(incoming_sequence.size());
        for (std::complex<float> e : incoming_sequence) {
            conjugated_incoming_sequence.emplace_back(conj(e));
        }

        std::vector<std::complex<float>> centralized_sequence_original =
                sequenceCentralizer(original_sequence, sequence_size);
        std::vector<std::complex<float>> centralized_sequence_incoming =
                sequenceCentralizer(conjugated_incoming_sequence, sequence_size);

        std::complex<float> dispersion_original =
                dispersion(centralized_sequence_original, sequence_size);
        std::complex<float> dispersion_incoming =
                dispersion(centralized_sequence_incoming, sequence_size);

        for (int i = 0; i < sequence_size; i++) {
            correlation_result +=
                    centralized_sequence_incoming[i] * centralized_sequence_original[i];
        }
        correlation_result *= (1.0f / static_cast<float>(sequence_size));
        correlation_result /= pow(dispersion_incoming * dispersion_original, 0.5);

        return std::norm(correlation_result);
    }
}
