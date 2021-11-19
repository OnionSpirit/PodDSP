#include <iostream>
#include <complex>
#include <vector>
#include <ctime>

#define RANDOM_NUMBER 1 + rand()%10
namespace Correlation{
    std::complex<float> complexSequenceCorrelation(std::vector<std::complex<float>> const &,
                                                   std::vector<std::complex<float>> const &);
    std::vector<std::complex<float>> sequenceCentralizer(std::vector<std::complex<float>> const &, int);
    std::complex<float> dispersion(std::vector<std::complex<float>> const &, int);
}

int main() {

    int sequence_length = 10;
    srand(time(nullptr));

    std::vector<std::complex<float>> original_sequence;
    std::vector<std::complex<float>> incoming_sequence;

    for(int i = 0; i < sequence_length; i++){
        original_sequence.emplace_back(RANDOM_NUMBER); original_sequence.back().imag(RANDOM_NUMBER);
        incoming_sequence.emplace_back(RANDOM_NUMBER); incoming_sequence.back().imag(RANDOM_NUMBER);
    }

    std::cout << std::norm(Correlation::complexSequenceCorrelation(original_sequence, original_sequence)) << std::endl;

    return 0;
}

namespace Correlation {
    std::complex<float> complexSequenceCorrelation(std::vector<std::complex<float>> const &original_sequence,
                                                   std::vector<std::complex<float>> const &incoming_sequence) {

        int sequence_size = (int) original_sequence.size();
        std::complex<float> correlation_result;
        std::vector<std::complex<float>> conjugated_incoming_sequence;

        for (std::complex<float> e: incoming_sequence) {
            conjugated_incoming_sequence.emplace_back(conj(e));
        }

        std::vector<std::complex<float>> centralized_sequence_original = sequenceCentralizer(original_sequence, sequence_size);
        std::vector<std::complex<float>> centralized_sequence_incoming = sequenceCentralizer(conjugated_incoming_sequence, sequence_size);

        std::complex<float> dispersion_original = dispersion(centralized_sequence_original, sequence_size);
        std::complex<float> dispersion_incoming = dispersion(centralized_sequence_incoming, sequence_size);

        for (int i = 0; i < sequence_size; i++) {
            correlation_result += centralized_sequence_incoming[i] * centralized_sequence_original[i];
        } correlation_result *= (1.0f / static_cast<float>(sequence_size)); correlation_result /= pow(dispersion_incoming * dispersion_original, 0.5);

        return correlation_result;
    }

    std::vector<std::complex<float>> sequenceCentralizer(std::vector<std::complex<float>> const &sequence, int sequence_size) {

        std::complex<float> average;
        std::vector<std::complex<float>> centralized_sequence;

        for(auto e : sequence){
            average += e;
        } average *= (1.0f / static_cast<float>(sequence_size));

        for(auto e : sequence){
            centralized_sequence.emplace_back(e - average);
        }

        return centralized_sequence;
    }

    std::complex<float> dispersion(std::vector<std::complex<float>> const &sequence, int sequence_size) {

        std::complex<float> D;

        for (std::complex<float> e: sequence) {
            D += std::norm(e);
        } D *= (1.0f / static_cast<float>(sequence_size));

        return D;
    }
}