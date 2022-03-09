#include "../include/poddsp.h"


namespace PodDSP {

    std::vector<std::complex<float>> complexSin(const float &freq,
                                                const int &count_of_samples,
                                                const float &initial_phase_deg) {

        float count_of_periods = freq;
        if ((float) count_of_samples <= 2 * freq) {
            throw std::invalid_argument(
                    ERROR_GEN"Kotelnikov theorem requires sampling freq more than doubled signal freq.");
        }

        float samples_per_period = static_cast<float>(count_of_samples) / static_cast<float>(count_of_periods);
        float sample_phase_offset = static_cast<float>(2 * M_PI) / samples_per_period;
        std::vector<std::complex<float>> generated_signal;
        std::complex<float> sample{0, 0};

        auto initial_phase_rad = static_cast<float>(2 * M_PI * initial_phase_deg / 360);

        for (int i = 0; i < count_of_samples; i++) {
            sample.real(static_cast<float>(cos(static_cast<double>(static_cast<float>(i) * sample_phase_offset + initial_phase_rad))));
            sample.imag(static_cast<float>(sin(static_cast<double>(static_cast<float>(i) * sample_phase_offset + initial_phase_rad))));
            generated_signal.emplace_back(sample);
        }
        return generated_signal;
    }

    std::vector<float> signalShelf(const std::vector<float> & signal,
                                   const float & shelf)
                                   noexcept {
        std::vector<float> res_arr;

        res_arr.reserve(signal.size());
        for(auto e : signal){
            res_arr.emplace_back(e + shelf);
        }
        return res_arr;
    }
}
