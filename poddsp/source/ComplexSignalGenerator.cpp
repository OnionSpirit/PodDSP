#include "../include/poddsp.h"


namespace poddsp {

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

    simpleSignal impulseGen(const int& width, const int& delay, const int& frame){

        if((delay + width) > frame)
            throw std::invalid_argument(ERROR_GEN "Impulse width and delay, bigger than frame");
        simpleSignal impulse;

        impulse.reserve(frame);
        for(int i = 0; i < delay; i++){
            impulse.emplace_back(0.0f);
        }
        for(int i = 0; i < width; i++){
            impulse.emplace_back(1.0f);
        }
        for(int i = delay+width; i < frame; i++){
            impulse.emplace_back(0.0f);
        }

        return impulse;
    }

    simpleSignal MeanderGen(const float& freq,
                             const int & count_of_samples,
                             const float & zero_phase,
                             bool is_simple){

        if ((float) count_of_samples <= 2 * freq) {
            throw std::invalid_argument(
                    ERROR_GEN"Kotelnikov theorem requires sampling freq more than doubled signal freq.");
        }

        simpleSignal res_arr;
        res_arr.reserve(count_of_samples);

        float samples_per_period = static_cast<float>(count_of_samples) / freq;

        simpleSignal example = impulseGen(static_cast<int>(samples_per_period) / 2, 0,
                                          static_cast<int>(samples_per_period) / 2);

        float arg = 1;
        for (int i = 0; i < 2 * static_cast<int>(freq); i++) {
            for (auto e: example) {
                res_arr.emplace_back(arg*e);
            }
            arg *= -1.0f;
        }


        return res_arr;
    }
}
