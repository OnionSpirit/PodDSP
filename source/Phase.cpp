#include "../include/poddsp.h"


namespace poddsp {

    std::complex<float> complexIntermediatePhaseCalculating(const std::vector<std::complex<float>> &samples_set) {
        int n = samples_set.size();
        if (n % 2 != 0) { throw std::invalid_argument(ERROR_PHASE "Samples set must be even"); }

        int k = n / 2 - 1;
        std::complex<float> res = 1.0f;

        for (int i = 0; auto j: samples_set) {
            if (i == k) { continue; }
            res *= (pow(j, (i + 1) % 2) * pow(conj(j), i % 2));
            i++;
        }
        return res;
    }


    float complexPhaseCalculating(const std::vector<std::complex<float>> &samples_S,
                                  const std::vector<std::complex<float>> &samples_X)
    noexcept {
        std::complex<float> phase = 0.0f;
        float deg_phase =0.0f;

        try {
            phase = J *
                    log(complexIntermediatePhaseCalculating(samples_S) /
                        complexIntermediatePhaseCalculating(samples_X));
            deg_phase = static_cast<float>(phase.real() * 180 / M_PI);
        }

        catch (std::exception& e) {
            e.what();
        }

        return deg_phase;
    }

    std::vector<std::complex<float>> complexPhaseChanger(const std::vector<std::complex<float>> &incoming_arr,
                                                         const float &additional_phase_deg,
                                                         const float &phase_attenuation_per_sample_deg)
    noexcept {

        std::vector<std::complex<float>> res_arr;
        auto additional_phase_rad = (float) (additional_phase_deg * M_PI / 180);
        auto phase_time_offset_rad = (float) (phase_attenuation_per_sample_deg TO_RAD);
        auto unique_sample_phase_offset = 0.0f;
        float arg;
        float current_phase;
        for (float i = 0.0f; auto e: incoming_arr) {
            current_phase = complexVectorPhase(e);
            if (current_phase < 0) current_phase += (2 * M_PI);
            arg = complexMagMeasurer(e);
            unique_sample_phase_offset += i * phase_time_offset_rad;
            while (fabs(unique_sample_phase_offset) > 2 * M_PI) {
                if (phase_attenuation_per_sample_deg > 0.0f) unique_sample_phase_offset -= 2 * M_PI;
                if (phase_attenuation_per_sample_deg < 0.0f) unique_sample_phase_offset += 2 * M_PI;
            }
            res_arr.emplace_back(std::complex<float>{
                    static_cast<float>(arg * cos(static_cast<double>(
                                                         current_phase
                                                         + additional_phase_rad
                                                         + unique_sample_phase_offset))),
                    static_cast<float>(arg * sin(static_cast<double>(
                                                         current_phase
                                                         + additional_phase_rad
                                                         + unique_sample_phase_offset)))});
            i++;
        }

        return res_arr;
    }

    float complexVectorPhase(const std::complex<float>& sample) noexcept {
        auto phase = static_cast<float>(atan2(static_cast<double>(sample.imag()), static_cast<double>(sample.real())));
        if(phase < 0)  phase += 2 * M_PI;
        return phase;
    }

    std::vector<float> complexSignalPhaseDependence(const std::vector<std::complex<float>> &signal){

        std::vector<float> res_arr;
        for (auto e: signal) {
            res_arr.emplace_back(complexVectorPhase(e));
        }
        return res_arr;
    }

    std::vector<float> phaseDependenceLining(const std::vector<float>& phase_graph) noexcept {

        auto len = phase_graph.size();
        float boost_up = 0;
        std::vector<float> t_s; t_s.resize(len);
        for(int i = 0; i < len; i++) {

            if (i > 0 and phase_graph[i] * 2 < phase_graph[i-1]) {

                boost_up += 2 * M_PI;
            }
            t_s[i] = phase_graph[i] + boost_up;
        }

        return t_s;
    }
}
