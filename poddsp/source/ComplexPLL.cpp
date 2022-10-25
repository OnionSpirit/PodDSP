#include "../include/poddsp.h"


namespace poddsp {

/// ToDo Добавить Точность для сглаживания графиков частоты и фазы, для эффективной работы гетеродина
    std::vector<std::complex<float>> complexPLL(const std::vector<std::complex<float>> &incoming_samples,
                                                float spec_freq) noexcept {

        std::vector<float> phase_dependence;
        std::vector<float> diff_phase_dependence;
        std::vector<float> second_diff_phase_dependence;
        std::vector<std::complex<float>> res_arr;

        phase_dependence = complexSignalPhaseDependence(incoming_samples);
        phase_dependence = phaseDependenceLining(phase_dependence);
        diff_phase_dependence = smoother(differentiation(phase_dependence));
        second_diff_phase_dependence = differentiation(diff_phase_dependence);
        res_arr = complexCFOCompensator(incoming_samples, signalMedValue(second_diff_phase_dependence));

        {
            bool freq_move_type = false;
            phase_dependence = phaseDependenceLining(complexSignalPhaseDependence(res_arr));
            phase_dependence = phaseModulationSkipEraser(phase_dependence);
            auto freq = fabs(phase_dependence.back()TO_DEG/360);

            if (spec_freq != freq and spec_freq != 0) {
                if (spec_freq > freq) freq_move_type = true;
                int freq_diff = fabs(spec_freq - freq);
                res_arr = Heterodyne((float)freq_diff, res_arr, freq_move_type);
            }
        }

        return res_arr;
    }

    std::vector<std::complex<float>> complexCFOCompensator(const std::vector<std::complex<float>> &incoming_arr,
                                                           const float &phase_attenuation_per_sample_rad) noexcept {

        return complexPhaseChanger(incoming_arr,  0.0f, -static_cast<float>(phase_attenuation_per_sample_rad TO_DEG));
    }
/// Todo Просмотреть
    std::vector<float> differentiation(const std::vector<float> &incoming_dependence) noexcept {

        if(incoming_dependence.size() < 2){
            return incoming_dependence;
        }
        std::vector<float> res_dependence;
        float curr_diff = 0.0f;
        for (int i = 0; i < incoming_dependence.size(); i++) {
            if (i + 1 == incoming_dependence.size()) {
                curr_diff = res_dependence.back();
                res_dependence.emplace_back(curr_diff);
                break;
            }
            curr_diff = ((incoming_dependence[i + 1] - incoming_dependence[i]));
//            if (res_dependence.size() > 1 && fabsf(curr_diff) > res_dependence.back()*100) curr_diff = res_dependence.back();
            res_dependence.emplace_back(curr_diff);
        }
        return res_dependence;
    }
}
