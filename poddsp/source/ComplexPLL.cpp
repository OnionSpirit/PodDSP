#include "../include/poddsp.h"


namespace poddsp {

/// ToDo Добавить Точность для сглаживания графиков частоты и фазы, для эффективной работы гетеродина
    std::vector<std::complex<float>> complexPLL(const float spec_freq, const std::vector<std::complex<float>> &incoming_samples,
                                                int count_of_processing) {

        std::vector<float> phase_dependence;
        std::vector<float> diff_phase_dependence;
        std::vector<float> second_diff_phase_dependence;
        std::vector<std::complex<float>> res_arr;

        phase_dependence = complexSignalPhaseDependence(incoming_samples);
//        PlotConstructor::drawPlot(phase_dependence, "фаза в сигнале");

        diff_phase_dependence = differentiation(phase_dependence);

//        PlotConstructor::drawPlot(diff_phase_dependence, "частота в сигнале");

        second_diff_phase_dependence = differentiation(diff_phase_dependence);
//        PlotConstructor::drawPlot(second_diff_phase_dependence, "рост частоты в сигнале");

        res_arr = complexCFOCompensator(incoming_samples, signalMedValue(second_diff_phase_dependence));

        phase_dependence = complexSignalPhaseDependence(incoming_samples);
        diff_phase_dependence = differentiation(phase_dependence);

        {
            auto freq = signalMedValue(diff_phase_dependence);
            PlotConstructor::drawPlot(diff_phase_dependence, "частота в сигнале");
            if ( freq > spec_freq ) {
                res_arr = Heterodyne( freq - spec_freq, res_arr, false);
            } else if (freq < spec_freq) {
                res_arr = Heterodyne( spec_freq - freq, res_arr);
            }
            std::cout << freq << std::endl;
        }

        return res_arr;

//        if (count_of_processing == 1) {
//            return res_arr;
//        } else {
//            return complexPLL(spec_freq, res_arr, count_of_processing - 1);
//        }
    }

    std::vector<std::complex<float>> complexCFOCompensator(const std::vector<std::complex<float>> &incoming_arr,
                                                           const float &phase_attenuation_per_sample_rad) noexcept {

        return complexPhaseChanger(incoming_arr,  0.0f, -static_cast<float>(phase_attenuation_per_sample_rad TO_DEG));
    }
/// Todo Просмотреть
    std::vector<float> differentiation(const std::vector<float> &incoming_dependence){

        if(incoming_dependence.size() < 2){
            throw std::invalid_argument("Too short array");
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
