#include "../include/poddsp.h"


namespace PodDSP {

    std::vector<std::complex<float>> complexPLL(const std::vector<std::complex<float>> &incoming_samples,
                                                int region_width_for_differentiation,
                                                int count_of_processing) {

        if (region_width_for_differentiation == 1) {
            region_width_for_differentiation = static_cast<int>(incoming_samples.size());
        }

        if ((incoming_samples.size() / region_width_for_differentiation * region_width_for_differentiation) !=
            incoming_samples.size()) {
            throw std::invalid_argument(ERROR_PLL "Can't reach integer count of iterations with current region "
                                        "for differentiation width and incoming array size");
        }

        std::vector<std::complex<float>> work_arr;
        std::vector<float> phase_dependence;
        std::vector<float> diff_phase_dependence;
        std::vector<float> second_diff_phase_dependence;
        std::vector<std::complex<float>> res_arr;

        for (int j = 0; j < static_cast<int>(incoming_samples.size() / region_width_for_differentiation); j++) {
            for (int i = 0; i < region_width_for_differentiation; i++) {
                work_arr.emplace_back(incoming_samples[i]);
            }
//            for (auto e: work_arr) {
//                phase_dependence.emplace_back(complexVectorPhase(e));
//            }

            phase_dependence = complexSignalPhaseDependence(work_arr);

//            PlotConstructor::drawPlot(phase_dependence, "Фаза сигнала");

            diff_phase_dependence = differentiation(phase_dependence);

//            PlotConstructor::drawPlot(diff_phase_dependence, "Частота сигнала");

            second_diff_phase_dependence = differentiation(diff_phase_dependence);

//            PlotConstructor::drawPlot(second_diff_phase_dependence, "Производная частоты сигнала");

            work_arr = complexCFOCompensator(work_arr, signalMedValue(second_diff_phase_dependence));

            for (auto e: work_arr) {
                res_arr.emplace_back(e);
            }

            work_arr = std::vector<std::complex<float>>();
//            phase_dependence = std::vector<float>();
        }
        if (count_of_processing == 1) {
            return res_arr;
        } else {
            return complexPLL(res_arr, region_width_for_differentiation, count_of_processing - 1);
        }
    }

    std::vector<std::complex<float>> complexCFOCompensator(const std::vector<std::complex<float>> &incoming_arr,
                                                           const float &phase_attenuation_per_sample_rad) noexcept {

        return complexPhaseChanger(incoming_arr,0.0f, -static_cast<float>(phase_attenuation_per_sample_rad TO_DEG));
    }

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
            if (res_dependence.size() > 1 && fabsf(curr_diff) > res_dependence.back()*100) curr_diff = res_dependence.back();
            res_dependence.emplace_back(curr_diff);
        }
        return res_dependence;
    }
}
