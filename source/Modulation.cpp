#include "../include/poddsp.h"


namespace poddsp {

    std::vector<std::complex<float>> complexMagModulator(const std::vector<std::complex<float>> & carrier_buffer,
                                                         const std::vector<float> & info_signal,
                                                         const float & modulation_depth) {
        if(fabsf(modulation_depth) > 1.0f){
            throw std::invalid_argument(ERROR_MAG_MODULATOR"modulation depth is too far from [0; 1] interval");
        }
        std::vector<std::complex<float>> res_arr;

        std::vector<std::complex<float>> information_signal_complex;
        std::vector<float> information_signal;
        std::vector<std::complex<float>> info_signal_complex;
        info_signal_complex.reserve(info_signal.size());
        for(auto e : info_signal){
            info_signal_complex.emplace_back(std::complex<float>{e, 0});
        }
        complexSignalResampler(info_signal_complex, information_signal_complex, static_cast<int>(carrier_buffer.size()));
        information_signal.reserve(information_signal_complex.size());
        for(auto e : information_signal_complex){
            information_signal.emplace_back(e.real());
        }

        if(modulation_depth == 0) return carrier_buffer;
        information_signal = signalShelf(information_signal, signalMaxValue(information_signal) / modulation_depth);

        float mag;
        float current_phase;
        for(int i = 0; i < carrier_buffer.size(); i++){
            mag = complexMagMeasurer(carrier_buffer[i]);
            mag *= complexMagMeasurer(information_signal[i]);
            current_phase = complexVectorPhase(carrier_buffer[i]);
            current_phase += complexVectorPhase(information_signal[i]);
            res_arr.emplace_back(std::complex<float>{
                    mag * static_cast<float>(cos(static_cast<double>(current_phase))),
                    mag * static_cast<float>(sin(static_cast<double>(current_phase)))});
        }
        return res_arr;
    }

    std::vector<float> complexMagDemodulator(const std::vector<std::complex<float>> &modulated_signal) noexcept{

        std::vector<float> res_arr;
        res_arr.reserve(modulated_signal.size());
        for(auto e : modulated_signal){
            res_arr.emplace_back(complexMagMeasurer(e));
        }
        res_arr = signalShelf(res_arr, -signalMedValue(res_arr));
        return res_arr;
    }

    std::vector<std::complex<float>> complexBPSKModulator(const std::vector<std::complex<float>> &carrier,
                                                          const std::vector<bool> &info_seq){
        std::vector<float> inter_res_arr;
        std::vector<std::complex<float>> res_arr;
        int periods_counter = 0;

        {
            std::vector<float> phase_dependence = complexSignalPhaseDependence(carrier);
            float curr_diff;

            for (int i = 0; i < phase_dependence.size(); i++) {
                if (i + 1 == phase_dependence.size()) {
                    break;
                }
                curr_diff = ((phase_dependence[i + 1] - phase_dependence[i]));
                if (fabsf(curr_diff) > 3+M_PI/2) periods_counter++;
            }
        }

        if(periods_counter != info_seq.size() - 1){
            throw std::invalid_argument(ERROR_BPSK_MODUlATOR "Modulation sequence must be equal to count of periods in carrier frame");
        }

        auto modulation_step = static_cast<float>(carrier.size()) / static_cast<float>(periods_counter);


        for(int i = 0; i < info_seq.size(); i++){
            if(!info_seq[i]){
                for(int j = 0; j < static_cast<int>(modulation_step); j++){
                    inter_res_arr.emplace_back(-carrier[i * static_cast<int>(modulation_step) + j].real());
//                        res_arr.emplace_back(std::complex<float>{-carrier[i * static_cast<int>(modulation_step) + j].real(),
//                                                                 -carrier[i * static_cast<int>(modulation_step) + j - static_cast<int>(modulation_step)/4].imag()});
                }
            } else {
                for(int j = 0; j < static_cast<int>(modulation_step); j++){
                    inter_res_arr.emplace_back(carrier[i * static_cast<int>(modulation_step) + j].real());
//                    res_arr.emplace_back(std::complex<float>{carrier[i * static_cast<int>(modulation_step) + j].real(),
//                                                             carrier[i * static_cast<int>(modulation_step) + j - static_cast<int>(modulation_step)/4].imag()});
                }
            }
        }

        for(int i = 0; i < info_seq.size() - 2; i++){
            for(int j = 0; j < static_cast<int>(modulation_step); j++){
                if(i - 1 < 0) continue;
                res_arr.emplace_back(std::complex<float>{inter_res_arr[i * static_cast<int>(modulation_step) + j],
                                                         inter_res_arr[i * static_cast<int>(modulation_step) + j + static_cast<int>(modulation_step)/4]});
            }
        }

        return res_arr;
    }
}
