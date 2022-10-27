#include "../include/vssdsp.h"


namespace poddsp {

    /// ToDo сделать нормальную квадратурную составляющую сигнала
    std::vector<std::complex<float>> complexMagModulator(const std::vector<std::complex<float>> & carrier_buffer,
                                                         const std::vector<float> & info_signal,
                                                         const float & modulation_depth) {
        if(fabsf(modulation_depth) > 1.0f){
            throw std::invalid_argument(ERROR_MAG_MODULATOR"modulation depth is too far from [0; 1] interval");
        }
        std::vector<std::complex<float>> res_arr;
        std::vector<float> intrm_arr;

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
        for(int i = 0; i < carrier_buffer.size(); i++) {
            mag = complexVectorMagnitude(carrier_buffer[i]);
            mag *= complexVectorMagnitude(information_signal[i]);
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
            res_arr.emplace_back(complexVectorMagnitude(e));
        }
        res_arr = signalShelf(res_arr, -signalMedValue(res_arr));
        return res_arr;
    }

    std::vector<std::complex<float>> complexBPSKModulator(const std::vector<bool> &info_seq, const int & samples_per_symbol){
        std::vector<float> inter_res_arr;
        std::vector<std::complex<float>> res_arr;

        int periods_counter = static_cast<int>(info_seq.size());

        std::vector<float> carrier = projection::takeProjection(
                complexSin(static_cast<float>(periods_counter),
                           periods_counter*samples_per_symbol));

        auto modulation_step = static_cast<float>(carrier.size()) / static_cast<float>(periods_counter);

        for(int i = 0; i < periods_counter; i++){
            if(!info_seq[i]){
                for (int j = 0; j < static_cast<int>(modulation_step); j++) {
                    if (i * static_cast<int>(modulation_step) + j > carrier.size()) break;
                    inter_res_arr.emplace_back(-carrier[i * static_cast<int>(modulation_step) + j]);
                }
            } else{
                for (int j = 0; j < static_cast<int>(modulation_step); j++) {
                    if (i * static_cast<int>(modulation_step) + j > carrier.size()) break;
                    inter_res_arr.emplace_back(carrier[i * static_cast<int>(modulation_step) + j]);
                }
            }
            for(auto arr = quadro_cast(inter_res_arr); auto e : arr){
                res_arr.emplace_back(e);
            }
            inter_res_arr = std::vector<float>();
        }
//        res_arr = quadro_cast(inter_res_arr);
        return res_arr;
    }
}
