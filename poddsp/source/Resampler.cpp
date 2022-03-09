#include "../include/poddsp.h"


namespace PodDSP {

    /// Изменяет количество комплексных семплов сигнала.
    /// Парамертры:
    /// 1) Входной массив семплов,
    /// 2) Массив для выходных данных,
    /// 3) Количество семплов в новом массиве.
    void complexSignalResampler(const std::vector<std::complex<float>> &input_arr,
                                        std::vector<std::complex<float>> &output_arr,
                                        const int &new_samples_count)
                                        noexcept {
        unsigned int h_len = 1;
        float r = static_cast<float>(new_samples_count) / static_cast<float>(input_arr.size());
        float bw = 0.0001f;
        float slsl = 1.0f;
        unsigned int npfb = 32;

        resamp_crcf q = resamp_crcf_create(r, h_len, bw, slsl, npfb);

        auto n = (unsigned int) ceilf(r);
        std::complex<float> output_buff[n];
        unsigned int num_written;
        for (auto e: input_arr) {
            resamp_crcf_execute(q, e, output_buff, &num_written);
            for (int i = 0; i < num_written; i++) {
                output_arr.emplace_back(output_buff[i]);
            }
        }
        resamp_crcf_destroy(q);
    }
}
