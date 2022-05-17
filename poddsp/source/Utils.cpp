#include "../include/poddsp.h"


namespace poddsp {

    float complexMagMeasurer(const std::complex<float> &sample) noexcept {
        float res;
        res = sqrt(pow(sample.real(), 2) + pow(sample.imag(), 2));
        return res;
    }

    float signalMaxValue(const std::vector<float> &arr) noexcept {
        float res = arr[0];
        for (auto e: arr) {
            if (res <= e) {
                res = e;
            }
        }
        return res;
    }

    float signalMinValue(const std::vector<float> &arr) noexcept {
        float res = arr[0];
        for (auto e: arr) {
            if (res >= e) {
                res = e;
            }
        }
        return res;
    }

    float signalMedValue(const std::vector<float> &arr) noexcept {
        float res = 0.0f;
        for (auto e: arr) {
            res += e;
        }
        res /= static_cast<float>(arr.size());
        return res;
    }

    std::vector<float> forwardFFT(const std::vector<float> &an_seq) noexcept {

        auto length = static_cast<int>(an_seq.size());
        std::vector<std::complex<float>> intrm_arr;
        intrm_arr.resize(length);

        auto forward_FFT =  fftwf_plan_dft_1d(length, (fftwf_complex *)(intrm_arr.data()),
                                             (fftwf_complex *)(intrm_arr.data()), FFTW_FORWARD, FFTW_ESTIMATE);

        for(int i = 0; i < length; i++){
            intrm_arr[i] = -std::complex<float>(an_seq[i], 0);
        }

        fftwf_execute(forward_FFT);
        fftwf_destroy_plan(forward_FFT);
        fftwf_cleanup();

        std::vector<float> res_arr;

        for(int i = 0; i < length; i++){
            res_arr.emplace_back(static_cast<float>(intrm_arr[i].imag()));
        }

        return res_arr;
    }

    std::vector<float> backwardFFT(const std::vector<float> &an_seq) noexcept {

        auto length = static_cast<int>(an_seq.size());
        std::vector<std::complex<float>> intrm_arr;
        intrm_arr.resize(length);

        auto backward_FFT =  fftwf_plan_dft_1d(length, (fftwf_complex *)(intrm_arr.data()),
                                              (fftwf_complex *)(intrm_arr.data()), FFTW_BACKWARD, FFTW_ESTIMATE);

        for(int i = 0; i < length; i++){
            intrm_arr[i] = std::complex<float>(an_seq[i], 0);
        }

        fftwf_execute(backward_FFT);
        fftwf_destroy_plan(backward_FFT);
        fftwf_cleanup();

        std::vector<float> res_arr;

        for(int i = 0; i < length; i++){
            res_arr.emplace_back(static_cast<float>(intrm_arr[i].imag()));
        }

        return res_arr;
    }

/// На вход принимает массив с сигналом и сдвигает его на 90 градусов
    std::vector<float> transformHilbert(const std::vector<float>& signal) noexcept{

        int count_of_samples = static_cast<int>(signal.size());
        std::vector<std::complex<float>> buff;
        buff.resize(count_of_samples);

        auto forward_FFT =  fftwf_plan_dft_1d(count_of_samples, (fftwf_complex *)(buff.data()),
                                              (fftwf_complex *)(buff.data()), FFTW_FORWARD, FFTW_ESTIMATE);

        auto backward_FFT =  fftwf_plan_dft_1d(count_of_samples, (fftwf_complex *)(buff.data()),
                                               (fftwf_complex *)(buff.data()), FFTW_BACKWARD, FFTW_ESTIMATE);


        for(int i = 0; i < count_of_samples; i++){
            buff[i] = std::complex<float>{signal[i] , 0};
        }

        fftwf_execute(forward_FFT);

        {
            auto Im_one = std::complex<float>{0, 1};

            for (int i = -count_of_samples/2; i < count_of_samples/2; i++) {
                buff[i + count_of_samples/2] = Im_one * simpleSgn(static_cast<float>(i)) * buff[i + count_of_samples/2];
                buff[i + count_of_samples/2] /= static_cast<float>(count_of_samples);
            }
        }

        fftwf_execute(backward_FFT);

        fftwf_destroy_plan(forward_FFT);
        fftwf_destroy_plan(backward_FFT);

        fftwf_cleanup();

        poddsp::simpleSignal res_arr; res_arr.reserve(count_of_samples);
        for(auto e : buff){
            res_arr.emplace_back(e.real());
        }

        return res_arr;
    }

    std::vector<std::complex<float>> quadro_cast(const std::vector<float> & signal) noexcept{

        std::vector<float> quadro_part = transformHilbert(signal);
        std::vector<std::complex<float>> res_arr; res_arr.reserve(signal.size());

        for(int i = 0; i < signal.size(); i++){
            res_arr.emplace_back(std::complex<float> {signal[i], quadro_part[i]});
        }
        return res_arr;
    }
}