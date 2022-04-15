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

    template<typename T>
    std::complex<T> complexSgn(std::complex<T> sample){

        auto arg = complexMagMeasurer(sample);
        if(arg == 0)
            return sample;
        return sample/fabsf(arg);
    }

    template<typename T>
    T simpleSgn(T n) {
        return (((n<0)*-1) | (n>0));
    }

    std::vector<float> forwardFFT(const std::vector<float> &an_seq) {

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

    std::vector<float> backwardFFT(const std::vector<float> &an_seq) {

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
    std::vector<float> transformHilbert(const std::vector<float>& signal){

        int count_of_samples = static_cast<int>(signal.size());
        std::vector<std::complex<float>> intrm_arr; intrm_arr.reserve(count_of_samples);
        std::vector<float> spectrum = forwardFFT(signal);

        {
            auto Im_one = std::complex<float>{0, 1};

            for (int i = -count_of_samples/2; i < count_of_samples/2; i++) {
                intrm_arr.emplace_back(-Im_one * simpleSgn(static_cast<float>(i))* spectrum[i + count_of_samples/2]);
                intrm_arr.back() /= static_cast<float>(count_of_samples);
            }
        }
        simpleSignal res_arr; res_arr.reserve(count_of_samples);
/// нет записи
        for(auto e : intrm_arr){
            res_arr.emplace_back(e.real());
        }

        return res_arr;
    }
}