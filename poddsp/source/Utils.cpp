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

    s_sig_t signalNormalizing(const s_sig_t & s_sig) noexcept {
        auto max = signalMaxValue(s_sig);
        s_sig_t ret;
        for(auto e : s_sig) {
            ret.emplace_back(e/max);
        }
        return ret;
    }

    c_sig_t signalNormalizing(const c_sig_t & c_sig) noexcept {
        auto i_max = signalMaxValue(projection::takeProjection(c_sig));
        auto q_max = signalMaxValue(projection::takeProjection(c_sig, projection::imaginary_projection));
        c_sig_t ret;
        for(auto e : c_sig) {
            ret.emplace_back(std::complex<float>{e.real() / i_max, e.imag() / q_max});
        }
        return ret;
    }

    c_sig_t forwardFFT(const c_sig_t &an_seq) noexcept {

        auto length = static_cast<int>(an_seq.size());
        std::vector<std::complex<float>> intrm_arr;
        intrm_arr.resize(length);

        auto forward_FFT =  fftwf_plan_dft_1d(length, (fftwf_complex *)(intrm_arr.data()),
                                              (fftwf_complex *)(intrm_arr.data()), FFTW_FORWARD, FFTW_ESTIMATE);

        intrm_arr = an_seq;
//        for(int i = 0; i < length; i++){
//            intrm_arr[i] = -std::complex<float>(an_seq[i], 0);
//        }

        fftwf_execute(forward_FFT);
        fftwf_destroy_plan(forward_FFT);
        fftwf_cleanup();

//        std::vector<float> res_arr;
//
//        for(int i = 0; i < length; i++){
//            res_arr.emplace_back(static_cast<float>(intrm_arr[i].imag()));
//        }

        return intrm_arr;
    }

    c_sig_t backwardFFT(const c_sig_t &an_seq) noexcept {

        auto length = static_cast<int>(an_seq.size());
        std::vector<std::complex<float>> intrm_arr;
        intrm_arr.resize(length);

        auto backward_FFT =  fftwf_plan_dft_1d(length, (fftwf_complex *)(intrm_arr.data()),
                                               (fftwf_complex *)(intrm_arr.data()), FFTW_BACKWARD, FFTW_ESTIMATE);
        intrm_arr = an_seq;

        fftwf_execute(backward_FFT);
        fftwf_destroy_plan(backward_FFT);
        fftwf_cleanup();

//        std::vector<float> res_arr;
//
//        for(int i = 0; i < length; i++){
//            res_arr.emplace_back(static_cast<float>(intrm_arr[i].imag()));
//        }

        return intrm_arr;
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

        poddsp::s_sig_t res_arr; res_arr.reserve(count_of_samples);
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

    s_sig_t smoother(const s_sig_t& sig) noexcept {
        auto res = sig;
        if (sig.empty()) return res;

        for (int i =0; i < sig.size(); i++) {
            if (i > 0 and i+1 < sig.size() and sig[i-1] < sig[i] and sig[i+1] < sig[i]) {
                res[i] = (sig[i-1] + sig[i+1])/2;
            }
        }

        return res;
    }


    std::vector<float> AWGN_generator(size_t len) noexcept {

        double temp1;
        double temp2;
        double result;
        int p;

        std::vector<float> res_arr;
        res_arr.reserve(len);


        srand(time(0));
        for(int i = 0; i < len; i++){
            p = 1;

            while (p > 0) {
                temp2 = (rand() / ((double) RAND_MAX));

                if (temp2 == 0) p = 1;
                else p = -1;

            }

            temp1 = cos((2.0 * (double) M_PI) * rand() / ((double) RAND_MAX));
            result = sqrt(-2.0 * log(temp2)) * temp1;

            res_arr.emplace_back(result);
            usleep(100);
        }

        {/// Нормировка
            float max_val = 0.0f;
            for (auto e: res_arr) {
                if (e > max_val) max_val = e;
            }
            for (auto &e: res_arr) {
                e /= max_val;
            }
        }

        return res_arr;
    }

    namespace projection {

        std::vector<float> takeProjection(const std::vector<std::complex<float>> &arr,
                                          const type_of_projection &type) {
            std::vector<float> res_arr;
            switch (type) {
                case 0:
                    for (auto e: arr) {
                        res_arr.emplace_back(e.real());
                    }
                    return res_arr;
                case 1:
                    for (auto e: arr) {
                        res_arr.emplace_back(e.imag());
                    }
                    return res_arr;
                default:
                    throw std::invalid_argument(ERROR_PLOT"Incorrect projection type");
            }
        }
    }
}