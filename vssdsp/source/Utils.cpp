#include "../include/vssdsp.h"


namespace vssdsp {

    float complexVectorMagnitude(const std::complex<float> &sample) noexcept {
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

    s_sig_t MagnitudeSpectrum(const c_sig_t & sig, int move_cnt) noexcept {
        c_sig_t fft_res = forwardFFT(sig);
        auto fft_mag_ratio = signalMaxValue(projection::takeProjection(backwardFFT(fft_res))) / signalMaxValue(projection::takeProjection(sig));
        s_sig_t res; res.resize(fft_res.size() + move_cnt);
        for (int i =0; auto& e : fft_res) {
            res[i + move_cnt] = complexVectorMagnitude(e) / fft_mag_ratio;
            i++;
        }
        return res;
    }

    s_sig_t PhaseSpectrum(const c_sig_t & sig, int move_cnt) noexcept {
        c_sig_t fft_res = forwardFFT(sig);
        auto fft_mag_ratio = signalMaxValue(projection::takeProjection(backwardFFT(fft_res))) / signalMaxValue(projection::takeProjection(sig));
        s_sig_t res; res.resize(fft_res.size() + move_cnt);
        for (int i =0; auto& e : fft_res) {
            res[i + move_cnt] = complexVectorPhase(e) / fft_mag_ratio;
            i++;
        }
        return res;
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

        vssdsp::s_sig_t res_arr; res_arr.reserve(count_of_samples);
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

    c_sig_t cutoff(const c_sig_t& sig, float lvl) noexcept {
        auto res = sig;

        for (auto &e : res) {
            if(fabsf(e.real()) > lvl) {
                if(e.real() > 0) {
                    e.real(lvl);
                } else {
                    e.real(-lvl);
                }
            }
            if(fabsf(e.imag()) > lvl) {
                if(e.imag() > 0) {
                    e.imag(lvl);
                } else {
                    e.imag(-lvl);
                }
            }
        }
        return res;
    }

    s_sig_t cutoff(const s_sig_t& sig, float lvl) noexcept {
        auto res = sig;

        for (auto &e : res) {
            if(fabsf(e) > lvl) {
               if(e > 0) {
                   e = lvl;
               }
                if(e < 0) {
                    e = -lvl;
                }

            }
        }

        return res;
    }

    s_sig_t amplifier(const s_sig_t& sig, float add_lvl) noexcept {
        auto res = sig;

        for(auto& e : res) {
            e *= add_lvl;
        }
        return res;
    }

    c_sig_t amplifier(const c_sig_t& sig, float add_lvl) noexcept {
        auto res = sig;

        for(auto& e : res) {
            e.real(e.real() * add_lvl);
            e.imag(e.imag() * add_lvl);
        }
        return res;
    }

    s_sig_t moreThen(const s_sig_t& sig, float lw_lvl) noexcept {
        s_sig_t res;

        for(auto &e : sig) {
            if(e > lw_lvl) {
                res.emplace_back(e);
            }
        }
        if(res.empty()) res.emplace_back(lw_lvl);
        return res;
    }

    s_sig_t lessThen(const s_sig_t& sig, float lw_lvl) noexcept {
        s_sig_t res;

        for(auto &e : sig) {
            if(e < lw_lvl) {
                res.emplace_back(e);
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

    float findModeWithEps(const vssdsp::s_sig_t& dep, float eps) noexcept {
        using namespace vssdsp;
        std::vector<std::vector<float>> mode_st;
        auto top = signalMaxValue(dep);
        auto bot = signalMinValue(dep);
        auto lvl_range = top - bot;
        size_t iteration_count = 0;
        if ( 0.0f != eps ) iteration_count = (size_t)(lvl_range / eps);
        else return 0.0f;
        mode_st.resize(iteration_count);
        for ( int i =0; i < iteration_count and top > bot; i++ ) {
            for ( int j =0; j < dep.size(); j++ ) {
                auto el = dep[j];
                if ( (el >= (top - eps)) and (el <= (top + eps)) ) {
                    mode_st[i].emplace_back(el);
                }
            }
            top -= 2*eps;
        }
        size_t longest_id =0;
        size_t max_len =0;
        for ( int i =0; i < mode_st.size(); i++ ) {
            auto id_len = mode_st[i].size();
            if(id_len > max_len) {
                max_len = id_len;
                longest_id = i;
            }
            if ( id_len == max_len ) {
                if ( signalMedValue(mode_st[i]) > signalMedValue(mode_st[longest_id]) ) {
                    max_len = id_len;
                    longest_id = i;
                }
            }
        }
        if ( not mode_st.empty() and not mode_st[longest_id].empty()) {
            return signalMedValue(mode_st[longest_id]);
        } else {
            return signalMaxValue(dep);
        }

    }

    float findModeWithEps(const vssdsp::c_sig_t& c_dep, float eps) noexcept {
        using namespace vssdsp;
        auto dep = projection::takeProjection(c_dep);
        return findModeWithEps(dep, eps);
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