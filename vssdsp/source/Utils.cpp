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

    s_sig_t smooth(const s_sig_t& input, size_t window =-1) noexcept {
        if(window == -1) window = (size_t)((float)input.size()/20.0f);

        auto n = input.size();
        s_sig_t output; output.resize(n);
        int i, j, z, k1, k2, hw;
        double tmp;
        if (fmod(window, 2) == 0) window++;
        hw = (window - 1) / 2;
        output[0] = input[0];

        for (i = 1; i < n; i++) {
            tmp = 0;
            if (i < hw) {
                k1 = 0;
                k2 = 2 * i;
                z = k2 + 1;
            } else if ((i + hw) > (n - 1)) {
                k1 = i - n + i + 1;
                k2 = n - 1;
                z = k2 - k1 + 1;
            } else {
                k1 = i - hw;
                k2 = i + hw;
                z = window;
            }

            for (j = k1; j <= k2; j++) {
                tmp = tmp + input[j];
            }
            output[i] = tmp / z;
        }
        return output;
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

    s_sig_t genAWGN(size_t len, float mag) noexcept {

        s_sig_t res_arr; res_arr.resize(len);

        std::random_device rd{};
        std::mt19937 gen{rd()};

        std::normal_distribution<float> gss(0, mag);
        for(auto &e : res_arr) {
            e += gss(gen);
        }
        return res_arr;
    }

    float findMode(const vssdsp::s_sig_t& dep) noexcept {
        using namespace vssdsp;

        std::vector<std::vector<float>> mode_st; mode_st.resize(10);
        float Mode = 0.0f;
        auto top = signalMaxValue(dep);
        auto bot = signalMinValue(dep);
        auto lvl_range = top - bot;
        auto eps = lvl_range / 20.0f; top -= eps;
        for ( int i =0; i < 10 and top > bot; i++ ) {
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
                continue;
            }
            if ( id_len == max_len ) {
                if ( signalMedValue(mode_st[i]) > signalMedValue(mode_st[longest_id]) ) {
                    max_len = id_len;
                    longest_id = i;
                }
            }
        }
        if ( not mode_st[longest_id].empty()) {
            s_sig_t res;
            auto exp_v = signalMedValue(mode_st[longest_id]);
            auto st_d = sqrtf(dispersion(mode_st[longest_id]));
            for(auto& e : mode_st[longest_id]) {
                if((e >= exp_v - st_d) and e <= (exp_v + st_d)) {
                    res.emplace_back(e);
                }
            }
            if(res.empty()) {
                Mode = exp_v;
            } else {
                Mode = signalMedValue(res);
            }
        } else {
            Mode = signalMaxValue(dep);
        }
        return Mode;
    }

    float findMode(const vssdsp::c_sig_t& c_dep) noexcept {

        using namespace vssdsp;
        auto dep = projection::takeProjection(c_dep);
        return findMode(dep);
    }

    std::vector<size_t> topIds(const s_sig_t& dep, size_t win) noexcept {
        if(win == -1) win = (size_t)((float)dep.size()/20.0f);
        auto diff = smooth(differentiation(dep), win);
        std::vector<size_t> top_ids;
        for (int i =0; i < diff.size(); i++) {
            if (i-1 >= 0 and diff[i] <= 0 and diff[i-1] > 0) {
                top_ids.emplace_back(i);
            }
        }
        return top_ids;
    }

    std::vector<size_t> botIds(const s_sig_t& dep, size_t win) noexcept {
        if(win == -1) win = (size_t)((float)dep.size()/20.0f);
        auto diff = smooth(differentiation(dep), win);
        std::vector<size_t> bot_ids;
        for (int i =0; i < diff.size(); i++) {
            if (i-1 >= 0 and diff[i] >= 0 and diff[i-1] < 0) {
                bot_ids.emplace_back(i);
            }
        }
        return bot_ids;
    }

    /// Поиск настоящей амплитуды сигнала обрезанного границами уровня канала.
    /// Параметры: 1) Массив отсчётов комплексного сигнала.
    /// (!!!) Не может рассчитать амплитуду корректно, если обрезано более четверти амплитуды.
    /// Среднее значения допустимого SNR - 30 Db.
    float OriginalMagnitudeFind(const c_sig_t& complex_signal, size_t win) noexcept {

        s_sig_t mags;
        for(auto &e : complex_signal) {
            mags.emplace_back(complexVectorMagnitude(e));
        }
        auto top = signalMaxValue(mags), bot = signalMinValue(mags);
        auto eps = (top - bot)/10.0f;
        auto mags_interm = s_sig_t();
        if(win == -1) win = (size_t)((float)mags.size() / 20.0f);
        for (float & mag : mags) {
            if (mag <= (bot + 3.0f * eps)) {
                continue;
            }
            mags_interm.emplace_back(mag);
        } mags = mags_interm;
        mags = smooth(mags,win);

        s_sig_t tops;
        for(auto top_id : topIds(mags, win)) {
            tops.emplace_back(mags[top_id]);
        }
        return signalMedValue(tops);
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