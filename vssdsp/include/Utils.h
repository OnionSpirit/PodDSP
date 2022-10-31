#pragma once

namespace vssdsp {

    typedef std::vector<std::complex<float>> c_sig_t;
    typedef std::vector<float> s_sig_t;
    typedef std::complex<float> c_smp;
    typedef float s_smp;
    const std::complex<float> J {0,1};


//    template<typename T>
//    std::vector<std::complex<float>> castToComplexSeq(const std::vector<T> &incoming_arr) {
//
//        std::vector<std::complex<float>> res_arr;
//        res_arr.reserve(incoming_arr.size());
//        if (typeid(T) == typeid(std::complex<float>) ||
//            typeid(T) == typeid(std::complex<double>) ||
//            typeid(T) == typeid(std::complex<int>) ||
//            typeid(T) == typeid(std::complex<char>)) {
//            for (auto e: incoming_arr) {
//                res_arr.emplace_back(e);
//            }
//            return res_arr;
//        }
//        if (typeid(T) == typeid(float) ||
//            typeid(T) == typeid(double) ||
//            typeid(T) == typeid(int) ||
//            typeid(T) == typeid(char)) {
//            for (auto e: incoming_arr) {
//                res_arr.emplace_back(e);
//                res_arr.back().imag(0.0f);
//            }
//        } else { throw std::invalid_argument(ERROR_CAST "Cant cast this type to complex array"); }
//        return res_arr;
//    }

    float complexVectorMagnitude(const std::complex<float> &) noexcept;

    float signalMaxValue(const std::vector<float> &) noexcept;

    float signalMinValue(const std::vector<float> &) noexcept;

    float signalMedValue(const std::vector<float> &) noexcept;

    s_sig_t signalNormalizing(const s_sig_t &) noexcept;

    c_sig_t signalNormalizing(const c_sig_t &) noexcept;

    c_sig_t forwardFFT(const c_sig_t &) noexcept;

    s_sig_t MagnitudeSpectrum(const c_sig_t &, int =0) noexcept;

    s_sig_t PhaseSpectrum(const c_sig_t &, int =0) noexcept;

    c_sig_t backwardFFT(const c_sig_t &) noexcept;

    template<typename T>
    std::complex<T> complexSgn(std::complex<T> sample){

        auto arg = complexVectorMagnitude(sample);
        if(arg == 0)
            return sample;
        return sample/fabsf(arg);
    }

    template<typename T>
    T simpleSgn(T n) {
        return (((n<0)*-1) | (n>0));
    }

    std::vector<float> transformHilbert(const std::vector<float>&) noexcept;

    std::vector<std::complex<float>> quadro_cast(const std::vector<float> &) noexcept;

    std::vector<float> AWGN_generator(size_t len) noexcept;

    s_sig_t smoother(const s_sig_t&) noexcept;

    c_sig_t cutoff(const c_sig_t&, float) noexcept;

    s_sig_t cutoff(const s_sig_t&, float) noexcept;

    s_sig_t amplifier(const s_sig_t&, float) noexcept;

    c_sig_t amplifier(const c_sig_t&, float) noexcept;

    s_sig_t moreThen(const s_sig_t&, float) noexcept;

    s_sig_t lessThen(const s_sig_t&, float) noexcept;

    float findModeWithEps(const vssdsp::s_sig_t&, float =0.01f) noexcept;

    float findModeWithEps(const vssdsp::c_sig_t&, float =0.01f) noexcept;

    float OriginalMagnitudeFind(const c_sig_t &) noexcept;

    namespace projection {

        enum type_of_projection {
            real_projection = 0,
            imaginary_projection = 1
        };

        /// Создание проекции трёхмерного графика на одну из плоскостей. Праметры:
        /// 1) Массив комплексных данных,
        /// 2) Тип проекции согласно vssdsp::projection::type_of_projection, по умолчанию real_projection.
        std::vector<float> takeProjection(const std::vector<std::complex<float>> &,
                                                 const type_of_projection & = type_of_projection::real_projection);
    }
}