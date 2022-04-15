#pragma once

namespace poddsp {

    using complexSignal = std::vector<std::complex<float>>;
    using simpleSignal = std::vector<float>;
    using complexSample = std::complex<float>;


    template<typename T>
    std::vector<std::complex<float>> castToComplexSeq(const std::vector<T> &incoming_arr) {

        std::vector<std::complex<float>> res_arr;
        res_arr.reserve(incoming_arr.size());
        if (typeid(T) == typeid(std::complex<float>) ||
            typeid(T) == typeid(std::complex<double>) ||
            typeid(T) == typeid(std::complex<int>) ||
            typeid(T) == typeid(std::complex<char>)) {
            for (auto e: incoming_arr) {
                res_arr.emplace_back(e);
            }
            return res_arr;
        }
        if (typeid(T) == typeid(float) ||
            typeid(T) == typeid(double) ||
            typeid(T) == typeid(int) ||
            typeid(T) == typeid(char)) {
            for (auto e: incoming_arr) {
                res_arr.emplace_back(e);
                res_arr.back().imag(0.0f);
            }
        } else { throw std::invalid_argument(ERROR_CAST "Cant cast this type to complex array"); }
        return res_arr;
    }

    float complexMagMeasurer(const std::complex<float> &) noexcept;

    float signalMaxValue(const std::vector<float> &) noexcept;

    float signalMinValue(const std::vector<float> &) noexcept;

    float signalMedValue(const std::vector<float> &) noexcept;

    std::vector<float> forwardFFT(const std::vector<float> &);

    std::vector<float> backwardFFT(const std::vector<float> &);

    template<typename T>
    std::complex<T> complexSgn(std::complex<T>);

    template<typename T>
    T simpleSgn(T);

    std::vector<float> transformHilbert(const std::vector<float>&);
}