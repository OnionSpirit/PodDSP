#pragma once

namespace poddsp {
    template<typename T>
    std::vector<std::complex<float>> castToComplexSeq(const std::vector<T> &incoming_arr) {
        if (typeid(T) == typeid(std::complex<float>) ||
            typeid(T) == typeid(std::complex<double>) ||
            typeid(T) == typeid(std::complex<int>) ||
            typeid(T) == typeid(std::complex<char>))
            return incoming_arr;
        std::vector<std::complex<float>> res_arr;
        for (auto e: incoming_arr) {
            res_arr.emplace_back(std::complex<float>{e, 0});
        }
        return res_arr;
    }

    float complexMagMeasurer(const std::complex<float> &) noexcept;

    float signalMaxValue(const std::vector<float> &) noexcept;

    float signalMinValue(const std::vector<float> &) noexcept;

    float signalMedValue(const std::vector<float> &) noexcept;
}