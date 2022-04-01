#include "../include/poddsp.h"


namespace poddsp{
    
    float complexMagMeasurer(const std::complex<float> & sample) noexcept{
        float res;
        res = sqrt(pow(sample.real(),2) + pow(sample.imag(),2));
        return res;
    }

    float signalMaxValue(const std::vector<float> & arr) noexcept {
        float res = arr[0];
        for(auto e : arr){
            if(res <= e){
                res = e;
            }
        }
        return res;
    }
    float signalMinValue(const std::vector<float> & arr) noexcept {
        float res = arr[0];
        for(auto e : arr){
            if(res >= e){
                res = e;
            }
        }
        return res;
    }
    float signalMedValue(const std::vector<float> & arr) noexcept {
        float res = 0.0f;
        for(auto e : arr){
            res += e;
        }
        res /= static_cast<float>(arr.size());
        return res;
    }

    std::vector<float> complexSignalPhaseDependence(const std::vector<std::complex<float>> &signal){

        std::vector<float> res_arr;
        for (auto e: signal) {
            res_arr.emplace_back(complexVectorPhase(e));
        }
        return res_arr;
    }
}