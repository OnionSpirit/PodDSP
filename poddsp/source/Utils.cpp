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

    std::vector<float> FFT(const std::vector<float> &an_seq) {

        std::vector<float> res_arr;
        int arr_length = static_cast<int>(an_seq.size());
        res_arr.reserve(arr_length);

        auto* input_arr = new double[arr_length];
        for(int i = 0; i < arr_length; i++) {
            input_arr[i] = an_seq.at(i);
        }

        auto* output_arr = new double[arr_length];

        FFTAnalysis(input_arr, output_arr, arr_length, arr_length);

        for(auto i = 0; i < arr_length; i++){
            res_arr.emplace_back(output_arr[i]);
        }


        delete [] input_arr;
        delete [] output_arr;

        res_arr.resize(arr_length/2);

        return res_arr;
    }

    void FFTAnalysis(double *AVal, double *FTvl, int Nvl, int Nft) {
        int i, j, n, m, Mmax, Istp;
        double Tmpr, Tmpi, Wtmp, Theta;
        double Wpr, Wpi, Wr, Wi;
        double *Tmvl;

        n = Nvl * 2; Tmvl = new double[n];

        for (i = 0; i < n; i+=2) {
            Tmvl[i] = 0;
            Tmvl[i+1] = AVal[i/2];
        }

        i = 1; j = 1;
        while (i < n) {
            if (j > i) {
                Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
                Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
            }
            i = i + 2; m = Nvl;
            while ((m >= 2) && (j > m)) {
                j = j - m; m = m >> 1;
            }
            j = j + m;
        }

        Mmax = 2;
        while (n > Mmax) {
            Theta = -2*M_PI / Mmax; Wpi = sin(Theta);
            Wtmp = sin(Theta / 2); Wpr = Wtmp * Wtmp * 2;
            Istp = Mmax * 2; Wr = 1; Wi = 0; m = 1;

            while (m < Mmax) {
                i = m; m = m + 2; Tmpr = Wr; Tmpi = Wi;
                Wr = Wr - Tmpr * Wpr - Tmpi * Wpi;
                Wi = Wi + Tmpr * Wpi - Tmpi * Wpr;

                while (i < n) {
                    j = i + Mmax;
                    Tmpr = Wr * Tmvl[j] - Wi * Tmvl[j-1];
                    Tmpi = Wi * Tmvl[j] + Wr * Tmvl[j-1];

                    Tmvl[j] = Tmvl[i] - Tmpr; Tmvl[j-1] = Tmvl[i-1] - Tmpi;
                    Tmvl[i] = Tmvl[i] + Tmpr; Tmvl[i-1] = Tmvl[i-1] + Tmpi;
                    i = i + Istp;
                }
            }

            Mmax = Istp;
        }

        for (i = 0; i < Nft; i++) {
            j = i * 2; FTvl[i] = 2*sqrt(pow(Tmvl[j],2) + pow(Tmvl[j+1],2))/Nvl;
        }

        delete []Tmvl;
    }

    float squareZeroPhaseSpectralFunc(float t){
        auto jr = std::complex<float>{1, 0};
        auto j = std::complex<float>{0, 1};
        return ((jr / (j * t)) * (jr - jr * exp(-j * t))).real();
    }

    float squareQuadroPhaseSpectralFunc(float t){
        auto jr = poddsp::complexSample{1, 0};
        auto j = poddsp::complexSample{0, 1};
        return ((exp(-j * t * 0.5f) / (j * t)) * (jr - jr * exp(-j * t))).real();
    }

    std::vector<float> sampleMath(int definition, float eps, float(*func)(float)){

        poddsp::simpleSignal res_arr;
        auto t = float();

        t -= eps*(static_cast<float>(definition)/2.0f);

        for(int i = 0; i < definition; i++){
            t += eps;
            res_arr.emplace_back(func(t));
        }
        return res_arr;
    }
}