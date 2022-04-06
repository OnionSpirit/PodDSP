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

    std::vector<float> FFT(const std::vector<float> & an_seq) {

        int i, j, n, m, Mmax, Istp, Nvl;
        float Tmpr, Tmpi, Wtmp, Theta;
        float Wpr, Wpi, Wr, Wi;

        std::vector<float> res_arr;
        std::vector<float> Tmvl;

        Nvl = static_cast<int>(an_seq.size());
//
        if(((Nvl << 1)>>1) != Nvl)
            Nvl--;
//
        n = Nvl * 2;
        Tmvl.reserve(n);
        res_arr.reserve(Nvl);

        for (i = 0; i < n; i+=2) {
            Tmvl.emplace_back(0.0f);
            Tmvl.emplace_back(an_seq[i / 2]);
        }

        i = 1; j = 1;
        while (i < n) {
            if (j > i) {
                std::swap(Tmvl[i],Tmvl[j]);
                std::swap(Tmvl[i+1],Tmvl[j+1]);
//                Tmpr = Tmvl[i]; Tmvl[i] = Tmvl[j]; Tmvl[j] = Tmpr;
//                Tmpr = Tmvl[i+1]; Tmvl[i+1] = Tmvl[j+1]; Tmvl[j+1] = Tmpr;
            }
            i = i + 2; m = Nvl;
            while ((m >= 2) && (j > m)) {
                j = j - m; m = m >> 1;
            }
            j = j + m;
        }

        Mmax = 2;
        while (n > Mmax) {

            Theta = static_cast<float>(-(2 * M_PI) / Mmax);
            Wpi = static_cast<float>(sin(static_cast<double>(Theta)));
            Wtmp = static_cast<float>(sin(static_cast<double>(Theta / 2.0f)));
            Wpr = Wtmp * Wtmp * 2;
            Istp = Mmax * 2;
            Wr = 1;
            Wi = 0;
            m = 1;

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


        for (i = 0; i < Nvl; i++) {
            j = i * 2;
            res_arr.emplace_back(2.0f * static_cast<float>(sqrt(pow(Tmvl[j], 2)) + pow(Tmvl[j + 1], 2)) / static_cast<float>(Nvl));
        }

        size_t tmvl_size = Tmvl.size();
        size_t res_arr_size = res_arr.size();


        return res_arr;
    }
}