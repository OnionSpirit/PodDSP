#include "../include/vssdsp.h"

namespace vssdsp {

    inline c_sig_t Heterodyne(float freq, const c_sig_t & converting_signal, bool rise) noexcept {

        auto len = converting_signal.size();
        auto ref_freq = complexSin(freq, len);
        c_sig_t res_arr; res_arr.resize(len);

        if (rise) {
            for(int i = 0; i < len; i++) {

                res_arr[i] = converting_signal[i] * ref_freq[i];
            }
        } else {
            for(int i = 0; i < len; i++) {

                res_arr[i] = converting_signal[i] / ref_freq[i];
            }
        }

        return res_arr;
    }
}