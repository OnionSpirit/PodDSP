#pragma once

#include <liquid/liquid.h>

namespace poddsp {
    void complexSignalResampler(const std::vector<std::complex<float>> &,
                                std::vector<std::complex<float>> &,
                                const int &)
    noexcept;
}
