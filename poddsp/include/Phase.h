#pragma once

namespace poddsp {

/// Вращение фазы сигнала против часовой стрелки + CFO gen
    std::vector<std::complex<float>> complexPhaseChanger(const std::vector<std::complex<float>> &incoming_arr,
                                                         const float &additional_phase_deg,
                                                         const float &phase_attenuation_per_sample_deg = 0)
    noexcept;

    std::complex<float> complexIntermediatePhaseCalculating(const std::vector<std::complex<float>> &);
/// Вычисление разности фаз между первым и вторым комплексным сигналом (в градусах)
    float complexPhaseCalculating(const std::vector<std::complex<float>> &,
                                  const std::vector<std::complex<float>> &)
    noexcept;

    std::vector<float> complexSignalPhaseDependence(const std::vector<std::complex<float>> &);

/// Вычисление фазы комплексного вектора
    float complexVectorPhase(const std::complex<float> &) noexcept;

}
