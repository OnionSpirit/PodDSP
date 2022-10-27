#pragma once

namespace vssdsp {

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
/// Строит линейный график фазы
    std::vector<float> phaseDependenceLining(const std::vector<float>& phase_graph) noexcept;

/// Рассчитывает фактическую набежавшую фазу сигнала
/// вне зависимости от модуляционных фазовых скачков,
/// принимает полную линейную фазовую зависимость участка сигнала и опционально уровень скачка-ошибки
    s_sig_t phaseModulationSkipEraser(const s_sig_t&, float =0) noexcept;

/// Убирает скачки при переходе из 360 в 0, делает зависимость не ограниченной от -П до П
    std::vector<float> phaseDependenceLining(const std::vector<float>& phase_graph) noexcept;
}
