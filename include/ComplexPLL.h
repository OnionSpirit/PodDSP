#pragma once

namespace PodDSP {

    std::vector<std::complex<float>> complexCFOCompensator(const std::vector<std::complex<float>> &,
                                                           const float &)
                                                           noexcept;

    std::vector<float> differentiation(const std::vector<float> &);

/// ФАПЧ для комплексного сигнала.
/// Параметры:
/// 1) Массив сигнала для обработки,
/// 2) Ширина участка сигнала для взятия производной изменения аргумента семплов, по умолчанию вся длина сигнала,
/// 3) Количество проходов ФАПЧ по всей длине сигнала, по умолчанию 1.
/// Длину участка дифференцирования, в количестве семплов,
/// рекомендуется указать максимально близкой к предполагаемому
/// количеству семплов приходящихся
/// на период самой низкочастотной составляющей исходного сигнала,
/// чем точнее к действительности тем точнее работа ФАПЧ.
/// Ширина окна дифференцирования так же должна быть кратна количеству семплов всего обрабатываемого сигнала.
/// Количество проходов ФАПЧ не следует выбирать слишком большим при использовании угловых модуляций.
/// При использовании на АМ сигнале глубина его модуляции не должна превышать 90%.
    std::vector<std::complex<float>> complexPLL(const std::vector<std::complex<float>> &,
                                                int = 1,
                                                int = 1);
}
