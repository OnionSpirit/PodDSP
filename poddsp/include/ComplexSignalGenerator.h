#pragma once

namespace poddsp {
/// Генератор комплексного синуса. Параметры
/// 1) Количество периодов во фрейме (частота),
/// 2) Общее количество семплов,
/// 3) Начальная фаза сигнала
    std::vector<std::complex<float>> complexSin(const float &,
                                                const int &,
                                                const float & = 0);

    std::vector<float> signalShelf(const std::vector<float> &,
                                   const float &)
                                   noexcept;
/// Генератор одиночно импульса единичного уровня. Параметры
/// 1) Ширина импульса,
/// 2) Задержка,
/// 3) Величина всего колличества отсчётов.
/// Ширина и задержка импульса не должны превышать полное колличество отсчётов.
    std::vector<float> impulseGen(const int &,
                                  const int &,
                                  const int &);
/// Генератор Менандра в комплексной и обычной форме. Параметры
/// 1) Количество периодов во фрейме (частота),
/// 2) Общее количество семплов,
/// 3) Начальная фаза сигнала
    complexSignal MeanderGen(const float &,
                             const int &,
                             const float & = 0);

    simpleSignal MeanderGen(const float &,
                            const int &,
                            const float & = 0,
                            bool = false);
}
