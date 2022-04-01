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
}
