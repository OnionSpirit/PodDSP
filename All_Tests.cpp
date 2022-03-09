#include <gtest/gtest.h>
#include <poddsp.h>

#define RANDOM_NUMBER 1 + rand()%10


TEST(complex_functions, complex_correlation_calculation){

    int sequence_length = 10;
    srand(time(nullptr));

    std::vector<std::complex<float>> original_sequence;
    std::vector<std::complex<float>> incoming_sequence;

    for(int i = 0; i < sequence_length; i++){
        original_sequence.emplace_back(RANDOM_NUMBER); original_sequence.back().imag(RANDOM_NUMBER);
        incoming_sequence.emplace_back(RANDOM_NUMBER); incoming_sequence.back().imag(RANDOM_NUMBER);
    }

    std::cout << PodDSP::complexSequenceCorrelation(original_sequence, original_sequence) << std::endl;
}

TEST(complex_functions, generating_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 10000;
    PodDSP::PlotConstructor::drawPlot(PodDSP::complexSin(freq, count_of_samples), "Sine " + std::to_string(freq));
}

TEST(complex_functions, resampling_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 1000;
    auto new_count_of_samples = 2500;

    std::vector<std::complex<float>> original_long_seq = PodDSP::complexSin(freq, count_of_samples);
    std::vector<std::complex<float>> original_short_seq = PodDSP::complexSin(freq, new_count_of_samples);
    std::vector<std::complex<float>> resampled_short_seq;

    PodDSP::complexSignalResampler(original_long_seq, resampled_short_seq, new_count_of_samples);


/// Draw 3M plots

    PodDSP::PlotConstructor::drawPlot(original_long_seq,
                    (std::to_string(count_of_samples) + " samples"));
    PodDSP::PlotConstructor::drawPlot(original_short_seq,
                    (std::to_string(new_count_of_samples) + " samples"));
    PodDSP::PlotConstructor::drawPlot(resampled_short_seq,
                    (std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Real plots

    std::vector<float> original_long_seq_Re = PodDSP::PlotConstructor::makeProjection(original_long_seq);
    PodDSP::PlotConstructor::drawPlot(original_long_seq_Re, ("RE " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Re = PodDSP::PlotConstructor::makeProjection(original_short_seq);
    PodDSP::PlotConstructor::drawPlot(original_short_seq_Re,
             ("RE " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Re = PodDSP::PlotConstructor::makeProjection(resampled_short_seq);
    PodDSP::PlotConstructor::drawPlot(resampled_short_seq_Re,
             ("RE " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Imaginary plots

    std::vector<float> original_long_seq_Im =
            PodDSP::PlotConstructor::makeProjection(original_long_seq,
                                                    PodDSP::PlotConstructor::type_of_projection::imaginary_projection);
    PodDSP::PlotConstructor::drawPlot(original_long_seq_Im,
             ("IM " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Im =
            PodDSP::PlotConstructor::makeProjection(original_short_seq,
                                                    PodDSP::PlotConstructor::type_of_projection::imaginary_projection);
    PodDSP::PlotConstructor::drawPlot(original_short_seq_Im,
             ("IM " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Im =
            PodDSP::PlotConstructor::makeProjection(resampled_short_seq,
                                                    PodDSP::PlotConstructor::type_of_projection::imaginary_projection);
    PodDSP::PlotConstructor::drawPlot(resampled_short_seq_Im,
             ("IM " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));

}

TEST(complex_functions, phase_calculator){

    auto freq = 5.0f;
    auto count_of_samples = 100;
    auto count_of_periods = 4;
    std::vector<std::complex<float>> sine;
    std::vector<std::complex<float>> cos;

    for (int i = 0; i < 362; i++) {
        sine = PodDSP::complexSin(freq, count_of_samples);
        cos = PodDSP::complexSin(freq, count_of_samples, i);
        auto phase = PodDSP::complexPhaseCalculating(sine, cos);
        std::cout << phase << std::endl;
    }
}

TEST(complex_functions, phase_rotation){

    auto freq = 5.0f;
    auto count_of_samples = 10000;
    auto phase_offset_in_time = 0.0001f;
    auto phase_common_offset = 0.0f;
    std::vector<std::complex<float>> sine;
    std::vector<std::complex<float>> rotated_sine;
    std::vector<std::complex<float>> rotated_with_offset_in_time_sine;
    std::vector<std::complex<float>> disoffseted_sine;
    sine = PodDSP::complexSin(freq, count_of_samples);
    rotated_sine = PodDSP::complexPhaseChanger(sine, phase_common_offset);
    rotated_with_offset_in_time_sine = PodDSP::complexPhaseChanger(sine, phase_common_offset, phase_offset_in_time);
    disoffseted_sine = PodDSP::complexPhaseChanger(rotated_with_offset_in_time_sine, -phase_common_offset, -phase_offset_in_time);

    PodDSP::PlotConstructor::drawPlot(sine, "OriginalSine");
    PodDSP::PlotConstructor::drawPlot(rotated_sine, "RotatedSine");
    PodDSP::PlotConstructor::drawPlot(rotated_with_offset_in_time_sine, "RotatedwithOffSine");
    PodDSP::PlotConstructor::drawPlot(disoffseted_sine, "DisoffsetedSine");
}

TEST(complex_functions, complexPLL) {

    auto freq = 10.0f;
    auto count_of_samples = 2000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> sine = PodDSP::complexSin(freq, count_of_samples,2);
//    PodDSP::PlotConstructor::drawPlot(sine, "Эталон");

    std::vector<std::complex<float>> sine_wpo = PodDSP::complexPhaseChanger(sine, 0.0f, phase_attenuation_per_sample_deg);
//    PodDSP::PlotConstructor::drawPlot(sine_wpo, "Эталон с ошибкой");

    std::vector<std::complex<float>> sine_wpo_fixed = PodDSP::complexPLL(sine_wpo, count_of_samples / (int)freq, 3);
//    PodDSP::PlotConstructor::drawPlot(sine_wpo_fixed, "Работа ФАПЧ с ошибочным сигналом");

    std::vector<float> sine_real = PodDSP::PlotConstructor::makeProjection(sine,
                                                                           PodDSP::PlotConstructor::type_of_projection::real_projection);
    PodDSP::PlotConstructor::drawPlot(sine_real, "Эталон (Проекция действительной части)");

    std::vector<float> sine_wpo_real = PodDSP::PlotConstructor::makeProjection(sine_wpo,
                                                                               PodDSP::PlotConstructor::type_of_projection::real_projection);
    PodDSP::PlotConstructor::drawPlot(sine_wpo_real, "Эталон с ошибкой (Проекция действительной части)");

    std::vector<float> sine_wpo_fixed_real = PodDSP::PlotConstructor::makeProjection(sine_wpo_fixed,
                                                                                     PodDSP::PlotConstructor::type_of_projection::real_projection);
    PodDSP::PlotConstructor::drawPlot(sine_wpo_fixed_real, "Эталон с ошибкой, исправленной ФАПЧ (Проекция действительной части)");

    std::cout << std::norm(PodDSP::complexSequenceCorrelation(sine, sine_wpo_fixed));
}

TEST(complex_functions, complexPLL_with_modulated_signal){
    auto freq = 40.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 4000;
    auto info_count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> complex_carrier = PodDSP::complexSin(freq, count_of_samples, -90);
    std::vector<float> complex_carrier_real = PodDSP::PlotConstructor::makeProjection(complex_carrier,
                                                                                      PodDSP::PlotConstructor::type_of_projection::real_projection);
    std::vector<float> mag_modulation = PodDSP::PlotConstructor::makeProjection(
                    PodDSP::complexSin(info_freq, info_count_of_samples, -90));
    PodDSP::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = PodDSP::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_real = PodDSP::PlotConstructor::makeProjection(modulated_carrier);

    PodDSP::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Проекция действительной части)");
    PodDSP::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");

    std::vector<std::complex<float>> modulated_carrier_wpo = PodDSP::complexPhaseChanger(modulated_carrier, 0.0f,
                                                                                         phase_attenuation_per_sample_deg);
    std::vector<float> modulated_carrier_wpo_real = PodDSP::PlotConstructor::makeProjection(modulated_carrier_wpo);
    PodDSP::PlotConstructor::drawPlot(modulated_carrier_wpo_real, "Амплитудно модулированный сигнал с ошибкой (Проекция действительной части)");

    std::vector<std::complex<float>> modulated_carrier_wpo_fixed = PodDSP::complexPLL(modulated_carrier_wpo, count_of_samples / (int)info_freq, 3);
    std::vector<float> modulated_carrier_wpo_fixed_real = PodDSP::PlotConstructor::makeProjection(modulated_carrier_wpo_fixed);
    PodDSP::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_real, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция действительной части)");
    std::cout << std::norm(PodDSP::complexSequenceCorrelation(modulated_carrier, modulated_carrier_wpo_fixed));

    std::vector<float> mag_demodulation_result = PodDSP::complexMagDemodulator(modulated_carrier_wpo_fixed);
    PodDSP::PlotConstructor::drawPlot(mag_demodulation_result, "Информационный сигнал из демодулятора");
}

TEST(complex_functions, BPSK){
    auto freq = 10.0f;
    std::vector<bool> info_signal{false,
                                  true,
                                  false,
                                  false,
                                  false,
                                  false,
                                  true,
                                  true,
                                  true,
                                  false,
                                  true};
    auto count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> complex_carrier = PodDSP::complexSin(freq, count_of_samples, -90);
    std::vector<float> complex_carrier_real = PodDSP::PlotConstructor::makeProjection(complex_carrier,
                                                                                      PodDSP::PlotConstructor::type_of_projection::real_projection);

    std::vector<float> complex_carrier_imag = PodDSP::PlotConstructor::makeProjection(complex_carrier,
                                                                                      PodDSP::PlotConstructor::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated = PodDSP::complexBPSKModulator(complex_carrier, info_signal);
    std::vector<float> bpsk_modulated_real = PodDSP::PlotConstructor::makeProjection(bpsk_modulated,
                                                                                      PodDSP::PlotConstructor::type_of_projection::real_projection);

    std::vector<float> bpsk_modulated_imag = PodDSP::PlotConstructor::makeProjection(bpsk_modulated,
                                                                                     PodDSP::PlotConstructor::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo = PodDSP::complexPhaseChanger(bpsk_modulated, 0.0f,
                                                                                      phase_attenuation_per_sample_deg);

    std::vector<float> bpsk_modulated_wpo_real = PodDSP::PlotConstructor::makeProjection(bpsk_modulated_wpo,
                                                                                     PodDSP::PlotConstructor::type_of_projection::real_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo_fixed = PodDSP::complexPLL(bpsk_modulated_wpo, count_of_samples/(int)freq);

    std::vector<float> bpsk_modulated_wpo_fixed_real = PodDSP::PlotConstructor::makeProjection(bpsk_modulated_wpo_fixed,
                                                                                         PodDSP::PlotConstructor::type_of_projection::real_projection);

    std::vector<float> bpsk_phase_function = PodDSP::complexSignalPhaseDependence(bpsk_modulated);


//    PodDSP::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Действительная часть)");
//    PodDSP::PlotConstructor::drawPlot(complex_carrier_imag, "Несущий сигнал (Мнимая часть)");
//    PodDSP::PlotConstructor::drawPlot(PodDSP::complexSignalPhaseDependence(complex_carrier), "Изменение фазы в несущем сигнале");

    PodDSP::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
    PodDSP::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");
    PodDSP::PlotConstructor::drawPlot(bpsk_phase_function, "Изменение фазы в BPSK сигнале");
//
//    PodDSP::PlotConstructor::drawPlot(bpsk_modulated_wpo_real, "BPSK модулированный сигнал с ошибкой");
//    PodDSP::PlotConstructor::drawPlot(bpsk_modulated_wpo_fixed_real, "BPSK модулированный сигнал с исправленной ошибкой");

    std::cout << std::norm(PodDSP::complexSequenceCorrelation(bpsk_modulated, bpsk_modulated_wpo_fixed));
}
