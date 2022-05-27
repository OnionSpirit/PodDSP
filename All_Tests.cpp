#include <gtest/gtest.h>
#include <ctime>

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

    std::cout << poddsp::complexSequenceCorrelation(original_sequence, original_sequence) << std::endl;
}

TEST(complex_functions, generating_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 10000;
    poddsp::PlotConstructor::drawPlot(poddsp::complexSin(freq, count_of_samples), "Sine " + std::to_string(freq));
}

TEST(complex_functions, resampling_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 1000;
    auto new_count_of_samples = 2500;

    std::vector<std::complex<float>> original_long_seq = poddsp::complexSin(freq, count_of_samples);
    std::vector<std::complex<float>> original_short_seq = poddsp::complexSin(freq, new_count_of_samples);
    std::vector<std::complex<float>> resampled_short_seq;

    poddsp::complexSignalResampler(original_long_seq, resampled_short_seq, new_count_of_samples);


/// Draw 3M plots

    poddsp::PlotConstructor::drawPlot(original_long_seq,
                                      (std::to_string(count_of_samples) + " samples"));
    poddsp::PlotConstructor::drawPlot(original_short_seq,
                                      (std::to_string(new_count_of_samples) + " samples"));
    poddsp::PlotConstructor::drawPlot(resampled_short_seq,
                                      (std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Real plots

    std::vector<float> original_long_seq_Re = poddsp::PlotConstructor::makeProjection(original_long_seq);
    poddsp::PlotConstructor::drawPlot(original_long_seq_Re, ("RE " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Re = poddsp::PlotConstructor::makeProjection(original_short_seq);
    poddsp::PlotConstructor::drawPlot(original_short_seq_Re,
                                      ("RE " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Re = poddsp::PlotConstructor::makeProjection(resampled_short_seq);
    poddsp::PlotConstructor::drawPlot(resampled_short_seq_Re,
                                      ("RE " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Imaginary plots

    std::vector<float> original_long_seq_Im =
            poddsp::PlotConstructor::makeProjection(original_long_seq,
                                                    poddsp::PlotConstructor::type_of_projection::imaginary_projection);
    poddsp::PlotConstructor::drawPlot(original_long_seq_Im,
                                      ("IM " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Im =
            poddsp::PlotConstructor::makeProjection(original_short_seq,
                                                    poddsp::PlotConstructor::type_of_projection::imaginary_projection);
    poddsp::PlotConstructor::drawPlot(original_short_seq_Im,
                                      ("IM " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Im =
            poddsp::PlotConstructor::makeProjection(resampled_short_seq,
                                                    poddsp::PlotConstructor::type_of_projection::imaginary_projection);
    poddsp::PlotConstructor::drawPlot(resampled_short_seq_Im,
                                      ("IM " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));

}

TEST(complex_functions, phase_calculator){

    auto freq = 5.0f;
    auto count_of_samples = 100;
    auto count_of_periods = 4;
    std::vector<std::complex<float>> sine;
    std::vector<std::complex<float>> cos;

    for (int i = 0; i < 362; i++) {
        sine = poddsp::complexSin(freq, count_of_samples);
        cos = poddsp::complexSin(freq, count_of_samples, i);
        auto phase = poddsp::complexPhaseCalculating(sine, cos);
        std::cout << phase << std::endl;
    }
}

TEST(complex_functions, phase_rotation){

    auto freq = 5.0f;
    auto count_of_samples = 1000;
    auto phase_offset_in_time = 0.0001f;
    auto phase_common_offset = 0.0f;
    std::vector<std::complex<float>> sine;
    std::vector<std::complex<float>> rotated_sine;
    std::vector<std::complex<float>> rotated_with_offset_in_time_sine;
    std::vector<std::complex<float>> disoffseted_sine;
    sine = poddsp::complexSin(freq, count_of_samples);
    rotated_sine = poddsp::complexPhaseChanger(sine, phase_common_offset);
    rotated_with_offset_in_time_sine = poddsp::complexPhaseChanger(sine, phase_common_offset, phase_offset_in_time);
    disoffseted_sine = poddsp::complexPhaseChanger(rotated_with_offset_in_time_sine, -phase_common_offset, -phase_offset_in_time);

    poddsp::PlotConstructor::drawPlot(sine, "OriginalSine");
    poddsp::PlotConstructor::drawPlot(rotated_sine, "RotatedSine");
    poddsp::PlotConstructor::drawPlot(rotated_with_offset_in_time_sine, "RotatedwithOffSine");
    poddsp::PlotConstructor::drawPlot(disoffseted_sine, "DisoffsetedSine");
}

TEST(complex_functions, complexPLL) {

    auto freq = 10.0f;
    auto count_of_samples = 2000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> sine = poddsp::complexSin(freq, count_of_samples);
//    PodDSP::PlotConstructor::drawPlot(sine, "Эталон");

    std::vector<std::complex<float>> sine_wpo = poddsp::complexPhaseChanger(sine, 0.0f, phase_attenuation_per_sample_deg);
//    PodDSP::PlotConstructor::drawPlot(sine_wpo, "Эталон с ошибкой");

    std::vector<std::complex<float>> sine_wpo_fixed = poddsp::complexPLL(sine_wpo, 1);
//    PodDSP::PlotConstructor::drawPlot(sine_wpo_fixed, "Работа ФАПЧ с ошибочным сигналом");

    std::vector<float> sine_real = poddsp::PlotConstructor::makeProjection(sine,
                                                                           poddsp::PlotConstructor::type_of_projection::real_projection);
    poddsp::PlotConstructor::drawPlot(sine_real, "Эталон (Проекция действительной части)");

    std::vector<float> sine_wpo_real = poddsp::PlotConstructor::makeProjection(sine_wpo,
                                                                               poddsp::PlotConstructor::type_of_projection::real_projection);
    poddsp::PlotConstructor::drawPlot(sine_wpo_real, "Эталон с ошибкой (Проекция действительной части)");

    std::vector<float> sine_wpo_fixed_real = poddsp::PlotConstructor::makeProjection(sine_wpo_fixed,
                                                                                     poddsp::PlotConstructor::type_of_projection::real_projection);
    poddsp::PlotConstructor::drawPlot(sine_wpo_fixed_real, "Эталон с ошибкой, исправленной ФАПЧ (Проекция действительной части)");

    std::cout << std::norm(poddsp::complexSequenceCorrelation(sine, sine_wpo_fixed));
}

TEST(complex_functions, complexPLL_with_modulated_signal){
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = 0.001f;

    std::vector<std::complex<float>> complex_carrier = poddsp::complexSin(freq, count_of_samples, 0);
    std::vector<float> complex_carrier_real = poddsp::PlotConstructor::makeProjection(complex_carrier,
                                                                                      poddsp::PlotConstructor::type_of_projection::real_projection);
    std::vector<float> mag_modulation = poddsp::PlotConstructor::makeProjection(
            poddsp::complexSin(info_freq, info_count_of_samples, -90));
    poddsp::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = poddsp::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_real = poddsp::PlotConstructor::makeProjection(modulated_carrier, poddsp::PlotConstructor::type_of_projection::imaginary_projection);
    std::vector<float> modulated_carrier_imag = poddsp::PlotConstructor::makeProjection(modulated_carrier);

    poddsp::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Проекция действительной части)");
    poddsp::PlotConstructor::drawPlot(modulated_carrier, "Амплитудно модулированный сигнал");
    poddsp::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
    poddsp::PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    std::vector<std::complex<float>> modulated_carrier_wpo = poddsp::complexPhaseChanger(modulated_carrier, 0.0f,
                                                                                         phase_attenuation_per_sample_deg);
    std::vector<float> modulated_carrier_wpo_real = poddsp::PlotConstructor::makeProjection(modulated_carrier_wpo);
    poddsp::PlotConstructor::drawPlot(modulated_carrier_wpo_real, "Амплитудно модулированный сигнал с ошибкой (Проекция действительной части)");
//
    std::vector<std::complex<float>> modulated_carrier_wpo_fixed = poddsp::complexPLL(modulated_carrier_wpo, 1);
    std::vector<float> modulated_carrier_wpo_fixed_real = poddsp::PlotConstructor::makeProjection(modulated_carrier_wpo_fixed);
    poddsp::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_real, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция действительной части)");
//    std::cout << std::norm(poddsp::complexSequenceCorrelation(modulated_carrier, modulated_carrier_wpo_fixed));

    std::vector<float> mag_demodulation_result = poddsp::complexMagDemodulator(modulated_carrier_wpo_fixed);
    poddsp::PlotConstructor::drawPlot(mag_demodulation_result, "Информационный сигнал из демодулятора");
}

TEST(complex_functions, BPSK_Modulation){
    auto freq = 10.0f;

    std::vector<bool> info_signal;
    for(int i = 0; i < static_cast<int>(freq); i++){
        if(RANDOM_NUMBER >= 5){
            info_signal.emplace_back(true);
        }else{
            info_signal.emplace_back(false);
        }
    }

    std::vector<std::complex<float>> bpsk_modulated = poddsp::complexBPSKModulator(info_signal, 100);
    std::vector<float> bpsk_modulated_real = poddsp::PlotConstructor::makeProjection(bpsk_modulated);
    std::vector<float> bpsk_modulated_imag = poddsp::PlotConstructor::makeProjection(bpsk_modulated,
                                                                                     poddsp::PlotConstructor::type_of_projection::imaginary_projection);

    poddsp::PlotConstructor::drawPlot(bpsk_modulated, "BPSK модулированный сигнал");
//    poddsp::PlotConstructor::drawPlot(poddsp::complexSignalPhaseDependence(bpsk_modulated), "Зависимость фазы");
    poddsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
    poddsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");

    for(auto e : info_signal){
        std::cout << e << "\t";
    }
}

TEST(complex_functions, BPSK_with_PLL){

/// ToDo Fix BPSK logic fully

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
                                  false};
//    auto count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = 0.003f;

//    std::vector<std::complex<float>> complex_carrier = poddsp::complexSin(freq, count_of_samples, -90);
//    std::vector<float> complex_carrier_real = poddsp::PlotConstructor::makeProjection(complex_carrier,
//                                                                                      poddsp::PlotConstructor::type_of_projection::real_projection);
//
//    std::vector<float> complex_carrier_imag = poddsp::PlotConstructor::makeProjection(complex_carrier,
//                                                                                      poddsp::PlotConstructor::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated = poddsp::complexBPSKModulator(info_signal);
    std::vector<float> bpsk_modulated_real = poddsp::PlotConstructor::makeProjection(bpsk_modulated,
                                                                                     poddsp::PlotConstructor::type_of_projection::real_projection);

    std::vector<float> bpsk_modulated_imag = poddsp::PlotConstructor::makeProjection(bpsk_modulated,
                                                                                     poddsp::PlotConstructor::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo = poddsp::complexPhaseChanger(bpsk_modulated, 0.0f,
                                                                                      phase_attenuation_per_sample_deg);

    std::vector<float> bpsk_modulated_wpo_real = poddsp::PlotConstructor::makeProjection(bpsk_modulated_wpo,
                                                                                         poddsp::PlotConstructor::type_of_projection::real_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo_fixed = poddsp::complexPLL(bpsk_modulated_wpo);

    std::vector<float> bpsk_modulated_wpo_fixed_real = poddsp::PlotConstructor::makeProjection(bpsk_modulated_wpo_fixed,
                                                                                               poddsp::PlotConstructor::type_of_projection::real_projection);

//    std::vector<float> bpsk_phase_function = poddsp::complexSignalPhaseDependence(bpsk_modulated);


//    PodDSP::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Действительная часть)");
//    PodDSP::PlotConstructor::drawPlot(complex_carrier_imag, "Несущий сигнал (Мнимая часть)");
//    PodDSP::PlotConstructor::drawPlot(PodDSP::complexSignalPhaseDependence(complex_carrier), "Изменение фазы в несущем сигнале");

    poddsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
    poddsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");
//    poddsp::PlotConstructor::drawPlot(bpsk_phase_function, "Изменение фазы в BPSK сигнале");
//
//    poddsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_real, "BPSK модулированный сигнал с ошибкой");
//    poddsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_fixed_real, "BPSK модулированный сигнал с исправленной ошибкой");

    std::cout << std::norm(poddsp::complexSequenceCorrelation(bpsk_modulated, bpsk_modulated_wpo_fixed));
}

TEST(transform, FFT){

    auto count_of_samples = 2000;
    auto freq = 16.0f;

    poddsp::simpleSignal impulse = poddsp::MeanderGen(freq, count_of_samples, 0, true);


    std::vector<float> FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);
    FFT_analysis_result = poddsp::forwardFFT(impulse);


    poddsp::PlotConstructor::drawPlot(impulse, "Меандр");
    poddsp::PlotConstructor::drawPlot(FFT_analysis_result, "Спектр сигнала");

    FFT_analysis_result = poddsp::backwardFFT(FFT_analysis_result);

    poddsp::PlotConstructor::drawPlot(FFT_analysis_result, "Опять сигнал");


}

TEST(transform, specturm_plots){

    auto freq = 4.0f;
    auto count_of_samples = 2000;

    poddsp::simpleSignal Sine = poddsp::PlotConstructor::makeProjection(poddsp::complexSin(freq, count_of_samples));
    poddsp::PlotConstructor::drawPlot(Sine, "Исходный сигнал");

    auto Rotated_sine = poddsp::forwardFFT(Sine);
    poddsp::PlotConstructor::drawPlot(Rotated_sine, "Повёрнутый сигнал");
}


TEST(transform, hilbert_transform){
    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    poddsp::simpleSignal Sine = poddsp::PlotConstructor::makeProjection(poddsp::complexSin(freq, count_of_samples, 0));
    poddsp::simpleSignal Hilberted = poddsp::transformHilbert(Sine);

    poddsp::PlotConstructor::drawPlot(Sine, "Source sine");
    poddsp::PlotConstructor::drawPlot(Hilberted, "Hilbert transformed sine");

}

TEST(transform, quadro_cast){

    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    poddsp::complexSignal Sine = poddsp::complexSin(freq, count_of_samples, 0);
    poddsp::simpleSignal pSine = poddsp::PlotConstructor::makeProjection(Sine);

    poddsp::complexSignal CSine = poddsp::quadro_cast(pSine);
    poddsp::PlotConstructor::drawPlot(pSine, "Source sine complex");
    poddsp::PlotConstructor::drawPlot(Sine, "Source sine projection");
    poddsp::PlotConstructor::drawPlot(CSine, "Complex source sine");
}

TEST(transform, fftw_speed_test){
    constexpr auto count_of_samples = 8192;
    constexpr auto freq = 16.0f;
    constexpr auto count_of_processing = 1000000;

    poddsp::complexSignal FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);

    auto forward_FFT =  fftwf_plan_dft_1d(count_of_samples, (fftwf_complex *)(FFT_analysis_result.data()),
                                          (fftwf_complex *)(FFT_analysis_result.data()), FFTW_FORWARD, FFTW_MEASURE);

    for(auto e : poddsp::MeanderGen(freq, count_of_samples, 0, true)){
        FFT_analysis_result.emplace_back(std::complex<float>{e, 0});
    }

//    auto gen_error = 0.0f;
//    auto interm_clock = 0.0f;
    auto start_time = clock();
    for(int i = count_of_processing; i != 0; i--){
        fftwf_execute(forward_FFT);
    }

    auto total_time = clock() - start_time;

    fftwf_destroy_plan(forward_FFT);
    fftwf_cleanup();

    std::cout << "Total time in seconds: " << static_cast<float>(total_time) / 1000000.0f << std::endl;
}

TEST(generators, AWGN_generator){

    poddsp::simpleSignal noise = poddsp::AWGN_generator(1000);

    poddsp::complexSignal  complex_noise = poddsp::quadro_cast(poddsp::AWGN_generator(1000));

    poddsp::PlotConstructor::drawPlot(complex_noise, "АГБШ");
}