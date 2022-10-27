#include <gtest/gtest.h>
#include <ctime>

#include <vssdsp.h>

#define RANDOM_NUMBER 1 + rand()%10
#define DO_PLOTS false


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
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(poddsp::complexSin(freq, count_of_samples), "Sine " + std::to_string(freq));
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

//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_long_seq,
//                                      (std::to_string(count_of_samples) + " samples"));
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_short_seq,
//                                      (std::to_string(new_count_of_samples) + " samples"));
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(resampled_short_seq,
//                                      (std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Real plots

    std::vector<float> original_long_seq_Re = poddsp::projection::takeProjection(original_long_seq);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_long_seq_Re, ("RE " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Re = poddsp::projection::takeProjection(original_short_seq);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_short_seq_Re,
                                      ("RE " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Re = poddsp::projection::takeProjection(resampled_short_seq);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(resampled_short_seq_Re,
                                      ("RE " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Imaginary plots

//    std::vector<float> original_long_seq_Im =
//            poddsp::projection::takeProjection(original_long_seq,
//                                                    poddsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_long_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " samples"));
//
//    std::vector<float> original_short_seq_Im =
//            poddsp::projection::takeProjection(original_short_seq,
//                                                    poddsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(original_short_seq_Im,
//                                      ("IM " + std::to_string(new_count_of_samples) + " samples"));
//
//    std::vector<float> resampled_short_seq_Im =
//            poddsp::projection::takeProjection(resampled_short_seq,
//                                                    poddsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(resampled_short_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));

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

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine, "OriginalSine");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(rotated_sine, "RotatedSine");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(rotated_with_offset_in_time_sine, "RotatedwithOffSine");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(disoffseted_sine, "DisoffsetedSine");
}

TEST(complex_functions, complexPLL) {

    auto freq = 10.0f;
    auto count_of_samples = 2000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> sine = poddsp::complexSin(freq, count_of_samples);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine, "Эталон");

    std::vector<std::complex<float>> sine_wpo = poddsp::complexPhaseChanger(sine, 0.0f, phase_attenuation_per_sample_deg);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine_wpo, "Эталон с ошибкой");

    std::vector<std::complex<float>> sine_wpo_fixed = poddsp::complexPLL(sine_wpo, freq);
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine_wpo_fixed, "Работа ФАПЧ с ошибочным сигналом");

    std::vector<float> sine_real = poddsp::projection::takeProjection(sine,
                                                                           poddsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine_real, "Эталон (Проекция действительной части)");

    std::vector<float> sine_wpo_real = poddsp::projection::takeProjection(sine_wpo,
                                                                               poddsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine_wpo_real, "Эталон с ошибкой (Проекция действительной части)");

    std::vector<float> sine_wpo_fixed_real = poddsp::projection::takeProjection(sine_wpo_fixed,
                                                                                     poddsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(sine_wpo_fixed_real, "Эталон с ошибкой, исправленной ФАПЧ (Проекция действительной части)");

    std::cout << std::norm(poddsp::complexSequenceCorrelation(sine, sine_wpo_fixed));
}

TEST(complex_functions, complexPLL_with_modulated_signal){
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = 0.001f;

    std::vector<std::complex<float>> complex_carrier = poddsp::complexSin(freq, count_of_samples, 0);
    std::vector<float> complex_carrier_real = poddsp::projection::takeProjection(complex_carrier,
                                                                                      poddsp::projection::type_of_projection::real_projection);
    std::vector<float> mag_modulation = poddsp::projection::takeProjection(
            poddsp::complexSin(info_freq, info_count_of_samples));
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = poddsp::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_imag = poddsp::projection::takeProjection(modulated_carrier, poddsp::projection::type_of_projection::imaginary_projection);
    std::vector<float> modulated_carrier_real = poddsp::projection::takeProjection(modulated_carrier);

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    std::vector<std::complex<float>> modulated_carrier_wpo = poddsp::complexPhaseChanger(modulated_carrier, 0.0f,
                                                                                         phase_attenuation_per_sample_deg);
    std::vector<float> modulated_carrier_wpo_real = poddsp::projection::takeProjection(modulated_carrier_wpo);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_wpo_real, "Амплитудно модулированный сигнал с ошибкой (Проекция действительной части)");
//
    std::vector<std::complex<float>> modulated_carrier_wpo_fixed = poddsp::complexPLL(modulated_carrier_wpo, freq);
    std::vector<float> modulated_carrier_wpo_fixed_real = poddsp::projection::takeProjection(modulated_carrier_wpo_fixed);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_real, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция действительной части)");
    std::vector<float> modulated_carrier_wpo_fixed_imag = poddsp::projection::takeProjection(modulated_carrier_wpo_fixed, poddsp::projection::imaginary_projection);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_imag, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция мнимой части)");


    std::vector<float> mag_demodulation_result = poddsp::complexMagDemodulator(modulated_carrier_wpo_fixed);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(mag_demodulation_result, "Информационный сигнал из демодулятора");

    std::cout << "Correlation " << std::norm(poddsp::complexSequenceCorrelation(modulated_carrier_wpo_fixed, modulated_carrier)) << std::endl;
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
    std::vector<float> bpsk_modulated_real = poddsp::projection::takeProjection(bpsk_modulated);
    std::vector<float> bpsk_modulated_imag = poddsp::projection::takeProjection(bpsk_modulated,
                                                                                     poddsp::projection::type_of_projection::imaginary_projection);

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated, "BPSK модулированный сигнал");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(poddsp::complexSignalPhaseDependence(bpsk_modulated), "Зависимость фазы");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");

    for(auto e : info_signal){
        std::cout << e << "\t";
    }
}

TEST(complex_functions, BPSK_with_PLL){

/// ToDo Fix BPSK logic fully

//    auto freq = 11.0f;
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
                                  false,
                                  true,
                                  false,
                                  false,
                                  false,
                                  false,
                                  true,
                                  true,
                                  true,};
    auto freq = info_signal.size();
//    auto count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = -0.003f;

//    std::vector<std::complex<float>> complex_carrier = poddsp::complexSin(freq, count_of_samples, -90);
//    std::vector<float> complex_carrier_real = poddsp::projection::takeProjection(complex_carrier,
//                                                                                      poddsp::projection::type_of_projection::real_projection);
//
//    std::vector<float> complex_carrier_imag = poddsp::projection::takeProjection(complex_carrier,
//                                                                                      poddsp::projection::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated = poddsp::complexBPSKModulator(info_signal, 200);
    std::vector<float> bpsk_modulated_real = poddsp::projection::takeProjection(bpsk_modulated,
                                                                                poddsp::projection::type_of_projection::real_projection);

    std::vector<float> bpsk_modulated_imag = poddsp::projection::takeProjection(bpsk_modulated,
                                                                                poddsp::projection::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo = poddsp::complexPhaseChanger(bpsk_modulated, 0.0f,
                                                                                      phase_attenuation_per_sample_deg);

    std::vector<float> bpsk_modulated_wpo_real = poddsp::projection::takeProjection(bpsk_modulated_wpo,
                                                                                    poddsp::projection::type_of_projection::real_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo_fixed = poddsp::complexPLL(bpsk_modulated_wpo, freq);

    std::vector<float> bpsk_modulated_wpo_fixed_real = poddsp::projection::takeProjection(bpsk_modulated_wpo_fixed,
                                                                                          poddsp::projection::type_of_projection::real_projection);

//    std::vector<float> bpsk_phase_function = poddsp::complexSignalPhaseDependence(bpsk_modulated);


//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Действительная часть)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(complex_carrier_imag, "Несущий сигнал (Мнимая часть)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(PodDSP::complexSignalPhaseDependence(complex_carrier), "Изменение фазы в несущем сигнале");

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_phase_function, "Изменение фазы в BPSK сигнале");
//
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_real, "BPSK модулированный сигнал с ошибкой");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_fixed_real, "BPSK модулированный сигнал с исправленной ошибкой");

    std::cout << std::norm(poddsp::complexSequenceCorrelation(bpsk_modulated, bpsk_modulated_wpo_fixed));
}

TEST(transform, FFT){

    auto count_of_samples = 2000;
    auto freq = 10.0f;

    poddsp::s_sig_t impulse = poddsp::MeanderGen(freq, count_of_samples, 0, true);
//    poddsp::c_sig_t impulse_comp  = poddsp::complexSin(freq, count_of_samples, 0);
//    poddsp::s_sig_t impulseMod = poddsp::projection::takeProjection(poddsp::complexSin(freq/4, count_of_samples, 0));
//    impulse_comp = poddsp::complexMagModulator(impulse_comp, impulseMod);
//    poddsp::s_sig_t impulse = poddsp::projection::takeProjection(impulse_comp);


    poddsp::c_sig_t FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);
    FFT_analysis_result = poddsp::forwardFFT(poddsp::quadro_cast(impulse));


    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(impulse, "Меандр");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(FFT_analysis_result, "Спектр сигнала");

    FFT_analysis_result = poddsp::backwardFFT(FFT_analysis_result);

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(FFT_analysis_result, "Опять сигнал");


}

TEST(transform, specturm_plots){

    auto freq = 4.0f;
    auto count_of_samples = 2000;

    poddsp::c_sig_t Sine = poddsp::complexSin(freq, count_of_samples);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(Sine, "Исходный сигнал");

    auto Rotated_sine = poddsp::forwardFFT(Sine);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(poddsp::projection::takeProjection(Rotated_sine, poddsp::projection::imaginary_projection),
                                                              "Повёрнутый сигнал");
}


TEST(transform, hilbert_transform){
    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    poddsp::s_sig_t Sine = poddsp::projection::takeProjection(poddsp::complexSin(freq, count_of_samples, 0));
    poddsp::s_sig_t Hilberted = poddsp::transformHilbert(Sine);

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(Sine, "Source sine");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(Hilberted, "Hilbert transformed sine");

}

TEST(transform, quadro_cast){

    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    poddsp::c_sig_t Sine = poddsp::complexSin(freq, count_of_samples, 0);
    poddsp::s_sig_t pSine = poddsp::projection::takeProjection(Sine);

    poddsp::c_sig_t CSine = poddsp::quadro_cast(pSine);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(pSine, "Source sine complex");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(Sine, "Source sine projection");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(CSine, "Complex source sine");
}

TEST(transform, fftw_speed_test){
    constexpr auto count_of_samples = 8192;
    constexpr auto freq = 16.0f;
    constexpr auto count_of_processing = 1000000;

    poddsp::c_sig_t FFT_analysis_result;
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

    poddsp::s_sig_t noise = poddsp::AWGN_generator(1000);

    poddsp::c_sig_t  complex_noise = poddsp::quadro_cast(poddsp::AWGN_generator(1000));

    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(complex_noise, "АГБШ");
}

TEST(complex_functions, CFO_search){
    using namespace poddsp;
    using namespace std;

    FILE *file;

    file=fopen("../2685_f.iq", "rb");
    if(fseek(file, 0, SEEK_END) < 0){
        perror("LOH");
    }
    auto sz = ftell(file)/8;
    c_sig_t arr, parr;
    arr.reserve(sz);
    parr.reserve(sz/10000);
    rewind(file);
    float ar[2];


    for(int i = 0; i < sz; i++) {
        fread(ar, 4, 2, file);
        arr.emplace_back(std::complex<float>(ar[0], ar[1]));
    }
    fclose(file);

    for(int i = 0; i < sz/10000; i++){
        parr.emplace_back(arr.at(i));
    }



//        PlotConstructor::drawPlot(PlotConstructor::makeProjection(parr, PlotConstructor::real_projection));

    s_sig_t phase_dependence;
    s_sig_t diff_phase_dependence;
    s_sig_t second_diff_phase_dependence;
    c_sig_t res_arr;

    phase_dependence = complexSignalPhaseDependence(arr);
    phase_dependence = poddsp::phaseDependenceLining(phase_dependence);
//        PlotConstructor::drawPlot(phase_dependence, "фаза в сигнале");

    diff_phase_dependence = differentiation(phase_dependence);
//        PlotConstructor::drawPlot(diff_phase_dependence, "частота в сигнале");

    second_diff_phase_dependence = differentiation(diff_phase_dependence);
        PlotConstructor::drawPlot(second_diff_phase_dependence, "рост частоты в сигнале");

    std::cout << signalMedValue(second_diff_phase_dependence) << std::endl;
}

TEST(complex_functions, do_DC_offset_count) {

    constexpr auto count_of_samples = 4096;
    constexpr auto freq = 8.0f;

    constexpr auto dc_offset = 2.0f;

    poddsp::s_sig_t sig = poddsp::projection::takeProjection(poddsp::complexSin(freq, count_of_samples));
    sig = poddsp::signalShelf(sig, dc_offset);

    ASSERT_EQ(poddsp::signalMedValue(sig), dc_offset);
}

TEST(complex_functions, heterodyne) {
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;

    auto move_freq = 4.0f;

    std::vector<std::complex<float>> complex_carrier = poddsp::complexSin(freq, count_of_samples, 0);
    std::vector<float> complex_carrier_real = poddsp::projection::takeProjection(complex_carrier,
                                                                                      poddsp::projection::type_of_projection::real_projection);
    std::vector<float> mag_modulation = poddsp::projection::takeProjection(
            poddsp::complexSin(info_freq, info_count_of_samples, -90));
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = poddsp::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_real = poddsp::projection::takeProjection(modulated_carrier, poddsp::projection::type_of_projection::imaginary_projection);
    std::vector<float> modulated_carrier_imag = poddsp::projection::takeProjection(modulated_carrier);

//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Проекция действительной части)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier, "Амплитудно модулированный сигнал");
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
//    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    poddsp::c_sig_t het_sig = poddsp::Heterodyne(move_freq, modulated_carrier, false);
    std::vector<float> het_sig_real = poddsp::projection::takeProjection(het_sig, poddsp::projection::type_of_projection::imaginary_projection);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(het_sig_real, "Гетеродинированный Сигнал (Проекция действительной части)");

    auto demo_het_sig = poddsp::complexMagDemodulator(het_sig);
    if constexpr (DO_PLOTS) poddsp::PlotConstructor::drawPlot(demo_het_sig, "Снятый с демодулятора сигнал от гетеродинированного сигнала");

    std::cout
    << std::norm(poddsp::complexSequenceCorrelation(poddsp::quadro_cast(demo_het_sig), poddsp::quadro_cast(mag_modulation)))
    << std::endl;
}



TEST(other, cutoff_compensation) {
    using namespace poddsp;
    #define PlotConstructor if constexpr (DO_PLOTS) PlotConstructor

    constexpr auto count_of_samples = 10000;
    constexpr auto freq = 10.0f;
    constexpr auto SigMag = 1.9f;
    constexpr auto NsMag = 0.8f;
    constexpr auto CutLvl = 1.5f;

    c_sig_t complex_signal = complexSin(freq, count_of_samples);
    complex_signal = amplifier(complex_signal, SigMag);
    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Original Sin");

    auto noise = AWGN_generator(count_of_samples);
    noise = amplifier(noise, NsMag);
    PlotConstructor::drawPlot(noise, "Noise");

    for(int i =0; auto &e : complex_signal) {
        auto n_s = noise[i];
        e.real(e.real() + n_s);
        e.imag(e.imag() + n_s);
        i++;
    }

    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Noised Sin");

    complex_signal = cutoff(complex_signal, CutLvl);
    auto max_taken_val = signalMaxValue(projection::takeProjection(complex_signal));
    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Cut Sin");

    s_sig_t mags;
    for(auto &e : complex_signal) {
        mags.emplace_back(complexVectorMagnitude(e));
    }
    mags = smoother(mags);
    auto eps = 1.0f/(20.0f * signalMaxValue(mags) * dispersion(mags)); /*eps = 0.01f;*/
    auto mags_interm = s_sig_t();
    for (float & mag : mags) {
        if (mag < (max_taken_val + 10.0f*eps)) {
            continue;
        }
        mags_interm.emplace_back(mag);
    } mags = mags_interm;
    std::cout << "Eps: " << eps << std::endl;
    std::cout << "Magnitude is: " << findModeWithEps(mags, eps) << std::endl;
    std::cout << "SNR: " << 20 * logf(SigMag / NsMag) << std::endl;
    std::cout << "Mag cut percent: " << ((SigMag - CutLvl) / SigMag)*100 << "%" << std::endl;
    PlotConstructor::drawPlot(mags, "Samples Mags");
}