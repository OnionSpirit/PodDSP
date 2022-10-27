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

    std::cout << vssdsp::complexSequenceCorrelation(original_sequence, original_sequence) << std::endl;
}

TEST(complex_functions, generating_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 10000;
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(vssdsp::complexSin(freq, count_of_samples), "Sine " + std::to_string(freq));
}

TEST(complex_functions, resampling_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 1000;
    auto new_count_of_samples = 2500;

    std::vector<std::complex<float>> original_long_seq = vssdsp::complexSin(freq, count_of_samples);
    std::vector<std::complex<float>> original_short_seq = vssdsp::complexSin(freq, new_count_of_samples);
    std::vector<std::complex<float>> resampled_short_seq;

    vssdsp::complexSignalResampler(original_long_seq, resampled_short_seq, new_count_of_samples);


/// Draw 3M plots

//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_long_seq,
//                                      (std::to_string(count_of_samples) + " samples"));
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_short_seq,
//                                      (std::to_string(new_count_of_samples) + " samples"));
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(resampled_short_seq,
//                                      (std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Real plots

    std::vector<float> original_long_seq_Re = vssdsp::projection::takeProjection(original_long_seq);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_long_seq_Re, ("RE " + std::to_string(count_of_samples) + " samples"));

    std::vector<float> original_short_seq_Re = vssdsp::projection::takeProjection(original_short_seq);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_short_seq_Re,
                                      ("RE " + std::to_string(new_count_of_samples) + " samples"));

    std::vector<float> resampled_short_seq_Re = vssdsp::projection::takeProjection(resampled_short_seq);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(resampled_short_seq_Re,
                                      ("RE " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Imaginary plots

//    std::vector<float> original_long_seq_Im =
//            vssdsp::projection::takeProjection(original_long_seq,
//                                                    vssdsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_long_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " samples"));
//
//    std::vector<float> original_short_seq_Im =
//            vssdsp::projection::takeProjection(original_short_seq,
//                                                    vssdsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(original_short_seq_Im,
//                                      ("IM " + std::to_string(new_count_of_samples) + " samples"));
//
//    std::vector<float> resampled_short_seq_Im =
//            vssdsp::projection::takeProjection(resampled_short_seq,
//                                                    vssdsp::projection::type_of_projection::imaginary_projection);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(resampled_short_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));

}

TEST(complex_functions, phase_calculator){

    auto freq = 5.0f;
    auto count_of_samples = 100;
    auto count_of_periods = 4;
    std::vector<std::complex<float>> sine;
    std::vector<std::complex<float>> cos;

    for (int i = 0; i < 362; i++) {
        sine = vssdsp::complexSin(freq, count_of_samples);
        cos = vssdsp::complexSin(freq, count_of_samples, i);
        auto phase = vssdsp::complexPhaseCalculating(sine, cos);
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
    sine = vssdsp::complexSin(freq, count_of_samples);
    rotated_sine = vssdsp::complexPhaseChanger(sine, phase_common_offset);
    rotated_with_offset_in_time_sine = vssdsp::complexPhaseChanger(sine, phase_common_offset, phase_offset_in_time);
    disoffseted_sine = vssdsp::complexPhaseChanger(rotated_with_offset_in_time_sine, -phase_common_offset, -phase_offset_in_time);

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine, "OriginalSine");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(rotated_sine, "RotatedSine");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(rotated_with_offset_in_time_sine, "RotatedwithOffSine");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(disoffseted_sine, "DisoffsetedSine");
}

TEST(complex_functions, complexPLL) {

    auto freq = 10.0f;
    auto count_of_samples = 2000;
    auto phase_attenuation_per_sample_deg = -0.003f;

    std::vector<std::complex<float>> sine = vssdsp::complexSin(freq, count_of_samples);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine, "Эталон");

    std::vector<std::complex<float>> sine_wpo = vssdsp::complexPhaseChanger(sine, 0.0f, phase_attenuation_per_sample_deg);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine_wpo, "Эталон с ошибкой");

    std::vector<std::complex<float>> sine_wpo_fixed = vssdsp::complexPLL(sine_wpo, freq);
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine_wpo_fixed, "Работа ФАПЧ с ошибочным сигналом");

    std::vector<float> sine_real = vssdsp::projection::takeProjection(sine,
                                                                           vssdsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine_real, "Эталон (Проекция действительной части)");

    std::vector<float> sine_wpo_real = vssdsp::projection::takeProjection(sine_wpo,
                                                                               vssdsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine_wpo_real, "Эталон с ошибкой (Проекция действительной части)");

    std::vector<float> sine_wpo_fixed_real = vssdsp::projection::takeProjection(sine_wpo_fixed,
                                                                                     vssdsp::projection::type_of_projection::real_projection);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(sine_wpo_fixed_real, "Эталон с ошибкой, исправленной ФАПЧ (Проекция действительной части)");

    std::cout << std::norm(vssdsp::complexSequenceCorrelation(sine, sine_wpo_fixed));
}

TEST(complex_functions, complexPLL_with_modulated_signal){
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = 0.001f;

    std::vector<std::complex<float>> complex_carrier = vssdsp::complexSin(freq, count_of_samples, 0);
    std::vector<float> complex_carrier_real = vssdsp::projection::takeProjection(complex_carrier,
                                                                                      vssdsp::projection::type_of_projection::real_projection);
    std::vector<float> mag_modulation = vssdsp::projection::takeProjection(
            vssdsp::complexSin(info_freq, info_count_of_samples));
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = vssdsp::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_imag = vssdsp::projection::takeProjection(modulated_carrier, vssdsp::projection::type_of_projection::imaginary_projection);
    std::vector<float> modulated_carrier_real = vssdsp::projection::takeProjection(modulated_carrier);

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    std::vector<std::complex<float>> modulated_carrier_wpo = vssdsp::complexPhaseChanger(modulated_carrier, 0.0f,
                                                                                         phase_attenuation_per_sample_deg);
    std::vector<float> modulated_carrier_wpo_real = vssdsp::projection::takeProjection(modulated_carrier_wpo);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_wpo_real, "Амплитудно модулированный сигнал с ошибкой (Проекция действительной части)");
//
    std::vector<std::complex<float>> modulated_carrier_wpo_fixed = vssdsp::complexPLL(modulated_carrier_wpo, freq);
    std::vector<float> modulated_carrier_wpo_fixed_real = vssdsp::projection::takeProjection(modulated_carrier_wpo_fixed);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_real, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция действительной части)");
    std::vector<float> modulated_carrier_wpo_fixed_imag = vssdsp::projection::takeProjection(modulated_carrier_wpo_fixed, vssdsp::projection::imaginary_projection);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_imag, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция мнимой части)");


    std::vector<float> mag_demodulation_result = vssdsp::complexMagDemodulator(modulated_carrier_wpo_fixed);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(mag_demodulation_result, "Информационный сигнал из демодулятора");

    std::cout << "Correlation " << std::norm(vssdsp::complexSequenceCorrelation(modulated_carrier_wpo_fixed, modulated_carrier)) << std::endl;
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

    std::vector<std::complex<float>> bpsk_modulated = vssdsp::complexBPSKModulator(info_signal, 100);
    std::vector<float> bpsk_modulated_real = vssdsp::projection::takeProjection(bpsk_modulated);
    std::vector<float> bpsk_modulated_imag = vssdsp::projection::takeProjection(bpsk_modulated,
                                                                                     vssdsp::projection::type_of_projection::imaginary_projection);

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated, "BPSK модулированный сигнал");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(vssdsp::complexSignalPhaseDependence(bpsk_modulated), "Зависимость фазы");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");

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

//    std::vector<std::complex<float>> complex_carrier = vssdsp::complexSin(freq, count_of_samples, -90);
//    std::vector<float> complex_carrier_real = vssdsp::projection::takeProjection(complex_carrier,
//                                                                                      vssdsp::projection::type_of_projection::real_projection);
//
//    std::vector<float> complex_carrier_imag = vssdsp::projection::takeProjection(complex_carrier,
//                                                                                      vssdsp::projection::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated = vssdsp::complexBPSKModulator(info_signal, 200);
    std::vector<float> bpsk_modulated_real = vssdsp::projection::takeProjection(bpsk_modulated,
                                                                                vssdsp::projection::type_of_projection::real_projection);

    std::vector<float> bpsk_modulated_imag = vssdsp::projection::takeProjection(bpsk_modulated,
                                                                                vssdsp::projection::type_of_projection::imaginary_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo = vssdsp::complexPhaseChanger(bpsk_modulated, 0.0f,
                                                                                      phase_attenuation_per_sample_deg);

    std::vector<float> bpsk_modulated_wpo_real = vssdsp::projection::takeProjection(bpsk_modulated_wpo,
                                                                                    vssdsp::projection::type_of_projection::real_projection);

    std::vector<std::complex<float>> bpsk_modulated_wpo_fixed = vssdsp::complexPLL(bpsk_modulated_wpo, freq);

    std::vector<float> bpsk_modulated_wpo_fixed_real = vssdsp::projection::takeProjection(bpsk_modulated_wpo_fixed,
                                                                                          vssdsp::projection::type_of_projection::real_projection);

//    std::vector<float> bpsk_phase_function = vssdsp::complexSignalPhaseDependence(bpsk_modulated);


//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Действительная часть)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(complex_carrier_imag, "Несущий сигнал (Мнимая часть)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(vssdsp::complexSignalPhaseDependence(complex_carrier), "Изменение фазы в несущем сигнале");

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_phase_function, "Изменение фазы в BPSK сигнале");
//
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_real, "BPSK модулированный сигнал с ошибкой");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(bpsk_modulated_wpo_fixed_real, "BPSK модулированный сигнал с исправленной ошибкой");

    std::cout << std::norm(vssdsp::complexSequenceCorrelation(bpsk_modulated, bpsk_modulated_wpo_fixed));
}

TEST(transform, FFT){

    auto count_of_samples = 2000;
    auto freq = 10.0f;

    vssdsp::s_sig_t impulse = vssdsp::MeanderGen(freq, count_of_samples, 0, true);
//    vssdsp::c_sig_t impulse_comp  = vssdsp::complexSin(freq, count_of_samples, 0);
//    vssdsp::s_sig_t impulseMod = vssdsp::projection::takeProjection(vssdsp::complexSin(freq/4, count_of_samples, 0));
//    impulse_comp = vssdsp::complexMagModulator(impulse_comp, impulseMod);
//    vssdsp::s_sig_t impulse = vssdsp::projection::takeProjection(impulse_comp);


    vssdsp::c_sig_t FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);
    FFT_analysis_result = vssdsp::forwardFFT(vssdsp::quadro_cast(impulse));


    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(impulse, "Меандр");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(FFT_analysis_result, "Спектр сигнала");

    FFT_analysis_result = vssdsp::backwardFFT(FFT_analysis_result);

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(FFT_analysis_result, "Опять сигнал");


}

TEST(transform, specturm_plots){

    auto freq = 4.0f;
    auto count_of_samples = 2000;

    vssdsp::c_sig_t Sine = vssdsp::complexSin(freq, count_of_samples);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(Sine, "Исходный сигнал");

    auto Rotated_sine = vssdsp::forwardFFT(Sine);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(vssdsp::projection::takeProjection(Rotated_sine, vssdsp::projection::imaginary_projection),
                                                              "Повёрнутый сигнал");
}


TEST(transform, hilbert_transform){
    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    vssdsp::s_sig_t Sine = vssdsp::projection::takeProjection(vssdsp::complexSin(freq, count_of_samples, 0));
    vssdsp::s_sig_t Hilberted = vssdsp::transformHilbert(Sine);

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(Sine, "Source sine");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(Hilberted, "Hilbert transformed sine");

}

TEST(transform, quadro_cast){

    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    vssdsp::c_sig_t Sine = vssdsp::complexSin(freq, count_of_samples, 0);
    vssdsp::s_sig_t pSine = vssdsp::projection::takeProjection(Sine);

    vssdsp::c_sig_t CSine = vssdsp::quadro_cast(pSine);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(pSine, "Source sine complex");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(Sine, "Source sine projection");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(CSine, "Complex source sine");
}

TEST(transform, fftw_speed_test){
    constexpr auto count_of_samples = 8192;
    constexpr auto freq = 16.0f;
    constexpr auto count_of_processing = 1000000;

    vssdsp::c_sig_t FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);

    auto forward_FFT =  fftwf_plan_dft_1d(count_of_samples, (fftwf_complex *)(FFT_analysis_result.data()),
                                          (fftwf_complex *)(FFT_analysis_result.data()), FFTW_FORWARD, FFTW_MEASURE);

    for(auto e : vssdsp::MeanderGen(freq, count_of_samples, 0, true)){
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

    vssdsp::s_sig_t noise = vssdsp::AWGN_generator(1000);

    vssdsp::c_sig_t  complex_noise = vssdsp::quadro_cast(vssdsp::AWGN_generator(1000));

    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(complex_noise, "АГБШ");
}

TEST(complex_functions, CFO_search){
    using namespace vssdsp;
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
    phase_dependence = vssdsp::phaseDependenceLining(phase_dependence);
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

    vssdsp::s_sig_t sig = vssdsp::projection::takeProjection(vssdsp::complexSin(freq, count_of_samples));
    sig = vssdsp::signalShelf(sig, dc_offset);

    ASSERT_EQ(vssdsp::signalMedValue(sig), dc_offset);
}

TEST(complex_functions, heterodyne) {
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;

    auto move_freq = 4.0f;

    std::vector<std::complex<float>> complex_carrier = vssdsp::complexSin(freq, count_of_samples, 0);
    std::vector<float> complex_carrier_real = vssdsp::projection::takeProjection(complex_carrier,
                                                                                      vssdsp::projection::type_of_projection::real_projection);
    std::vector<float> mag_modulation = vssdsp::projection::takeProjection(
            vssdsp::complexSin(info_freq, info_count_of_samples, -90));
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

    std::vector<std::complex<float>> modulated_carrier = vssdsp::complexMagModulator(complex_carrier, mag_modulation, 0.5f);
    std::vector<float> modulated_carrier_real = vssdsp::projection::takeProjection(modulated_carrier, vssdsp::projection::type_of_projection::imaginary_projection);
    std::vector<float> modulated_carrier_imag = vssdsp::projection::takeProjection(modulated_carrier);

//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Проекция действительной части)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier, "Амплитудно модулированный сигнал");
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
//    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    vssdsp::c_sig_t het_sig = vssdsp::Heterodyne(move_freq, modulated_carrier, false);
    std::vector<float> het_sig_real = vssdsp::projection::takeProjection(het_sig, vssdsp::projection::type_of_projection::imaginary_projection);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(het_sig_real, "Гетеродинированный Сигнал (Проекция действительной части)");

    auto demo_het_sig = vssdsp::complexMagDemodulator(het_sig);
    if constexpr (DO_PLOTS) vssdsp::PlotConstructor::drawPlot(demo_het_sig, "Снятый с демодулятора сигнал от гетеродинированного сигнала");

    std::cout
    << std::norm(vssdsp::complexSequenceCorrelation(vssdsp::quadro_cast(demo_het_sig), vssdsp::quadro_cast(mag_modulation)))
    << std::endl;
}



TEST(other, cutoff_compensation) {
    using namespace vssdsp;
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