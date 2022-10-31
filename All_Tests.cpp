#include <gtest/gtest.h>
#include <ctime>

#include <vssdsp.h>

#define RANDOM_NUMBER 1 + rand()%10
#define PlotConstructor  if (DO_PLOTS) PlotConstructor

const auto DO_PLOTS = false;
using namespace vssdsp;


TEST(complex_functions, complex_correlation_calculation){

    int sequence_length = 10;
    srand(time(nullptr));

     c_sig_t original_sequence;
     c_sig_t incoming_sequence;

    for(int i = 0; i < sequence_length; i++){
        original_sequence.emplace_back(RANDOM_NUMBER); original_sequence.back().imag(RANDOM_NUMBER);
        incoming_sequence.emplace_back(RANDOM_NUMBER); incoming_sequence.back().imag(RANDOM_NUMBER);
    }

    std::cout << complexSequenceCorrelation(original_sequence, original_sequence) << std::endl;
}

TEST(complex_functions, generating_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 10000;
    PlotConstructor::drawPlot(complexSin(freq, count_of_samples), "Sine " + std::to_string(freq));
}

TEST(complex_functions, resampling_complex_signal){

    auto freq = 10.0f;
    auto count_of_samples = 1000;
    auto new_count_of_samples = 2500;

    c_sig_t original_long_seq = complexSin(freq, count_of_samples);
    c_sig_t original_short_seq = complexSin(freq, new_count_of_samples);
    c_sig_t resampled_short_seq;

    complexSignalResampler(original_long_seq, resampled_short_seq, new_count_of_samples);


/// Draw 3M plots

//     PlotConstructor::drawPlot(original_long_seq,
//                                      (std::to_string(count_of_samples) + " samples"));
//     PlotConstructor::drawPlot(original_short_seq,
//                                      (std::to_string(new_count_of_samples) + " samples"));
//     PlotConstructor::drawPlot(resampled_short_seq,
//                                      (std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Real plots

     s_sig_t original_long_seq_Re = projection::takeProjection(original_long_seq);
    PlotConstructor::drawPlot(original_long_seq_Re, ("RE " + std::to_string(count_of_samples) + " samples"));

     s_sig_t original_short_seq_Re = projection::takeProjection(original_short_seq);
    PlotConstructor::drawPlot(original_short_seq_Re,
                                      ("RE " + std::to_string(new_count_of_samples) + " samples"));

     s_sig_t resampled_short_seq_Re = projection::takeProjection(resampled_short_seq);
    PlotConstructor::drawPlot(resampled_short_seq_Re,
                                      ("RE " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));


/// Draw 2M Imaginary plots

//     s_sig_t original_long_seq_Im =
//            projection::takeProjection(original_long_seq,
//                                                    projection::type_of_projection::imaginary_projection);
//     PlotConstructor::drawPlot(original_long_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " samples"));
//
//     s_sig_t original_short_seq_Im =
//            projection::takeProjection(original_short_seq,
//                                                    projection::type_of_projection::imaginary_projection);
//     PlotConstructor::drawPlot(original_short_seq_Im,
//                                      ("IM " + std::to_string(new_count_of_samples) + " samples"));
//
//     s_sig_t resampled_short_seq_Im =
//            projection::takeProjection(resampled_short_seq,
//                                                    projection::type_of_projection::imaginary_projection);
//     PlotConstructor::drawPlot(resampled_short_seq_Im,
//                                      ("IM " + std::to_string(count_of_samples) + " to " + std::to_string(new_count_of_samples) + " samples"));

}

TEST(complex_functions, phase_calculator){

    auto freq = 5.0f;
    auto count_of_samples = 100;
    auto count_of_periods = 4;
    c_sig_t sine;
    c_sig_t cos;

    for (int i = 0; i < 362; i++) {
        sine = complexSin(freq, count_of_samples);
        cos = complexSin(freq, count_of_samples, i);
        auto phase = complexPhaseCalculating(sine, cos);
        std::cout << phase << std::endl;
    }
}

TEST(complex_functions, phase_rotation){

    auto freq = 5.0f;
    auto count_of_samples = 1000;
    auto phase_offset_in_time = 0.0001f;
    auto phase_common_offset = 0.0f;
    c_sig_t sine;
    c_sig_t rotated_sine;
    c_sig_t rotated_with_offset_in_time_sine;
    c_sig_t disoffseted_sine;
    sine = complexSin(freq, count_of_samples);
    rotated_sine = complexPhaseChanger(sine, phase_common_offset);
    rotated_with_offset_in_time_sine = complexPhaseChanger(sine, phase_common_offset, phase_offset_in_time);
    disoffseted_sine = complexPhaseChanger(rotated_with_offset_in_time_sine, -phase_common_offset, -phase_offset_in_time);

    PlotConstructor::drawPlot(sine, "OriginalSine");
    PlotConstructor::drawPlot(rotated_sine, "RotatedSine");
    PlotConstructor::drawPlot(rotated_with_offset_in_time_sine, "RotatedwithOffSine");
    PlotConstructor::drawPlot(disoffseted_sine, "DisoffsetedSine");
}

TEST(complex_functions, complexPLL) {

    auto freq = 10.0f;
    auto count_of_samples = 2000;
    auto phase_attenuation_per_sample_deg = -0.003f;

     c_sig_t sine = complexSin(freq, count_of_samples);
//     PlotConstructor::drawPlot(sine, "Эталон");

     c_sig_t sine_wpo = complexPhaseChanger(sine, 0.0f, phase_attenuation_per_sample_deg);
//     PlotConstructor::drawPlot(sine_wpo, "Эталон с ошибкой");

     c_sig_t sine_wpo_fixed = complexPLL(sine_wpo, freq);
//     PlotConstructor::drawPlot(sine_wpo_fixed, "Работа ФАПЧ с ошибочным сигналом");

     s_sig_t sine_real = projection::takeProjection(sine,
                                                                           projection::type_of_projection::real_projection);
     PlotConstructor::drawPlot(sine_real, "Эталон (Проекция действительной части)");

     s_sig_t sine_wpo_real = projection::takeProjection(sine_wpo,
                                                                               projection::type_of_projection::real_projection);
     PlotConstructor::drawPlot(sine_wpo_real, "Эталон с ошибкой (Проекция действительной части)");

     s_sig_t sine_wpo_fixed_real = projection::takeProjection(sine_wpo_fixed,
                                                                                     projection::type_of_projection::real_projection);
     PlotConstructor::drawPlot(sine_wpo_fixed_real, "Эталон с ошибкой, исправленной ФАПЧ (Проекция действительной части)");

    std::cout << std::norm(complexSequenceCorrelation(sine, sine_wpo_fixed));
}

TEST(complex_functions, complexPLL_with_modulated_signal){
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;
    auto phase_attenuation_per_sample_deg = 0.001f;

     c_sig_t complex_carrier = complexSin(freq, count_of_samples, 0);
     s_sig_t complex_carrier_real = projection::takeProjection(complex_carrier,
                                                                                      projection::type_of_projection::real_projection);
     s_sig_t mag_modulation = projection::takeProjection(
            complexSin(info_freq, info_count_of_samples));
     PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

     c_sig_t modulated_carrier = complexMagModulator(complex_carrier, mag_modulation, 0.5f);
     s_sig_t modulated_carrier_imag = projection::takeProjection(modulated_carrier, projection::type_of_projection::imaginary_projection);
     s_sig_t modulated_carrier_real = projection::takeProjection(modulated_carrier);

     PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
     PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

     c_sig_t modulated_carrier_wpo = complexPhaseChanger(modulated_carrier, 0.0f,
                                                                                         phase_attenuation_per_sample_deg);
     s_sig_t modulated_carrier_wpo_real = projection::takeProjection(modulated_carrier_wpo);
     PlotConstructor::drawPlot(modulated_carrier_wpo_real, "Амплитудно модулированный сигнал с ошибкой (Проекция действительной части)");
//
     c_sig_t modulated_carrier_wpo_fixed = complexPLL(modulated_carrier_wpo, freq);
     s_sig_t modulated_carrier_wpo_fixed_real = projection::takeProjection(modulated_carrier_wpo_fixed);
     PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_real, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция действительной части)");
     s_sig_t modulated_carrier_wpo_fixed_imag = projection::takeProjection(modulated_carrier_wpo_fixed, projection::imaginary_projection);
     PlotConstructor::drawPlot(modulated_carrier_wpo_fixed_imag, "Амплитудно модулированный сигнал с ошибкой, исправленной ФАПЧ(Проекция мнимой части)");


     s_sig_t mag_demodulation_result = complexMagDemodulator(modulated_carrier_wpo_fixed);
     PlotConstructor::drawPlot(mag_demodulation_result, "Информационный сигнал из демодулятора");

    std::cout << "Correlation " << std::norm(complexSequenceCorrelation(modulated_carrier_wpo_fixed, modulated_carrier)) << std::endl;
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

     c_sig_t bpsk_modulated = complexBPSKModulator(info_signal, 100);
     s_sig_t bpsk_modulated_real = projection::takeProjection(bpsk_modulated);
     s_sig_t bpsk_modulated_imag = projection::takeProjection(bpsk_modulated,
                                                                                     projection::type_of_projection::imaginary_projection);

     PlotConstructor::drawPlot(bpsk_modulated, "BPSK модулированный сигнал");
//     PlotConstructor::drawPlot(complexSignalPhaseDependence(bpsk_modulated), "Зависимость фазы");
     PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
     PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");

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

//     c_sig_t complex_carrier = complexSin(freq, count_of_samples, -90);
//     s_sig_t complex_carrier_real = projection::takeProjection(complex_carrier,
//                                                                                      projection::type_of_projection::real_projection);
//
//     s_sig_t complex_carrier_imag = projection::takeProjection(complex_carrier,
//                                                                                      projection::type_of_projection::imaginary_projection);

     c_sig_t bpsk_modulated = complexBPSKModulator(info_signal, 200);
     s_sig_t bpsk_modulated_real = projection::takeProjection(bpsk_modulated,
                                                                                projection::type_of_projection::real_projection);

     s_sig_t bpsk_modulated_imag = projection::takeProjection(bpsk_modulated,
                                                                                projection::type_of_projection::imaginary_projection);

     c_sig_t bpsk_modulated_wpo = complexPhaseChanger(bpsk_modulated, 0.0f,
                                                                                      phase_attenuation_per_sample_deg);

     s_sig_t bpsk_modulated_wpo_real = projection::takeProjection(bpsk_modulated_wpo,
                                                                                    projection::type_of_projection::real_projection);

     c_sig_t bpsk_modulated_wpo_fixed = complexPLL(bpsk_modulated_wpo, freq);

     s_sig_t bpsk_modulated_wpo_fixed_real = projection::takeProjection(bpsk_modulated_wpo_fixed,
                                                                                          projection::type_of_projection::real_projection);

//     s_sig_t bpsk_phase_function = complexSignalPhaseDependence(bpsk_modulated);


//     PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Действительная часть)");
//     PlotConstructor::drawPlot(complex_carrier_imag, "Несущий сигнал (Мнимая часть)");
//     PlotConstructor::drawPlot(complexSignalPhaseDependence(complex_carrier), "Изменение фазы в несущем сигнале");

     PlotConstructor::drawPlot(bpsk_modulated_real, "BPSK модулированный сигнал (Действительная часть)");
//     PlotConstructor::drawPlot(bpsk_modulated_imag, "BPSK модулированный сигнал (Мнимая часть)");
//     PlotConstructor::drawPlot(bpsk_phase_function, "Изменение фазы в BPSK сигнале");
//
     PlotConstructor::drawPlot(bpsk_modulated_wpo_real, "BPSK модулированный сигнал с ошибкой");
     PlotConstructor::drawPlot(bpsk_modulated_wpo_fixed_real, "BPSK модулированный сигнал с исправленной ошибкой");

    std::cout << std::norm(complexSequenceCorrelation(bpsk_modulated, bpsk_modulated_wpo_fixed));
}

TEST(transform, FFT){

    auto count_of_samples = 2000;
    auto freq = 10.0f;

    s_sig_t impulse = MeanderGen(freq, count_of_samples, 0, true);
//    c_sig_t impulse_comp  = complexSin(freq, count_of_samples, 0);
//    s_sig_t impulseMod = projection::takeProjection(complexSin(freq/4, count_of_samples, 0));
//    impulse_comp = complexMagModulator(impulse_comp, impulseMod);
//    s_sig_t impulse = projection::takeProjection(impulse_comp);


    c_sig_t FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);
    FFT_analysis_result = forwardFFT(quadro_cast(impulse));


     PlotConstructor::drawPlot(impulse, "Меандр");
     PlotConstructor::drawPlot(FFT_analysis_result, "Спектр сигнала");

    FFT_analysis_result = backwardFFT(FFT_analysis_result);

     PlotConstructor::drawPlot(FFT_analysis_result, "Опять сигнал");


}

TEST(transform, specturm_plots){

    auto freq = 4.0f;
    auto count_of_samples = 2000;

    c_sig_t Sine = complexSin(freq, count_of_samples);
     PlotConstructor::drawPlot(Sine, "Исходный сигнал");

    auto Rotated_sine = forwardFFT(Sine);
     PlotConstructor::drawPlot(projection::takeProjection(Rotated_sine, projection::imaginary_projection),
                                                              "Повёрнутый сигнал");
}


TEST(transform, hilbert_transform){
    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    s_sig_t Sine = projection::takeProjection(complexSin(freq, count_of_samples, 0));
    s_sig_t Hilberted = transformHilbert(Sine);

     PlotConstructor::drawPlot(Sine, "Source sine");
     PlotConstructor::drawPlot(Hilberted, "Hilbert transformed sine");

}

TEST(transform, quadro_cast){

    const auto freq = 4.0f;
    const auto count_of_samples = 4000;

    c_sig_t Sine = complexSin(freq, count_of_samples, 0);
    s_sig_t pSine = projection::takeProjection(Sine);

    c_sig_t CSine = quadro_cast(pSine);
     PlotConstructor::drawPlot(pSine, "Source sine complex");
     PlotConstructor::drawPlot(Sine, "Source sine projection");
     PlotConstructor::drawPlot(CSine, "Complex source sine");
}

TEST(transform, fftw_speed_test){
    constexpr auto count_of_samples = 8192;
    constexpr auto freq = 16.0f;
    constexpr auto count_of_processing = 1000000;

    c_sig_t FFT_analysis_result;
    FFT_analysis_result.reserve(count_of_samples);

    auto forward_FFT =  fftwf_plan_dft_1d(count_of_samples, (fftwf_complex *)(FFT_analysis_result.data()),
                                          (fftwf_complex *)(FFT_analysis_result.data()), FFTW_FORWARD, FFTW_MEASURE);

    for(auto e : MeanderGen(freq, count_of_samples, 0, true)){
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

    s_sig_t noise = AWGN_generator(1000);

    c_sig_t  complex_noise = quadro_cast(AWGN_generator(1000));

     PlotConstructor::drawPlot(complex_noise, "АГБШ");
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
    phase_dependence = phaseDependenceLining(phase_dependence);
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

    s_sig_t sig = projection::takeProjection(complexSin(freq, count_of_samples));
    sig = signalShelf(sig, dc_offset);
    std::cout << dc_offset << '\t' << signalMedValue(sig) << std::endl;
}

TEST(complex_functions, heterodyne) {
    auto freq = 12.0f;
    auto info_freq = 4.0f;
    auto count_of_samples = 2000;
    auto info_count_of_samples = 1000;

    auto move_freq = 4.0f;

     c_sig_t complex_carrier = complexSin(freq, count_of_samples, 0);
     s_sig_t complex_carrier_real = projection::takeProjection(complex_carrier,
                                                                                      projection::type_of_projection::real_projection);
     s_sig_t mag_modulation = projection::takeProjection(
            complexSin(info_freq, info_count_of_samples, -90));
     PlotConstructor::drawPlot(mag_modulation, "Информационный сигнал");

     c_sig_t modulated_carrier = complexMagModulator(complex_carrier, mag_modulation, 0.5f);
     s_sig_t modulated_carrier_real = projection::takeProjection(modulated_carrier, projection::type_of_projection::imaginary_projection);
     s_sig_t modulated_carrier_imag = projection::takeProjection(modulated_carrier);

//     PlotConstructor::drawPlot(complex_carrier_real, "Несущий сигнал (Проекция действительной части)");
//     PlotConstructor::drawPlot(modulated_carrier, "Амплитудно модулированный сигнал");
     PlotConstructor::drawPlot(modulated_carrier_real, "Амплитудно модулированный сигнал (Проекция действительной части)");
//     PlotConstructor::drawPlot(modulated_carrier_imag, "Амплитудно модулированный сигнал (Проекция мнимой части)");

    c_sig_t het_sig = Heterodyne(move_freq, modulated_carrier, false);
     s_sig_t het_sig_real = projection::takeProjection(het_sig, projection::type_of_projection::imaginary_projection);
     PlotConstructor::drawPlot(het_sig_real, "Гетеродинированный Сигнал (Проекция действительной части)");

    auto demo_het_sig = complexMagDemodulator(het_sig);
     PlotConstructor::drawPlot(demo_het_sig, "Снятый с демодулятора сигнал от гетеродинированного сигнала");

    std::cout
    << std::norm(complexSequenceCorrelation(quadro_cast(demo_het_sig), quadro_cast(mag_modulation)))
    << std::endl;
}



TEST(other, cutoff_compensation) {

    constexpr auto count_of_samples = 100;
    constexpr auto freq = 1.0f;
    constexpr auto SigMag = 1.20f;
    constexpr auto NsMag = 0.2f;
    constexpr auto CutLvl = 1.00f;

    /// Signal gen
    c_sig_t complex_signal = complexSin(freq, count_of_samples);
    complex_signal = amplifier(complex_signal, SigMag);
    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Original Sin");

    /// Noise gen
    auto noise = AWGN_generator(count_of_samples);
    noise = amplifier(noise, NsMag);
    PlotConstructor::drawPlot(noise, "Noise");

    /// Noising
    for(int i =0; auto &e : complex_signal) {
        auto n_s = noise[i];
        e.real(e.real() + n_s);
        e.imag(e.imag() + n_s);
        i++;
    }

    /// Calculating Signal Original Magnitude
    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Noised Sin");
    complex_signal = cutoff(complex_signal, CutLvl);
    PlotConstructor::drawPlot(projection::takeProjection(complex_signal), "Cut Sin");
    auto MagFind = OriginalMagnitudeFind(complex_signal);
    std::cout << "Magnitude is: " << MagFind << std::endl;
    std::cout << "SNR is: " << (SigMag / NsMag) << std::endl;
    std::cout << "magnitude cut percent: " << ((SigMag - CutLvl) / SigMag) * 100 << "%" << std::endl;
    std::cout << "Error percent: " << (fabsf(SigMag - MagFind) / SigMag) * 100 << "%" << std::endl;
}