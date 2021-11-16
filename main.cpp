#include <iostream>
#include <complex>
#include <vector>
#include <ctime>

#define RANDOM_NUMBER 1 + rand()%10

std::complex<float> complexCorrelation(std::vector<std::complex<float>> const &, std::vector<std::complex<float>> const &);

int main() {

    int sequence_length = 10;
    srand(time(nullptr));

    std::vector<std::complex<float>> original_sequence;
    std::vector<std::complex<float>> incoming_sequence;

    for(int i = 0; i < sequence_length; i++){
        original_sequence.emplace_back(RANDOM_NUMBER); original_sequence.back().imag(RANDOM_NUMBER);
        incoming_sequence.emplace_back(RANDOM_NUMBER); incoming_sequence.back().imag(RANDOM_NUMBER);
    }

    std::cout << complexCorrelation(original_sequence, incoming_sequence) << std::endl;

    return 0;
}

std::complex<float> complexCorrelation(std::vector<std::complex<float>> const & original_sequence, std::vector<std::complex<float>> const & incoming_sequence){

    int sequence_size = (int)original_sequence.size();
    std::complex<float> correlation_result;
    std::complex<float> conjugated_incoming_sequence[sequence_size];


    std::complex<float> real_avg_original;
    std::complex<float> imaginary_avg_original;
    std::complex<float> dispersion_original;
    std::complex<float> centralized_sequence_original[sequence_size];

    for (std::complex<float> e : original_sequence){
        real_avg_original += e.real();
    } real_avg_original *= (1.0f / (float)sequence_size);

    float temp_original = 0;
    for (std::complex<float> e : original_sequence){
        temp_original += e.imag();
    } imaginary_avg_original.imag(temp_original); imaginary_avg_original *= (1.0f / (float)sequence_size);

    for (int i = 0; std::complex<float> e : original_sequence){
        centralized_sequence_original[i] = (e - (real_avg_original + imaginary_avg_original));
        i++;
    }

    for (std::complex<float> e : centralized_sequence_original){
        dispersion_original += norm(e);
    } dispersion_original *= (1.0f / (float)sequence_size);



    std::complex<float> real_avg_incoming;
    std::complex<float> imaginary_avg_incoming;
    std::complex<float> dispersion_incoming;
    std::complex<float> centralized_sequence_incoming[sequence_size];

    for(int i = 0; std::complex<float> e : incoming_sequence){
        conjugated_incoming_sequence[i] = conj(e);
        i++;
    }

    for (std::complex<float> e : conjugated_incoming_sequence){
        real_avg_incoming += e.real();
    } real_avg_incoming *= (1.0f / (float)sequence_size);

    float temp_incoming = 0;
    for (std::complex<float> e : conjugated_incoming_sequence){
        temp_incoming += e.imag();
    } imaginary_avg_incoming.imag(temp_incoming); imaginary_avg_incoming *= (1.0f / (float)sequence_size);

    for (int i = 0; std::complex<float> e : conjugated_incoming_sequence){
        centralized_sequence_incoming[i] = (e - (real_avg_incoming + imaginary_avg_incoming));
        i++;
    }

    for (std::complex<float> e : centralized_sequence_incoming){
        dispersion_incoming += norm(e);
    } dispersion_incoming *= (1.0f / (float)sequence_size);



    for(int i = 0; i < sequence_size; i++){
        correlation_result += centralized_sequence_incoming[i] * centralized_sequence_original[i];
    } correlation_result *= (1.0f / (float)sequence_size); correlation_result /= pow(dispersion_incoming * dispersion_original, 0.5);



    return correlation_result;
}
