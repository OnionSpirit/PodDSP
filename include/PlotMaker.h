#pragma once

namespace poddsp {

    class PlotConstructor {

/// Рисование графика в трёхмерной системе координат с комплексной осью
        static void drawComplexPlot(const std::vector<std::complex<float>> &,
                                    const std::string & = "NoTitle",
                                    int = 1)
                                    noexcept;

    public:
        enum type_of_projection {
            real_projection = 0,
            imaginary_projection = 1
        };

/// Рисование графика (Для трёхмерного графика, входной массив должен быть комплексным).
/// Параметры:
/// 1) Массив данных,
/// 2) Заголовок окна графика, по умолчанию NoTitle,
/// 3) Величина отступа от точки к точке по оси Х
        template<typename T>
        static void drawPlot(const std::vector<T> &data,
                             const std::string &title = "NoTitle",
                             int step = 1)
                             noexcept {
            std::ofstream datafile;
            if (typeid(T) == typeid(std::complex<float>) ||
                typeid(T) == typeid(std::complex<double>) ||
                typeid(T) == typeid(std::complex<int>) ||
                typeid(T) == typeid(std::complex<char>)) {
                std::vector<std::complex<float>> complex_data;
                for (auto e : data) {
                    complex_data.emplace_back(e);
                }
                PlotConstructor::drawComplexPlot(complex_data, title, step);
                return;
            }
            datafile.open("plotData.dat");
            for (int i = 0; auto e: data) {
                datafile << i << "\t" << e << std::endl;
                i += step;
            }
            auto gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set grid xtics ytics \n");
            fprintf(gp, "set title '%s' \n", title.c_str());
            fprintf(gp, "plot '%s' using 1:2 with lines \n", "plotData.dat");
            pclose(gp);

            if (remove("plotData.dat")) std::filesystem::filesystem_error(ERROR_PLOT"Temporary file removing failure",
                                                                          std::error_code());
        }

/// Создание проекции трёхмерного графика на одну из плоскостей. Праметры:
/// 1) Массив комплексных данных,
/// 2) Тип проекции согласно PodDSP::PlotConstructor::type_of_projection, по умолчанию real_projection.
        static std::vector<float> makeProjection(const std::vector<std::complex<float>> &,
                                                 const type_of_projection & = type_of_projection::real_projection);
    };
}
