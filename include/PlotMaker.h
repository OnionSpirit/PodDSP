#pragma once

namespace vssdsp {

    class PlotConstructor {

/// Рисование графика в трёхмерной системе координат с комплексной осью
        static void drawComplexPlot(const std::vector<std::complex<float>> &,
                                    const std::string & = "NoTitle",
                                    int = 1)
        noexcept;
    public:

/// Рисование графика (Для трёхмерного графика, входной массив должен быть комплексным).
/// Параметры:
/// 1) Массив данных,
/// 2) Заголовок окна графика, по умолчанию NoTitle,
/// 3) Величина отступа от точки к точке по оси Х
///ИЛИ ПРИ ИСПОЛЬЗОВАНИИ ПОЛЬЗОВАТЕЛЬСКОЙ ШКАЛЫ Х
/// 1) Массив первой величины (X),
/// 2) Массив второй величины,
/// 3) Заголовок окна графика, по умолчанию NoTitle,
/// 4) Величина отступа от точки к точке по оси Х
        template<typename T>
        static void drawPlot(const std::vector<T> &data,
                             const std::string &title = "NoTitle",
                             int step = 1) {
            std::ofstream datafile;
            if (typeid(T) == typeid(std::complex<float>) ||
                typeid(T) == typeid(std::complex<double>) ||
                typeid(T) == typeid(std::complex<int>) ||
                typeid(T) == typeid(std::complex<char>))
            {
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

            if (remove("plotData.dat")) std::cout << "Temporary file removing failure" << std::endl;
        }

        template<typename T>
        static void drawPlot(const std::vector<T> &data_y,
                             const std::vector<T> &data_x,
                             const std::string &title = "NoTitle",
                             int step = 1)
        noexcept {
            std::ofstream datafile;
            datafile.open("plotData.dat");
            const auto len = data_x.size() > data_y.size() ? data_y.size() : data_x.size();
            for (int i = 0; i < len; i++) {
                datafile << data_x[i] << "\t" << data_y[i] << std::endl;
                i += step;
            }
            datafile.close();

            auto gp = popen("gnuplot -persist", "w");
            fprintf(gp, "set grid xtics ytics \n");
            fprintf(gp, "set title '%s' \n", title.c_str());
            fprintf(gp, "plot '%s' using 1:2 with lines \n", "plotData.dat");
            pclose(gp);

            if (remove("plotData.dat")) std::cout << "Temporary file removing failure" << std::endl;
        }
    };
}
