#include "../include/vssdsp.h"


namespace vssdsp {

    void PlotConstructor::drawComplexPlot(const std::vector<std::complex<float>> &data,
                                          const std::string &title,
                                          int step)
    noexcept {
        std::ofstream datafile;
        datafile.open("plotData.dat");

        for (int i = 0; auto e: data) {
            datafile << i << "\t" << e.real() << "\t" << e.imag() << std::endl;
            i += step;
        }
        datafile.close();

        auto gp = popen("gnuplot -persist", "w");
        fprintf(gp, "set grid xtics ytics ztics \n");
        fprintf(gp, "set title '%s' \n", title.c_str());
        fprintf(gp, "splot '%s' using 2:3:1 with lines \n", "plotData.dat");
        pclose(gp);

        if (remove("plotData.dat")) std::cout << "Temporary file removing failure" << std::endl;
    }
}
