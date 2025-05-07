#pragma once
#include "nlohmann/json.hpp"
#include "Box.h"
#include <complex>
class BSpline
{
    public:
        BSpline() = delete;
        explicit BSpline(const nlohmann::json& input_file)
        : n_bspline{input_file.at("BSpline").at("n_bspline")}
        , order   {input_file.at("BSpline").at("order")}
        , R0      {input_file.at("BSpline").at("R0_r")}
        , eta     {input_file.at("BSpline").at("eta_r")}
        {
            degree = order - 1;
        }

        void buildLinearKnots(const Box& box);
        void buildComplexKnots();
        void buildR0();

        std::complex<double> ecs_x(double x) const;
        std::complex<double> ecs_w(double x, double w) const;
        std::complex<double> B(int i, std::complex<double> x) const;
        std::complex<double> dB(int i, std::complex<double> x) const;

        void dumpTo(const Box& box, const std::string& directory, int rank);

    private:
        int n_bspline{};
        int order{};
        int degree{};
        double R0{};
        double eta{};

        std::vector<double> knots{};
        std::vector<std::complex<double>> complex_knots{};
};
