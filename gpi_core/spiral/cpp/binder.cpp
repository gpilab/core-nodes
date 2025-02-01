#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include "WPGen.hpp"
#include "WPCompose.hpp"
#include "SpiralOrder.hpp"
#include <memory>

namespace py = pybind11;

using namespace wp;

PYBIND11_MODULE(whirledpeas, m) {
    ////////////////////////////////////////////////////////////////////////////
    py::class_<WPGen>(m, "WPGen")
        .def(py::init<double, double, double, double, double, double, double, double, double, bool, double>(),
        py::arg("fov") = 0.01 * 24.0, // mm
        py::arg("res") = 0.01 * 1.0, // mm
        py::arg("reqTau"),
        py::arg("sMax") = 150.0, // mT/m/sec
        py::arg("gMax") = 40.0, // mT/m
        py::arg("omegaMax") = 2.0 * M_PI,
        py::arg("gradRast") = 0.001 * 4,
        py::arg("dwell") = 0.001 * 2,
        py::arg("gamma") = 42.57,
        py::arg("addOuterRing") = false,
        py::arg("alphaVD") = 0,
        "WPGen class constructor that takes Tau as input")

        .def(py::init<double, double, int, double, double, double, double, double, double, bool, double>(),
            py::arg("fov") = 0.01 * 24.0, // mm
            py::arg("res") = 0.01 * 1.0, // mm
            py::arg("number_arms") = 16,
            py::arg("slew_max") = 150.0, //mT/m/sec
            py::arg("grad_max") = 40.0, //mT/m
            py::arg("omega_max") = 2.0 * M_PI,
            py::arg("grast") = 0.001 * 4,
            py::arg("dwell") = 0.001 * 2,
            py::arg("gamma") = 42.57,
            py::arg("addOuterRing") = false,
            py::arg("alphaVD") = 0,
            "WP Spiral class constructor that takes n arms as input")

        .def("GetSpiralReadPts", &WPGen::GetSpiralReadPts)
        .def("GetSpiralRampPts", &WPGen::GetSpiralRampPts)
        .def("GetTotalSpiralPts", &WPGen::GetTotalSpiralPts)
        .def("GetM0M1Pts", &WPGen::GetM0M1Pts)
        .def("GetSlewMax", &WPGen::GetSlewMax)
        .def("GetGradMax", &WPGen::GetGradMax)
        .def("GetNyqNumberArms", &WPGen::GetNyqNumberArms)
        .def("GetArcDur", &WPGen::GetArcDur)
        .def("GetFreqDur", &WPGen::GetFreqDur)
        .def("GetSlewDur", &WPGen::GetSlewDur)
        .def("GetGradDur", &WPGen::GetGradDur)
        .def("GetRampDur", &WPGen::GetRampDur)
        .def("GetActualTau", &WPGen::GetActualTau)
        .def("GetOmegaMax", &WPGen::GetOmegaMax)
        .def("GetSnrFactor", &WPGen::GetSnrFactor)
        .def("GetGridMtxSize", &WPGen::GetGridMtxSize)
        .def("GetOutRingDur", &WPGen::GetOutRingDur)
        .def("GetKspacePts", &WPGen::GetKspacePts)
        .def("ComputeBaseSpiral", [](WPGen& self)
        {
            std::vector<double> spiralX;
            std::vector<double> spiralY;

            self.ComputeBaseSpiral(spiralX, spiralY);

            const int totalPts = spiralX.size();
            py::array_t<double> spiral({2, totalPts});

            auto m0m1_ptr = spiral.mutable_unchecked<2>();
            
            for (int i = 0; i < totalPts; i++) {
                m0m1_ptr(0, i) = spiralX[i]; // Assign m0m1X[i] to the first dimension at index i
                m0m1_ptr(1, i) = spiralY[i]; // Assign m0m1Y[i] to the second dimension at index i
            }
            return spiral;
        })
        .def("ComputeRampDown", [](WPGen& self)
        {
            std::vector<double> rampX;
            std::vector<double> rampY;

            self.ComputeRampDown(rampX, rampY);

            const int totalPts = rampY.size();
            py::array_t<double> spiral({2, totalPts});

            auto m0m1_ptr = spiral.mutable_unchecked<2>();
            
            for (int i = 0; i < totalPts; i++) {
                m0m1_ptr(0, i) = rampX[i]; // Assign m0m1X[i] to the first dimension at index i
                m0m1_ptr(1, i) = rampY[i]; // Assign m0m1Y[i] to the second dimension at index i
            }
            return spiral;
        })
        .def("ComputeM0M1Trapezoids", [](WPGen& self)
        {
            std::vector<double> m0m1X;
            std::vector<double> m0m1Y;

            self.ComputeM0M1Trapezoids(m0m1X, m0m1Y);

            const int totalPts = m0m1X.size();
            py::array_t<double> m0m1({2, totalPts});

            auto m0m1_ptr = m0m1.mutable_unchecked<2>();
            
            for (int i = 0; i < totalPts; i++) {
                m0m1_ptr(0, i) = m0m1X[i]; // Assign m0m1X[i] to the first dimension at index i
                m0m1_ptr(1, i) = m0m1Y[i]; // Assign m0m1Y[i] to the second dimension at index i
            }
            return m0m1;
        })
        .def("ComputeKSPnSDC", [](WPGen& self) {

            std::vector<double> kspX;
            std::vector<double> kspY;
            std::vector<double> sdc_;

            self.ComputeKSPnSDC(kspX, kspY, sdc_);

            const int totalPts = kspX.size();
            py::array_t<double> ksp({2, totalPts});
            py::array_t<double> sdc(totalPts);

            auto ksp_ptr = ksp.mutable_unchecked<2>();
            for (int i = 0; i < totalPts; i++) {
                ksp_ptr(0, i) = kspX[i]; // Assign kspX[i] to the first dimension at index i
                ksp_ptr(1, i) = kspY[i]; // Assign kspY[i] to the second dimension at index i
            }

            auto sdc_ptr = sdc.mutable_unchecked<1>();
            // Copy sdc data to sdc numpy array
            for (int i = 0; i < totalPts; i++) {
                sdc_ptr(i) = sdc_[i]; // Assign sdc[i] to the first dimension at index i
            }

            return std::make_tuple(ksp, sdc);
        });

    ////////////////////////////////////////////////////////////////////////////
    // Define the WPCompose class

    ////////////////////////////////////////////////////////////////////////////

    py::enum_<WPCompose::SpiralDirect>(m, "SpiralDirect")
        .value("SPIRAL_OUT", WPCompose::SpiralDirect::SPIRAL_OUT)
        .value("SPIRAL_IN", WPCompose::SpiralDirect::SPIRAL_IN)
        .value("SPIRAL_INOUT", WPCompose::SpiralDirect::SPIRAL_INOUT)
        .value("SPIRAL_OUTIN", WPCompose::SpiralDirect::SPIRAL_OUTIN)
        .export_values();

    py::class_<WPCompose>(m, "WPCompose")
        .def(py::init<WPGen*, int, WPCompose::SpiralDirect>(),
            py::arg("wp1"),
            py::arg("numEchoes"),
            py::arg("spiralType"),
            "WPCompose class constructor for same Tau spiral")
        .def(py::init<WPGen*, WPGen*, int>(),
            py::arg("wp1"),
            py::arg("wp2"),
            py::arg("centerNullPts"),
            "WPCompose class constructor for dual Tau spiral")
        .def("composeGradients", [](WPCompose& self) 
        {
            // Compose gradients
            std::vector<double> gradX;
            std::vector<double> gradY;
            self.composeGradients(gradX, gradY);

            double* gX = gradX.data();
            const int totalPts = gradX.size();

            // Use the NumPy API to create a 2D array with the appropriate shape
            py::array_t<double> gradients({totalPts, 2});

            // Obtain a pointer to the underlying data
            auto gradients_ptr = gradients.mutable_unchecked<2>(); 

            // Copy gradX and gradY data to gradients numpy array along the two dimensions
            for (int i = 0; i < totalPts; i++) {
                gradients_ptr(i, 0) = gradX[i]; // Assign gradX[i] to the first dimension at index i
                gradients_ptr(i, 1) = gradY[i]; // Assign gradY[i] to the second dimension at index i
            }

            return gradients;

        })
        .def("composeKSPnSDC", [](WPCompose& self) 
        {
            // Compose K-space and SDC
            std::vector<double> kspX;
            std::vector<double> kspY;
            std::vector<double> sdc_;
            std::vector<int> kspLengths_;
            self.composeKSPnSDC(kspX,kspY, sdc_, false, kspLengths_);

            const int nEchoes = self.getNEchoes();
            const int totalKspPts = kspX.size() / nEchoes; // Size of the first dimension of ksp

            // Use the NumPy API to create a 3D array with the appropriate shape
            py::array_t<double> ksp_array({nEchoes, totalKspPts, 2});
            // Obtain a pointer to the underlying data
            auto ksp_ptr = ksp_array.mutable_unchecked<3>(); 

            // Copy ksp data to ksp numpy array along the three dimensions
            for (int echo = 0; echo < nEchoes; echo++) {
                for (int i = 0; i < totalKspPts; i++) {
                    ksp_ptr(echo, i, 0) = kspX[echo * totalKspPts + i]; 
                    ksp_ptr(echo, i, 1) = kspY[echo * totalKspPts + i]; 
                }
            }

            // Use the NumPy API to create a 2D array with the appropriate shape
            py::array_t<double> sdc(sdc_.size());
            // Obtain a pointer to the underlying data
            auto sdc_ptr = sdc.mutable_unchecked<1>();

            // Copy sdc data to sdc numpy array
            for (size_t i = 0; i < sdc_.size(); i++) {
                sdc_ptr(i) = sdc_[i]; // Assign sdc_[echo * totalKspPts + i] to the first dimension at index i
            }

            return std::make_tuple(ksp_array, sdc);

        })
        .def("composeTimeMap", [](WPCompose& self) 
        {
            // Compose time map with input and output
            dVector2D timeMapIn, timeMapOut;
            self.composeTimeMap(timeMapIn, timeMapOut);

            const int gridMtxSize = self.getGridMtxSize();
            auto convertToNumpyArray = [](const dVector2D& vec, int size) {
                if (vec.empty() || vec[0].empty()) {
                    return py::array_t<double>(py::none());
                }
                py::array_t<double> array({size, size});
                auto array_ptr = array.mutable_unchecked<2>();
                for (int i = 0; i < size; i++) {
                    for (int j = 0; j < size; j++) {
                        array_ptr(j, i) = vec[i][j];
                    }
                }
                return array;
            };

            py::array_t<double> timeMapIn_array = convertToNumpyArray(timeMapIn, gridMtxSize);
            py::array_t<double> timeMapOut_array = convertToNumpyArray(timeMapOut, gridMtxSize);

            return std::make_tuple(timeMapIn_array, timeMapOut_array);

        })
        .def("getTimeEcho1", &WPCompose::getTimeEcho1)
        .def("getNyqArms", &WPCompose::getNyqArms)
        .def("getGradMax", &WPCompose::getGradMax)
        .def("getDeltaTE", &WPCompose::getDeltaTE);

        py::enum_<SpiralArmOrder::OrderType>(m, "OrderType")
            .value("LINEAR", SpiralArmOrder::OrderType::LINEAR)
            .value("SKIP", SpiralArmOrder::OrderType::SKIPP)
            .value("TWO_WAY", SpiralArmOrder::OrderType::TWO_WAY)
            .value("MIXED", SpiralArmOrder::OrderType::MIXED)
            .value("GOLDEN", SpiralArmOrder::OrderType::GOLDEN)
            .export_values();

        py::class_<SpiralArmOrder>(m, "SpiralArmOrder")
            .def(py::init<int, int, float, SpiralArmOrder::OrderType, bool, bool>(),
            py::arg("number_arms"),
            py::arg("number_echo"),
            py::arg("angular_coverage"),
            py::arg("order_type"),
            py::arg("in_coherent_mech"),
            py::arg("debug"),
            "SpiralArmOrder class constructor")
            .def("composeAngles", [](SpiralArmOrder& self){

            dVector2D anglesArray;
            self.composeAngles(anglesArray);

            // Accessing arg1 and arg2 from the class constructor
            int nArms = anglesArray.size();
            int nEcho = anglesArray[0].size();

            py::array_t<double> angles({nEcho, nArms});
            // Obtain a pointer to the underlying data
            auto angles_ptr = angles.mutable_unchecked<2>();

            // Copy angles data to angles numpy array along the two dimensions
            for (int echo = 0; echo < nEcho; echo++) {
                for (int arm = 0; arm < nArms; arm++) {
                angles_ptr(echo, arm) = anglesArray[arm][echo];
                }
            }

            return angles;

            });

    m.doc() = "pybind11 example plugin"; // optional module docstring
}
