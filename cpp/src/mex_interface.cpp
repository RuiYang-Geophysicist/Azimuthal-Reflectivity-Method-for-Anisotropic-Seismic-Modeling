/**
 * @file mex_interface.cpp
 * @brief MATLAB MEX interface for AzRM C++ library (Traditional C API)
 *
 * Usage in MATLAB:
 *   [Rpp, Rpsv, Rpsh] = azrm_mex(freq, angles, thickness, density, stiffness, phi)
 *
 * @author Rui Yang, 2024
 */

#include "mex.h"
#include "matrix.h"
#include "azrm_core.hpp"
#include <cstring>

/**
 * @brief MEX gateway function (Traditional C API)
 *
 * [Rpp, Rpsv, Rpsh] = azrm_mex(freq, angles, thickness, density, stiffness, phi)
 */
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    // Check number of inputs
    if (nrhs < 6 || nrhs > 7) {
        mexErrMsgIdAndTxt("AZRM:nrhs",
            "Usage: [Rpp, Rpsv, Rpsh] = azrm_mex(freq, angles, thickness, density, stiffness, phi, [nthreads])");
    }

    // Check number of outputs
    if (nlhs > 3) {
        mexErrMsgIdAndTxt("AZRM:nlhs", "Too many output arguments");
    }

    // Validate inputs
    // freq (scalar)
    if (!mxIsDouble(prhs[0]) || mxIsComplex(prhs[0]) || mxGetNumberOfElements(prhs[0]) != 1) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "freq must be a real double scalar");
    }

    // angles (vector)
    if (!mxIsDouble(prhs[1]) || mxIsComplex(prhs[1])) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "angles must be a real double vector");
    }

    // thickness (vector)
    if (!mxIsDouble(prhs[2]) || mxIsComplex(prhs[2])) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "thickness must be a real double vector");
    }

    // density (vector)
    if (!mxIsDouble(prhs[3]) || mxIsComplex(prhs[3])) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "density must be a real double vector");
    }

    // stiffness (6x6xnL)
    if (!mxIsDouble(prhs[4]) || mxIsComplex(prhs[4])) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "stiffness must be a real double 3D array (6x6xnL)");
    }
    const mwSize *dims = mxGetDimensions(prhs[4]);
    mwSize ndims = mxGetNumberOfDimensions(prhs[4]);
    if (ndims < 2 || dims[0] != 6 || dims[1] != 6) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "stiffness must have dimensions 6x6xnL");
    }

    // phi (scalar)
    if (!mxIsDouble(prhs[5]) || mxIsComplex(prhs[5]) || mxGetNumberOfElements(prhs[5]) != 1) {
        mexErrMsgIdAndTxt("AZRM:invalidInput", "phi must be a real double scalar");
    }

    // Extract inputs
    double freq = mxGetScalar(prhs[0]);

    double *anglesPtr = mxGetPr(prhs[1]);
    int nAngles = static_cast<int>(mxGetNumberOfElements(prhs[1]));

    double *thicknessPtr = mxGetPr(prhs[2]);
    int nThickness = static_cast<int>(mxGetNumberOfElements(prhs[2]));

    double *densityPtr = mxGetPr(prhs[3]);
    int nDensity = static_cast<int>(mxGetNumberOfElements(prhs[3]));

    double *stiffnessPtr = mxGetPr(prhs[4]);
    int nLayers = (ndims >= 3) ? static_cast<int>(dims[2]) : 1;

    double phi = mxGetScalar(prhs[5]);

    // Optional: number of threads
    if (nrhs >= 7) {
        int nthreads = static_cast<int>(mxGetScalar(prhs[6]));
        azrm::azrm_set_num_threads(nthreads);
    }

    // Validate dimensions
    if (nThickness != nLayers || nDensity != nLayers) {
        mexErrMsgIdAndTxt("AZRM:dimMismatch",
            "thickness, density, and stiffness must have consistent layer count");
    }

    // Allocate output arrays
    plhs[0] = mxCreateDoubleMatrix(1, nAngles, mxCOMPLEX);
    if (nlhs >= 2) plhs[1] = mxCreateDoubleMatrix(1, nAngles, mxCOMPLEX);
    if (nlhs >= 3) plhs[2] = mxCreateDoubleMatrix(1, nAngles, mxCOMPLEX);

    // Allocate temporary arrays for C interface
    std::vector<double> Rpp_real(nAngles), Rpp_imag(nAngles);
    std::vector<double> Rpsv_real(nAngles), Rpsv_imag(nAngles);
    std::vector<double> Rpsh_real(nAngles), Rpsh_imag(nAngles);

    // Call C++ implementation
    azrm::azrm_compute_Rpp(
        freq,
        anglesPtr,
        nAngles,
        thicknessPtr,
        densityPtr,
        stiffnessPtr,
        nLayers,
        phi,
        Rpp_real.data(),
        Rpp_imag.data(),
        Rpsv_real.data(),
        Rpsv_imag.data(),
        Rpsh_real.data(),
        Rpsh_imag.data()
    );

    // Copy results to MATLAB output
    #if MX_HAS_INTERLEAVED_COMPLEX
    // R2018a+ interleaved complex
    mxComplexDouble *RppOut = mxGetComplexDoubles(plhs[0]);
    mxComplexDouble *RpsvOut = (nlhs >= 2) ? mxGetComplexDoubles(plhs[1]) : nullptr;
    mxComplexDouble *RpshOut = (nlhs >= 3) ? mxGetComplexDoubles(plhs[2]) : nullptr;

    for (int i = 0; i < nAngles; ++i) {
        RppOut[i].real = Rpp_real[i];
        RppOut[i].imag = Rpp_imag[i];
        if (RpsvOut) {
            RpsvOut[i].real = Rpsv_real[i];
            RpsvOut[i].imag = Rpsv_imag[i];
        }
        if (RpshOut) {
            RpshOut[i].real = Rpsh_real[i];
            RpshOut[i].imag = Rpsh_imag[i];
        }
    }
    #else
    // Legacy separate real/imag
    double *RppReal = mxGetPr(plhs[0]);
    double *RppImag = mxGetPi(plhs[0]);
    double *RpsvReal = (nlhs >= 2) ? mxGetPr(plhs[1]) : nullptr;
    double *RpsvImag = (nlhs >= 2) ? mxGetPi(plhs[1]) : nullptr;
    double *RpshReal = (nlhs >= 3) ? mxGetPr(plhs[2]) : nullptr;
    double *RpshImag = (nlhs >= 3) ? mxGetPi(plhs[2]) : nullptr;

    for (int i = 0; i < nAngles; ++i) {
        RppReal[i] = Rpp_real[i];
        RppImag[i] = Rpp_imag[i];
        if (RpsvReal) {
            RpsvReal[i] = Rpsv_real[i];
            RpsvImag[i] = Rpsv_imag[i];
        }
        if (RpshReal) {
            RpshReal[i] = Rpsh_real[i];
            RpshImag[i] = Rpsh_imag[i];
        }
    }
    #endif
}
