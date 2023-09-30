/* =============================================================================

 Copyright (C) 2009-2021 Valerii Sukhorukov. All Rights Reserved.

 Permission is hereby granted, free of charge, to any person obtaining a copy
 of this software and associated documentation files (the "Software"), to deal
 in the Software without restriction, including without limitation the rights
 to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 copies of the Software, and to permit persons to whom the Software is
 furnished to do so, subject to the following conditions:

 The above copyright notice and this permission notice shall be included in all
 copies or substantial portions of the Software.

 THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 SOFTWARE.

================================================================================
*/

/**
 * \file with_cuda.h
 * \brief Contains class Cuda.
 * \author Valerii Sukhorukov
 */

#ifndef UTILS_RANDOM_WITH_CUDA_H
#define UTILS_RANDOM_WITH_CUDA_H

//#include <cuda_runtime.h>
//#include <curand.h>

#include "../common/misc.h"
#include "core.h"

/// Pseugo-random number generation.
namespace utils::random {

/// \brief Random number factory based on Cuda rng library.
/// \tparam real Floating point type.
template<std::floating_point real> 
class Cuda
    : public Core<real> {

    // Ensure that the template parameter is a floating type:
    static_assert(
        std::is_floating_point<real>::value,
        "Class Cuda can only be instantiated with floating-point types."
    );

    using Core<real>::buffersize;

    real* rU01;    ///< Buffer array for storing random numbers (host).
    real* d_Rand;  ///< Buffer array for storing random numbers (device).
    
    int rU01_ind;   ///< Index of the current random number in \a rU01.
                                                
    std::mt19937      gCPU; ///< Random number generator using CPU.
    curandGenerator_t gGPU; ///< Random number generator using GPU.

public:

//    /// \brief Default constructor.
//    Cuda() = default;

    /// \brief Constructor setting the seed uncoupled from run index.
    /// \param seed Random number generator seed.
    /// \param runName Human-readable run index.
    /// \param msgr Output message processor.
    explicit Cuda(
        unsigned seed,
        const std::string& runName,
        Msgr& msgr);

    /// \brief Constructor setting the seed depending on run index.
    /// \param seedFname Name of the file contining seeds.
    /// \param runIndex Run index.
    /// \param msgr Output message processor.
    explicit Cuda(
        unsigned runIndex,
        Msgr& msgr);

    // The rule of five is triggered by the destructor, the defaults suffice:
    Cuda(const Cuda&) = delete;             ///< copy constructor
    Cuda& operator=(const Cuda&) = delete;  ///< copy assignment
    Cuda(Cuda&&) = delete;                  ///< move constructor
    Cuda& operator=(Cuda&&) = delete;       ///< move assignment
    ~Cuda();                                ///< destructor

    /// Initialize CUDA rng machinery.
    void initialize_CUDA_rng();
    
    /// A pseudo-random number with uniform distribution over [0,1).
    real r01u();

    /// \brief A pseudo-random unsigned int from the range [0, max-1].
    /// \param max Max boundary of the sampled range.
    uint uniformInt0(const uint max);

private:

    /// Populate the buffer array with a new butch of random numbers over (0,1].
    void prepare_uniform_real01();

    /// Check CUDA errors.
    void checkCudaErrors(const curandStatus_t& e);

    /// Check CUDA errors.
    void checkCudaErrors(const cudaError_t& e);
};


// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template<std::floating_point real>
Cuda<real>::
Cuda(
    const unsigned seed,
    const std::string& runName,
    Msgr& msgr)
    : Core<real> {seed, runName, msgr}
{
    gCPU.seed(this->seed);

    initialize_CUDA_rng();

    rU01_ind = -1;
    prepare_uniform_real01();
}


template<std::floating_point real>
Cuda<real>::
Cuda(
    const unsigned runInd,
    Msgr& msgr)
    : Core<real> {runInd, msgr}
{
    
    gCPU.seed(this->seed);
    
    initialize_CUDA_rng();
    
    rU01_ind = -1;
    prepare_uniform_real01();
}


template<std::floating_point real>
Cuda<real>::
~Cuda()
{
    checkCudaErrors(curandDestroyGenerator(gGPU));
    checkCudaErrors(cudaFree(d_Rand));
    free(rU01);
}


template<std::floating_point real> 
void Cuda<real>::
initialize_CUDA_rng()
{
     checkCudaErrors(curandCreateGenerator(&gGPU, CURAND_RNG_PSEUDO_MTGP32));
     checkCudaErrors(curandSetPseudoRandomGeneratorSeed(gGPU, this->seed));
     checkCudaErrors(cudaMalloc((void **)&d_Rand, buffersize * sizeof(real)));
     rU01 = (real *)malloc(buffersize * sizeof(real));
}


// Generates real random numbers with uniform distribution over (0,1]   (!!!)
template<> 
void Cuda<float>::
prepare_uniform_real01()
{
    checkCudaErrors(curandGenerateUniform(gGPU, (float *) d_Rand, buffersize));
    // synchronizes the GPU threads before copy
    checkCudaErrors(cudaMemcpy(rU01, d_Rand, buffersize * sizeof(float),
                    cudaMemcpyDeviceToHost));
}


template<> 
void Cuda<double>::
prepare_uniform_real01()
{
    checkCudaErrors(curandGenerateUniformDouble(gGPU, (double *) d_Rand,
                                                buffersize));
    // synchronizes the GPU threads before copy
    checkCudaErrors(cudaMemcpy(rU01, d_Rand, buffersize * sizeof(double),
                               cudaMemcpyDeviceToHost));
}


// returns int in the range [0,max-1]
template<std::floating_point real> 
uint Cuda<real>::
uniformInt0(
    const uint& max
)
{            
    auto ir {static_cast<uint>(r01u() * max)};
    while (ir >= max)
        ir = static_cast<uint>(r01u() * max);

    return ir;
}


// returns a random number with uniform distribution over [0,1)
template<std::floating_point real>
real Cuda<real>::r01u()
{    
    if (++rU01_ind == buffersize) {
        prepare_uniform_real01();
        rU01_ind = 0;
    }
    return rU01[rU01_ind];
}


template<std::floating_point real>
void Cuda<real>::
checkCudaErrors(
    const curandStatus_t& err
)
{
    if (err != CURAND_STATUS_SUCCESS)
        throw common::Exception(
            "CURAND error: " + std::to_string(err), this->msgr);
}


template<std::floating_point real>
void Cuda<real>::
checkCudaErrors(
    const cudaError_t& err
)
{
    if (err != cudaSuccess)
        throw common::Exception(
            "CUDA error: " + std::to_string(err), this->msgr);
}


}  // namespace utils::random

#endif // UTILS_RANDOM_WITH_CUDA_H
