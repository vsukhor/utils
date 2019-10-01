/* ==============================================================================

 Copyright (C) 2009-2019, Valerii Sukhorukov, <vsukhorukov@yahoo.com>

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

============================================================================== */

/**
* \file with_cuda.h
* Contains class Cuda.
* \author Valerii Sukhorukov
*/

#ifndef UTILS_RANDOM_WITH_CUDA_H
#define UTILS_RANDOM_WITH_CUDA_H

#include <cuda_runtime.h>
#include <curand.h>

#include "common/misc.h"
#include "core.h"

/// Library-wide.
namespace Utils {
/// \brief Pseugo-random number generation.
namespace Random {

using namespace Common;
using namespace Arrays;

/// \brief Random number factory based on Cuda rng library.
/// \tparam realT Floating point type.
template <typename realT> 
class Cuda
	: public Core<realT> {

	// Ensure that the template parameter is a floating type
	static_assert(std::is_floating_point<realT>::value,
				  "Class Core can only be instantiated with floating point types");

	using Core<realT>::buffersize;

	realT* rU01;			///< Buffer array for storing random numbers (host).
	realT* d_Rand;			///< Buffer array for storing random numbers (device).
	
	int	   rU01_ind;		///< Index of the current random number in \a rU01.
												
	std::mt19937		gCPU;		///< Random number generator using CPU.
	curandGenerator_t	gGPU;		///< Random number generator using GPU.

public:

//	/// \brief Default constructor.
//	Cuda() = default;

	/// \brief Constructor.
	/// \param seedFname Name of the file contining seeds.
	/// \param ii Run index.
	/// \param msgr Output message processor.
	explicit Cuda(
		const std::string& seedFname,
		const szt ii,
		Msgr& msgr);

	/// \brief Constructor.
	/// \param seed Random number generator seed.
	/// \param runName Human-readable run index.
	/// \param msgr Output message processor.
	explicit Cuda(
		const int seed,
		const std::string& runName,
		Msgr& msgr);

	// The rule of five is triggered by the destructor, the defaults suffice
    Cuda(const Cuda&) = delete;        		///< copy constructor
    Cuda& operator=(const Cuda&) = delete; 	///< copy assignment
    Cuda(Cuda&&) = delete;             		///< move constructor
    Cuda& operator=(Cuda&&) = delete;  		///< move assignment
	~Cuda();								///< destructor

	/// \brief Initialize CUDA rng machinery.
	void initialize_CUDA_rng();
	
	/// A pseudo-random number with uniform distribution over [0,1).
	realT r01u();

	/// \brief A pseudo-random unsigned int from the range [0, max-1].
	/// \param max Max boundary of the sampled range.
	uint uniformInt0(const uint max);

private:

	/// \brief Populate the buffer array \a rU01 with a new butch of random numbers over (0,1] (!!!).
	void prepare_uniform_real01();

	/// \brief Check CUDA errors.
	void checkCudaErrors(const curandStatus_t& e);
	/// \brief Check CUDA errors.
	void checkCudaErrors(const cudaError_t& e);
};

// IMPLEMENTATION ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

template <typename realT> 
Cuda<realT>::
Cuda( 	const std::string& seedFname,
		const szt ii,
		Msgr& msgr )
	: RandCore<realT> {msgr, seedFname, ii}
{
	
	gCPU.seed(this->seed);
	
    initialize_CUDA_rng();
	
	rU01_ind = -1;
	prepare_uniform_real01();
}

template <typename realT> 
Cuda<realT>::
Cuda( const int seed,
		const std::string& runName,
		Msgr& msgr )
	: Core<realT> {msgr, seed, runName}
{
	gCPU.seed(this->seed);
	
    initialize_CUDA_rng();
	
	rU01_ind = -1;
	prepare_uniform_real01();
}

template <typename realT> 
Cuda<realT>::
~Cuda()
{
    checkCudaErrors(curandDestroyGenerator(gGPU));
    checkCudaErrors(cudaFree(d_Rand));
    free(rU01);
}

template <typename realT> 
void Cuda<realT>::
initialize_CUDA_rng()
{
     checkCudaErrors(curandCreateGenerator(&gGPU, CURAND_RNG_PSEUDO_MTGP32));
     checkCudaErrors(curandSetPseudoRandomGeneratorSeed(gGPU, this->seed));
     checkCudaErrors(cudaMalloc((void **)&d_Rand, buffersize * sizeof(realT)));
     rU01 = (realT *)malloc(buffersize * sizeof(realT));
}
// Generates real random numbers with uniform distribution over (0,1]   (!!!)
template <> 
void Cuda<float>::
prepare_uniform_real01()
{
	checkCudaErrors(curandGenerateUniform(gGPU, (float *) d_Rand, buffersize));
	// synchronizes the GPU threads before copy
    checkCudaErrors(cudaMemcpy(rU01, d_Rand, buffersize * sizeof(float), cudaMemcpyDeviceToHost));
}

template <> 
void Cuda<double>::
prepare_uniform_real01()
{
	checkCudaErrors(curandGenerateUniformDouble( gGPU, (double *) d_Rand, buffersize));
	// synchronizes the GPU threads before copy
    checkCudaErrors(cudaMemcpy(rU01, d_Rand, buffersize * sizeof(double), cudaMemcpyDeviceToHost));
}

// returns int in the range [0,max-1]
template <typename realT> 
uint Cuda<realT>::
uniformInt0( const uint& max )
{			
	auto ir {static_cast<uint>(r01u() * max)};
	while (ir >= max)
		ir = static_cast<uint>(r01u() * max);

	return ir;
}

// returns a random number with uniform distribution over [0,1)
template<typename realT>
realT Cuda<realT>::r01u()
{	
	if (++rU01_ind == buffersize) {
		prepare_uniform_real01();
		rU01_ind = 0;
	}
	return rU01[rU01_ind];
}
template<typename realT>
void Cuda<realT>::
checkCudaErrors( const curandStatus_t& err )
{
	if (err != CURAND_STATUS_SUCCESS)
		throw Exceptions::Simple("CURAND error: " + STR(err), this->msgr);
}
template<typename realT>
void Cuda<realT>::
checkCudaErrors( const cudaError_t& err )
{
	if (err != cudaSuccess)
		throw Exceptions::Simple("CUDA error: " + STR(err), this->msgr);
}

}	// namespace Random
}	// namespace Utils

#endif // UTILS_RANDOM_WITH_CUDA_H
