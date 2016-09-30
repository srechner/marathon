/*
 * CudaInitFinalize.cpp
 *
 * Created on: Mar 23, 2016
 * Author: Steffen Rechner <steffen.rechner@informatik.uni-halle.de>
 *
 * This file is part of the marathon software.
 *
 * Copyright (c) 2016, Steffen Rechner
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is furnished
 * to do so, subject to the following conditions:
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
 * THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */

#include <cublas_v2.h>
#include <cublasXt.h>
#include <cuda_runtime.h>

namespace marathon {

	namespace cuda {

		// a handle for using cublas library
		cublasHandle_t cublasHandle;
		cublasXtHandle_t cublasXtHandle;
        bool cuda_init = false;

		extern bool cudaInit() {

			// check of cuda capable gpu is available
			int deviceCount = 0;
			cudaError_t error_id = cudaGetDeviceCount(&deviceCount);

			if (error_id != cudaSuccess || deviceCount == 0)
				return false;

            cuda_init = true;

			cublasCreate_v2(&cublasHandle);

			cublasXtCreate(&cublasXtHandle);

			// determine the number of gpu's
			int numDevices;
			cudaGetDeviceCount(&numDevices);

			// use all GPU's
			int devices[numDevices];
			for(int i=0; i<numDevices; i++)
				devices[i] = i;

			cublasXtDeviceSelect(cublasXtHandle, numDevices, devices);

			return true;
		}

		extern void cudaFinalize() {

            if(cuda_init) {
                cublasDestroy_v2(cublasHandle);
                cublasXtDestroy(cublasXtHandle);
                cuda_init = false;
            }
		}

	}
}


