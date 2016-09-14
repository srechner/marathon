/*
 * CudaWrapper.h
 *
 *  Created on: Mar 23, 2016
 *      Author: rechner
 */

#ifndef CUDAWRAPPER_H_
#define CUDAWRAPPER_H_

namespace marathon {

/**
 * Declare external functions that do the CUDA stuff. Mostly Wrapper around CUDA functions.
 */

    extern void myCudaMallocPitch(void **devPtr, size_t *pitch, size_t width,
                                  size_t height);

    extern void myCudaMalloc(void **devPtr, size_t size);

    extern void myCudaFree(void *devPtr);

    extern void myCudaMemcpyHostToDevice(void *dst, const void *src,
                                         size_t count);

    extern void myCudaMemset2D(void *devPtr, size_t pitch, int value,
                               size_t width, size_t height);

    extern void myCudaMemcpy2DHostToDevice(void *dst, size_t dpitch,
                                           const void *src, size_t spitch, size_t width, size_t height);

    extern void myCudaMemcpy2DDeviceToHost(void *dst, size_t dpitch,
                                           const void *src, size_t spitch, size_t width, size_t height);

    extern void myCudaMemcpy2DDeviceToDevice(void *dst, size_t dpitch,
                                             const void *src, size_t spitch, size_t width, size_t height);


}

#endif /* CUDAWRAPPER_H_ */
