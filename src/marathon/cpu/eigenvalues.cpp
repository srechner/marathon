/*
 * eigenvalues.h
 *
 *  Created on: Nov 23, 2014
 *      Author: steffen
 */

#ifndef EIGENVALUES_H_
#define EIGENVALUES_H_

// system includes
#include <iostream>
#include <cstdlib>
#include <algorithm>
#include <cmath>

// Project Includes
#include "../../../include/marathon/exceptions.h"
#include "../../../include/marathon/cpu/analyzer.h"
#include "../../../include/arpack++/arlssym.h"
#include "../../../include/arpack++/arlsmat.h"

/**
 * Computes the second largest value of the spectrum of transition matrix P.
 *
 * Let the sequence of eigenvalues be (-1 < l1 <= l2 <= ... <= l_N = 1).
 * The method uses ARPACK++ library to compute the two eigenvalues with largest
 * magnitude and returns the lesser one.
 */
template<typename T>
T marathon::cpu::secondLargestEigenvalue(const StateGraph* mc) {

	// Variables
	size_t omega;			// Number of states
	size_t numArcs;			// Number of transitions
	unsigned int i;

	omega = mc->getNumStates();

	// Check trivial cases
	if (omega == 1)
		return 1;

	// count number of non-zero elements in left upper part of transition matrix
	numArcs = 0;
	for (Transition t : mc->arcs) {
		if (t.u <= t.v)
			numArcs++;
	}

	// 0. Try to alloate memory
	T *p = new T[numArcs];
	if (p == nullptr) {
		throw BAD_HOST_MALLOC_EXCEPTION;
	}
	int *col = new int[numArcs];
	if (col == nullptr) {
		delete[] p;
		throw BAD_HOST_MALLOC_EXCEPTION;
	}
	int *row = new int[omega + 1];
	if (row == nullptr) {
		delete[] p;
		delete[] col;
		throw BAD_HOST_MALLOC_EXCEPTION;
	}

	// 1. Define Matrix
	ARluSymMatrix<T> B;
	for (i = 0; i < omega + 1; i++)
		row[i] = numArcs;

	i = 0;
	for (std::vector<Transition>::const_iterator it = mc->arcs.begin();
			it != mc->arcs.end(); ++it) {
		if (it->u <= it->v) {
			// symmetrize (see theory paper)
			T x = it->p.convert_to<T>();
			Rational y = mc->getWeight(it->u) / mc->getWeight(it->v);
			x *= sqrt(y.convert_to<T>());
			p[i] = x;
			col[i] = it->v;
			if (i < row[it->u])
				row[it->u] = i;
			i++;
		}
	}

	/*std::cout << "P: ";
	 for (i = 0; i < numArcs; i++)
	 std::cout << p[i] << " ";
	 std::cout << std::endl;

	 std::cout << "row: ";
	 for (i = 0; i <= omega; i++)
	 std::cout << row[i] << " ";
	 std::cout << std::endl;

	 std::cout << "col: ";
	 for (i = 0; i < numArcs; i++)
	 std::cout << col[i] << " ";
	 std::cout << std::endl;*/

	B.DefineMatrix(omega, numArcs, p, col, row, 'L');

	// 2. Define Problem

	// Symmetric Standard Eigenvalue Problem
	// Compute the two eigenvalues with largest magnitude

	ARluSymStdEig<T> prob;
	int nev;
	T tol = 1e-10f;
	int maxit = 1000 * omega;
	char which[3];

	if (omega == 2) {
		nev = 1;
		strcpy(which, "SM");
	} else {
		nev = 2;
		strcpy(which, "LM");
	}

	int ncv = std::min(2 * nev, (int) omega);

	//std::cout << ncv << std::endl;

	prob.DefineParameters(omega, nev, &B, &ARluSymMatrix<T>::MultMv, which, ncv,
			tol, maxit);

	// 3. Solve Problem
	prob.FindEigenvalues();

	/*std::cout << "found " << prob.ConvergedEigenvalues() << std::endl;
	 for (int i = 0; i < prob.ConvergedEigenvalues(); i++) {
	 std::cout << i << ": " << prob.Eigenvalue(i) << std::endl;
	 }*/

	// free memory
	delete[] p;
	delete[] col;
	delete[] row;

	return prob.Eigenvalue(0);
}

/**
 * Export Template Specialization
 */

template float marathon::cpu::secondLargestEigenvalue<float>(const StateGraph*);
template double marathon::cpu::secondLargestEigenvalue<double>(
		const StateGraph*);

#endif /* EIGENVALUES_H_ */
