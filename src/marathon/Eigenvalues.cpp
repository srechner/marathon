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
#include <climits>

#include "../../include/marathon/Eigenvalues.h"

#include "../../include/marathon/marathon.h"
#include "../../include/arpack++/arlssym.h"
#include "../../include/arpack++/arlsmat.h"

namespace marathon {
namespace eigenvalue {

/**
 * Computes the second largest value of the spectrum of transition matrix P.
 *
 * Let the sequence of eigenvalues be (1 >= l1 >= l2 >= ... >= lN > -1).
 * The method uses ARPACK++ library to compute the two eigenvalues with
 * largest magnitute and returns the lesser one.
 */
template<typename T>
T eigenvalue(const StateGraph* sg, eigenvalue_t opt) {

	// Variables
	size_t omega;			// Number of states
	size_t numArcs;			// Number of transitions
	unsigned int i;

	omega = sg->getNumStates();

	// Check trivial cases
	if (omega == 1)
		return 1;

	// count number of non-zero elements in left upper part of transition matrix
	numArcs = 0;
	for (Transition* t : sg->getArcs()) {
		if (t->u <= t->v)
			numArcs++;
	}

	// 0. Try to alloate memory
	T *p = new T[numArcs];
	if (p == nullptr) {
		throw std::runtime_error(
				"marathon::cpu::eigenvalues::exception: bad memory alloc!");
	}
	int *col = new int[numArcs];
	if (col == nullptr) {
		delete[] p;
		throw std::runtime_error(
				"marathon::cpu::eigenvalues::exception: bad memory alloc!");
	}
	int *row = new int[omega + 1];
	if (row == nullptr) {
		delete[] p;
		delete[] col;
		throw std::runtime_error(
				"marathon::cpu::eigenvalues::exception: bad memory alloc!");
	}

	// 1. Define Matrix
	ARluSymMatrix<T> B;
	for (i = 0; i < omega + 1; i++)
		row[i] = numArcs;

	// make a copy of the arc pointers
	std::vector<Transition*> copy(sg->getArcs());

	// sort transitions ascending
	std::sort(copy.begin(), copy.end(),
			[](const Transition* a, const Transition* b) -> bool
			{
				if(a->u == b->u)
				return a->v < b->v;
				else
				return a->u < b->u;
			});
	i = 0;
	for (const Transition* t : copy) {
		if (t->u <= t->v) {
			// symmetrize (see theory paper)
			T x = t->p.convert_to<T>();
			rational y = sg->getWeight(t->u) / sg->getWeight(t->v);
			x *= sqrt(y.convert_to<T>());
			p[i] = x;
			col[i] = t->v;
			if (i < row[t->u])
				row[t->u] = i;
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
	T tol = std::numeric_limits<T>::epsilon() * 10.0;// one power more than machine precision
	int maxit = 1000 * omega;
	char which[3];

	// decide which eigenvalues to compute
	if (omega == 2) {
		nev = 1;
		strcpy(which, "SM");
	} else if (opt == _2ndLargestMagnitude) {
		nev = 2;
		strcpy(which, "LM");
	} else if (opt == _2ndLargestAlgebraic) {
		nev = 2;
		strcpy(which, "LA");
	} else {
		throw std::runtime_error(
				"marathon::cpu::eigenvalue::exception: unknown option: "
						+ std::to_string(opt));
	}

	int ncv = std::min(2 * nev, (int) omega);

	/*std::cout << "omega = " << omega << std::endl;
	std::cout << "nev   = " << nev << std::endl;
	std::cout << "ncv   = " << ncv << std::endl;
	std::cout << "maxit = " << maxit << std::endl;
	std::cout << "tol   = " << tol << std::endl;*/

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

template float eigenvalue<float>(const StateGraph*, eigenvalue_t);
template double eigenvalue<double>(const StateGraph*, eigenvalue_t);

}
}

#endif /* EIGENVALUES_H_ */
