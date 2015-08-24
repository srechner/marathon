// system includes
#include <iostream>
#include <fstream>
#include <cmath>
#include <string>
#include <climits>
#include <pthread.h>

// project includes
#include "../../../include/marathon/chains/matching/chain_JSV04.h"
#include "../../../include/marathon/chains/matching/chain_JS89.h"
#include "../../../include/marathon/chains/sequences/switch_chain_bipartite.h"
#include "../../../include/marathon/chains/sequences/switch_chain_bipartite_fast.h"

// sampling includes
#include "../../../include/marathon.h"

using namespace std;
using namespace marathon;

time_t start;

/**
 * Data Packet for threaded Computation of MC-Properties
 */
typedef struct data_package {

	StateGraph *mc;
	uint line;				// line number
	string s;				// hashstring
	double eps;				// epsilon
	size_t omega;			// number of states
	Rational pimin;			// minimal stationary propability
	int t;					// total mixing time

	// lower and upper bounds
	double lower_bound_eigen;
	double upper_bound_eigen;
	double upper_bound_canon;
	double lambda_max;

	data_package(StateGraph* mc, uint line, string s, double eps) :
			mc(mc), line(line), s(s), eps(eps), omega(mc->getNumStates()), pimin(
					0), t(0), lower_bound_eigen(
					-1), upper_bound_eigen(-1), upper_bound_canon(-1), lambda_max(
					-1) {

		if (omega > 0) {
			pimin = mc->getMinimalStationary();
		}
	}

} thread_data;

std::ostream &operator<<(std::ostream &out, data_package const &d) {

	// get current time
	time_t now = time(NULL);
	time_t dur = now - start;
	uint sec = dur % 60;
	dur /= 60;
	uint min = dur % 60;
	dur /= 60;
	uint hours = dur;

	printf("(%ih:%02im:%02is)\t", hours, min, sec);
	cout << d.line << "\t";
	cout << d.s << "\t";
	cout << d.omega << "\t";
	cout << d.pimin << "\t";

	if (d.lambda_max != -1) {
		cout << d.lambda_max << "\t";
		if (d.lambda_max < 0)
			cout << "0\t";
		else
			cout << "1\t";
	} else
		cout << "NA\tNA\t";

	if (d.t >= 0)
		out << d.t << "\t";
	else
		out << "NA" << "\t";

	if (d.lower_bound_eigen != -1)
		out << d.lower_bound_eigen << "\t";
	else
		out << "NA" << "\t";

	if (d.upper_bound_eigen != -1)
		out << d.upper_bound_eigen << "\t";
	else
		out << "NA" << "\t";

	if (d.upper_bound_canon != -1)
		out << d.upper_bound_canon << "\t";
	else
		out << "NA" << "\t";

	return out;
}

void *mixingTimeThread(void *arg) {
	data_package *tdata = (data_package *) arg;

	StateGraph* mc = tdata->mc;
	double eps = tdata->eps;

	size_t omega = mc->getNumStates();

	// for small matrices do cpu work
	if (omega < 20) {
		tdata->t = marathon::cpu::totalMixingTime<double>(mc, eps);
	} else if (omega <= 20000) {

		// try to call Device function
		tdata->t = marathon::gpu::totalMixingTime<double>(mc, eps);

		// if not succesful: try hybrid
		if (tdata->t == -1) {

			tdata->t = marathon::hybrid::totalMixingTime<double>(mc, eps);
		}

		// if not succesful: try host (fallback solution)
		if (tdata->t == -1) {
			//std::cout << "try host" << std::endl;
			tdata->t = marathon::cpu::totalMixingTime<double>(mc, eps);

		}

	}

	pthread_exit(NULL);
}

void *canonicalPathThread(void *arg) {
	data_package *tdata = (data_package *) arg;

	StateGraph* mc = tdata->mc;
	double eps = tdata->eps;
	size_t omega = tdata->omega;

	if (omega > 1 && omega <= 20000) {
		// variables
		Rational pimin = 1;		// minimal stationary propability
		double x;				// temporary variable
		Rational congestion;	// canonical path congestion

		x = -log(tdata->pimin.convert_to<double>() * eps);
		congestion = marathon::cpu::pathCongestion(mc);
		tdata->upper_bound_canon = x * congestion.convert_to<double>();
	}

	pthread_exit(NULL);
}

void *eigenThread(void *arg) {
	data_package *tdata = (data_package *) arg;

	StateGraph* mc = tdata->mc;
	double eps = tdata->eps;

	size_t omega = mc->getNumStates();

	if (omega > 1 && omega <= 20000) {

		// variables
		Rational pimin = 1;		// minimal stationary propability
		double x;				// temporary variable

		x = -log(tdata->pimin.convert_to<double>() * eps);

		double lambda = marathon::cpu::secondLargestEigenvalue<double>(mc);

		if (fabs(lambda) < 1e-15)
			lambda = 0.0;
		tdata->lambda_max = lambda;

		lambda = fabs(lambda);

		tdata->lower_bound_eigen = 0.5 * (lambda / (1.0 - lambda))
				* -log(2 * eps);
		tdata->upper_bound_eigen = x / (1.0 - lambda);

	}

	pthread_exit(NULL);
}

int main(int argc, char** argv) {

	if (argc != 4) {
		cout
				<< "usage: bounds <file> <js89|jsv04|swapBip|swapBipFast|curveball> <epsilon>"
				<< endl;
		return 1;
	}

	ifstream file(argv[1]);
	if (!file.is_open()) {
		cerr << "error while reading input file: " << argv[1] << endl;
		return 1;
	}

	double eps = 1.0;

	try {
		eps = atof(argv[3]);
	} catch (...) {
		cerr << "problems while reading epsilon: " << argv[3] << endl;
		return 1;
	}

	// start time measurement
	start = time(NULL);

	// print headline
	cout << "time\t";
	cout << "line\t";
	cout << "instance\t";
	cout << "omega\t";
	cout << "pimin\t";
	cout << "lambda\t";
	cout << "[lambda>0]\t";
	cout << "mix\t";
	cout << "leigen\t";
	cout << "ueigen\t";
	cout << "ucanon\t";
	cout << endl;

	// init library
	marathon::gpu::init();
	marathon::hybrid::init();

	// read input file
	string s = "";
	size_t line = 1;
	std::vector<string> lines;

	while (std::getline(file, s)) {

		marathon::StateGraph *sg = nullptr;

		// check which chain is selected
		if (strcmp(argv[2], "js89") == 0)
			sg = new marathon::chain::matching::MatchingChain89(s);
		else if (strcmp(argv[2], "jsv04") == 0)
			sg = new marathon::chain::matching::MatchingChain04(s);
		else if (strcmp(argv[2], "swapBip") == 0)
			sg = new marathon::chain::sequence::SwapChainBipartite(s);
		else if (strcmp(argv[2], "swapBipFast") == 0)
			sg = new marathon::chain::sequence::SwapChainBipartiteFast(s);
		//else if (strcmp(argv[1], "curveball") == 0)
		//	sg = new marathon::chain::sequence::Curveball(s);
		else {
			std::cerr << "unknown chain specifier: " << argv[1] << std::endl;
			return 1;
		}

		// construct state graph
		sg->constructStatespace();

		// compute properties of state graph
		if (sg->getNumStates() > 1) {

			pthread_t tid[3];		// for asynchronous compation of mixing time

			// fill data package with default values
			data_package tdata(sg, line, s, eps);

			// compute useful information about mc
			pthread_create(&tid[0], NULL, mixingTimeThread, (void *) &tdata);
			pthread_create(&tid[1], NULL, canonicalPathThread, (void *) &tdata);
			pthread_create(&tid[2], NULL, eigenThread, (void *) &tdata);

			// wait for completion of threads
			pthread_join(tid[0], NULL);
			pthread_join(tid[1], NULL);
			pthread_join(tid[2], NULL);

			// print
			cout << tdata << endl;
		}

		line++;
		delete sg;

	}

	marathon::gpu::finalize();
	marathon::hybrid::finalize();

	return 0;
}
