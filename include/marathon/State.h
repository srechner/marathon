/*
 * State.h
 *
 *  Created on: Mar 3, 2016
 *      Author: rechner
 */

#ifndef INCLUDE_MARATHON_STATE_H_
#define INCLUDE_MARATHON_STATE_H_

#include <cstdlib>
#include <string>
#include <iostream>

namespace marathon {

/**
 * Abstract Base Class for States.
 */
class State {

public:

	virtual ~State() {

	}

	/**
	 * Virtual Hash Function for State Type.
	 */
	virtual size_t hash_value() const = 0;

	/**
	 * Compare this and s by structural properties.
	 *
	 * If this<s : return -1.
	 * If this==s: return 0.
	 * If this>s : return 1.
	 */
	virtual int compare_to(const State* s) const = 0;

	/**
	 * Return a string representation of the state.
	 */
	virtual std::string to_string() const = 0;

	/**
	 * To output into streams.
	 */
	friend inline std::ostream& operator<<(std::ostream& out, const State& s) {
		out << s.to_string();
		return out;
	}

	/**
	 * To output into streams.
	 */
	friend inline std::ostream& operator<<(std::ostream& out, const State* s) {
		out << s->to_string();
		return out;
	}

	/**
	 * Wrapper Class for the use in std::unordered_maps.
	 */
	class Hash {
	public:
		size_t operator()(State * x) const {
			const size_t res = x->hash_value();
			return res;
		}
	};

	/*
	 * Wrapper Class for the use in std::unordered_maps.
	 */
	class Equal {
	public:
		bool operator()(State * x1, State * x2) const {
			const bool res = x1->compare_to(x2) == 0;
			return res;
		}
	};

};

}

#endif /* INCLUDE_MARATHON_STATE_H_ */
