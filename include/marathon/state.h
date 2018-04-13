/*
 * Created on: Mar 3, 2016
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

#ifndef INCLUDE_MARATHON_STATE_H_
#define INCLUDE_MARATHON_STATE_H_

#include <cstdlib>
#include <string>
#include <iostream>
#include <memory>

namespace marathon {

/**
 * Abstract Base Class for States.
 */
    class State {

    public:

        /**
         * Create a deep copy of the state.
         */
        virtual std::unique_ptr<State> copy() const = 0;

        /**
         * Return a hash function of the state.
         */
        virtual size_t hashValue() const = 0;

        /**
         * Compare the state with another state.
         * @param s State object.
         * @return Zero, if this==s. Negative value, if this<s. Positive value if s<this.
         */
        virtual int compare(const State &s) const = 0;

        /**
         * Return a string representation of the state.
         */
        virtual std::string toString() const = 0;

        /**
         * Output into streams.
         */
        friend inline std::ostream &operator<<(std::ostream &out, const State &s) {
            out << s.toString();
            return out;
        }

        /**
         * Output into streams.
         */
        friend inline std::ostream &operator<<(std::ostream &out, const State *s) {
            out << s->toString();
            return out;
        }

        /**
         * Wrapper Class for the use in std::unordered_maps.
         */
        class Hash {
        public:
            size_t operator()(State *x) const {
                const size_t res = x->hashValue();
                return res;
            }
        };

        /*
         * Wrapper Class for the use in std::unordered_maps.
         */
        class Equal {
        public:
            bool operator()(State *x1, State *x2) const {
                const bool res = x1->compare(*x2) == 0;
                return res;
            }
        };

        /*
         * Wrapper Class for the use in std::maps.
         */
        class Less {
        public:
            bool operator()(State *x1, State *x2) const {
                const bool res = x1->compare(*x2) < 0;
                return res;
            }
        };

    };

}

#endif /* INCLUDE_MARATHON_STATE_H_ */
