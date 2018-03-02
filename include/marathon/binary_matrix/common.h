//
// Created by rechner on 22.05.17.
//

#ifndef MARATHON_COMMON_H
#define MARATHON_COMMON_H

#include <cstdlib>
#include <cstring>

namespace marathon {
    namespace binary_matrix {

        /**
         * Create the conjugate sequence to src.
         * @param dst Sequence of length len_dst.
         * @param src Sequence of length len_src.
         * @param len_dst Length of dst.
         * @param len_src Length of src.
         */
        inline
        void conjugate(int *dst, const int *src, const int len_dst, const int len_src) {
            memset(dst, 0, len_dst * sizeof(int));
            for (int i = 0; i < len_src; i++) {
                for (int l = 0; l < std::min(src[i], len_dst); l++) {
                    dst[l]++;
                }
            }
        }

        /**
         * Does p dominate q?
         * @param p Sequence of integers.
         * @param q Sequence of integers.
         * @param len Length of sequences.
         * @return True, if p dominates q. False, otherwise.
         */
        inline
        bool isDominating(const int* p, const int* q, const int len) {
            int sum_p = 0;
            int sum_q = 0;
            for(int i=0; i<len; i++) {
                sum_p += p[i];
                sum_q += q[i];
                if(sum_p < sum_q) {
                    return false;
                }
            }
            return true;
        }

    }
}

#endif //MARATHON_COMMON_H
