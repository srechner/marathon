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
         * Create the conjugate vector of src.
         * @param src Integer vector.
         * @param len Length of conjugated vector.
         * @return Conjugate vector of src.
         */
        inline
        std::vector<int>
        conjugate(const std::vector<int> &src, const size_t len) {
            std::vector<int> dst(len, 0);
            for (size_t i = 0; i < src.size(); i++) {
                for (size_t l = 0; l < std::min(static_cast<size_t>(src[i]), len); l++) {
                    dst[l]++;
                }
            }
            return dst;
        }

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
        bool isDominating(const std::vector<int> &p, const std::vector<int> &q) {

            if (p.size() != q.size())
                return false;

            int sum_p = 0;
            int sum_q = 0;
            for (int i = 0; i < p.size(); i++) {
                sum_p += p[i];
                sum_q += q[i];
                if (sum_p < sum_q) {
                    return false;
                }
            }
            return true;
        }

        /**
         * Does p dominate q?
         * @param p Sequence of integers.
         * @param q Sequence of integers.
         * @param len Length of sequences.
         * @return True, if p dominates q. False, otherwise.
         */
        inline
        bool isDominating(const int *p, const int *q, const int len) {
            int sum_p = 0;
            int sum_q = 0;
            for (int i = 0; i < len; i++) {
                sum_p += p[i];
                sum_q += q[i];
                if (sum_p < sum_q) {
                    return false;
                }
            }
            return true;
        }

    }
}

#endif //MARATHON_COMMON_H
