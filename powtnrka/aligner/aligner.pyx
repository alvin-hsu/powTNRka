# Copyright (C) 2024   Alvin Hsu

from typing import Optional, Tuple

from cpython.mem cimport PyMem_Malloc, PyMem_Realloc, PyMem_Free

from .aligner cimport *


cdef Direction UP = 0x01
cdef Direction LEFT = 0x02
cdef Direction DIAG = 0x04
cdef float NINF = -10000.0
cdef float EPS = 1e-8
cdef Base base_mapping[256]
base_mapping[ord('A')] = A
base_mapping[ord('C')] = C
base_mapping[ord('G')] = G
base_mapping[ord('T')] = T
base_mapping[ord('M')] = M
base_mapping[ord('R')] = R
base_mapping[ord('W')] = W
base_mapping[ord('S')] = S
base_mapping[ord('Y')] = Y
base_mapping[ord('K')] = K
base_mapping[ord('B')] = B
base_mapping[ord('D')] = D
base_mapping[ord('H')] = H
base_mapping[ord('V')] = V
base_mapping[ord('N')] = N
# Dummy for repeats/heterogeneous
for i in range(0, 10):
    base_mapping[48 + i] = N
base_mapping[ord('*')] = N


cdef Base *map_bases(str s):
    """
    Takes in a Python string s and allocates a Base * that corresponds to it.
    """
    cdef Base *retval = <Base *>PyMem_Malloc((len(s) + 1)*sizeof(Base))
    if not retval: raise MemoryError()
    for i, c in enumerate(s):
        retval[i] = base_mapping[ord(c)]
    retval[len(s)] = <Base>0
    return retval


cdef inline size_t is_close(float f1, float f2):
    return <size_t>(-EPS < f1 - f2 < EPS)


cdef class PowtnrkaAligner:
    # Scoring information
    cdef float match, degen_match, mismatch, gap_open, gap_ext, gap_rep, gap_incentive, het_score
    # Reference information
    cdef readonly str ref_str
    cdef size_t ref_len
    cdef Base *reference
    cdef readonly list nick_sites
    # Incentivizing gaps at cutsites
    cdef float *gap_incentives
    # Heterogeneous regions
    cdef readonly size_t num_hets
    cdef size_t *het_idxs
    # Repeat information
    cdef readonly size_t num_repeats
    cdef size_t *repeat_idxs
    cdef size_t *repeat_lens
    cdef Base **repeats
    cdef readonly list repeat_strs
    # Data storage - do NOT touch directly!
    cdef float *_score_data
    cdef size_t max_query_len
    # Scoring matrices (pointers that index into _score_data)
    cdef float *m_arr
    cdef float *i_arr
    cdef float *d_arr
    cdef float **r_arrs
    # Direction matrices for tracebacks
    cdef Direction *tb_arr  # Match matrix traceback

    def __init__(self, str reference, int num_repeats, list repeats, list nick_sites,
                 float match = 5.0, float degen_match = 1.0, float mismatch = -3.0,
                 float gap_open = -10.0, float gap_ext = -5.0, float gap_rep = -3.0,
                 float gap_incentive = 0.0, float het_score = 1.0, size_t init_query_len = 63):
        # Reference information
        self.ref_str = reference
        self.ref_len = len(reference)
        self.reference = map_bases(reference)
        self.nick_sites = nick_sites
        # Scoring information
        self.match, self.degen_match, self.mismatch = match, degen_match, mismatch
        self.gap_open, self.gap_ext, self.gap_rep = gap_open, gap_ext, gap_rep
        self.gap_incentive, self.het_score = gap_incentive, het_score
        # Set gap incentives for each index
        self.gap_incentives = <float *>PyMem_Malloc((len(reference) + 1) * sizeof(float))
        for i in range(len(reference) + 1):
            self.gap_incentives[i] = 0
        for idx in nick_sites:
            self.gap_incentives[idx] = gap_incentive
            self.gap_incentives[idx + 1] = gap_incentive
        # Heterogeneous regions
        self.num_hets = reference.count('*')
        self.het_idxs = <size_t *>PyMem_Malloc(self.num_hets * sizeof(size_t))
        i = 0
        for ref_i, c in enumerate(reference, 1):
            if c == '*':
                self.het_idxs[i] = ref_i
                i += 1
        # Repeat information
        self.num_repeats = num_repeats
        self.repeat_idxs = <size_t *>PyMem_Malloc(num_repeats * sizeof(size_t))
        for i, c in enumerate(reference, 1):
            if '0' <= c <= '9':
                self.repeat_idxs[ord(c) - 48] = i
        self.repeat_lens = <size_t *>PyMem_Malloc(num_repeats * sizeof(size_t))
        self.repeats = <Base **> PyMem_Malloc(num_repeats * sizeof(Base *))
        for i in range(num_repeats):
            self.repeats[i] = map_bases(repeats[i])
            self.repeat_lens[i] = len(repeats[i])
        self.repeat_strs = repeats
        # Initialize buffers
        self.max_query_len = init_query_len
        self._score_data = NULL
        self.m_arr = NULL
        self.i_arr = NULL
        self.d_arr = NULL
        self.r_arrs = <float **>PyMem_Malloc(num_repeats * sizeof(float *))
        self.tb_arr = NULL
        self.make_buffers(init_query_len)

    def __dealloc__(self):
        PyMem_Free(self.reference)
        PyMem_Free(self.gap_incentives)
        PyMem_Free(self.het_idxs)
        PyMem_Free(self.repeat_idxs)
        PyMem_Free(self.repeat_lens)
        for k in range(self.num_repeats):
            PyMem_Free(self.repeats[k])
        PyMem_Free(self.repeats)
        PyMem_Free(self.r_arrs)
        PyMem_Free(self.tb_arr)

    cdef void make_buffers(self, size_t max_query_len):
        # Calculate size needed for scoring matrices
        cdef size_t ref_dim = self.ref_len + 1
        cdef size_t qry_dim = max_query_len + 1
        cdef size_t mid_entries = 3 * ref_dim * qry_dim
        cdef size_t rep_entries = 0
        for k in range(self.num_repeats):
            rep_entries += self.repeat_lens[k] * qry_dim
        cdef size_t total_entries = mid_entries + rep_entries
        # Allocate and assign pointers to scoring matrices
        self._score_data = <float *>PyMem_Realloc(self._score_data, total_entries * sizeof(float))
        self.max_query_len = max_query_len
        self.m_arr = &self._score_data[0 * ref_dim * qry_dim]
        self.i_arr = &self._score_data[1 * ref_dim * qry_dim]
        self.d_arr = &self._score_data[2 * ref_dim * qry_dim]
        cdef size_t rep_idx = 3 * ref_dim * qry_dim
        for k in range(self.num_repeats):
            self.r_arrs[k] = &self._score_data[rep_idx]
            rep_idx += self.repeat_lens[k] * qry_dim
        # Allocate traceback arrays
        self.tb_arr =  <Direction *>PyMem_Realloc(self.tb_arr, mid_entries * sizeof(Direction))

    cdef inline size_t idx(self, size_t r_idx, size_t q_idx):
        return r_idx * (self.max_query_len + 1) + q_idx

    cdef float m_value(self, Base *query, size_t ref_i, size_t qry_i):
        # If in a repeat, then immediately return NINF
        for k in range(self.num_repeats):
            if ref_i == self.repeat_idxs[k]:
                return NINF
        # If in a heterogeneous region, then add heterogeneous penalty
        cdef float m_prev, h_prev
        for k in range(self.num_hets):
            if ref_i == self.het_idxs[k]:
                m_prev = self.m_arr[self.idx(ref_i - 1, qry_i - 1)]
                h_prev = self.m_arr[self.idx(ref_i,     qry_i - 1)]
                if m_prev >= h_prev:
                    self.tb_arr[self.idx(ref_i, qry_i)] = DIAG
                else:
                    self.tb_arr[self.idx(ref_i, qry_i)] = LEFT
                return max(m_prev, h_prev) + self.het_score
        # If right after a repeat, then consider all possible repeat values
        cdef size_t after_rep = 0
        cdef size_t rep_len
        cdef float best_score
        for k in range(self.num_repeats):
            if ref_i == self.repeat_idxs[k] + 1:
                after_rep = 1
                break
        if after_rep:
            rep_len = self.repeat_lens[k]
            best_score = self.r_arrs[k][self.idx(rep_len - 1, qry_i - 1)]
            for r_i in range(rep_len - 1):
                best_score = max(best_score,
                                 self.r_arrs[k][self.idx(r_i, qry_i - 1)] + self.gap_rep)
            return best_score
        # Otherwise, typical Gotoh-style alignment
        cdef size_t prev_idx = self.idx(ref_i - 1, qry_i - 1)
        cdef float max_prev = max(self.m_arr[prev_idx],
                                  self.i_arr[prev_idx],
                                  self.d_arr[prev_idx])
        cdef float score
        if self.reference[ref_i - 1] == query[qry_i - 1]:
            score = self.match
        elif self.reference[ref_i - 1] & query[qry_i - 1]:
            score = self.degen_match
        else:
            score = self.mismatch
        return max_prev + score

    cdef float i_value(self, size_t ref_i, size_t qry_i):
        # If repeat or heterogeneous region, return NINF
        for k in range(self.num_repeats):
            if ref_i == self.repeat_idxs[k]:
                return NINF
        for k in range(self.num_hets):
            if ref_i == self.het_idxs[k]:
                return NINF
        cdef size_t prev_idx = self.idx(ref_i, qry_i - 1)
        return max(self.m_arr[prev_idx] + self.gap_open + self.gap_incentives[ref_i],
                   self.i_arr[prev_idx] + self.gap_ext)

    cdef float d_value(self, size_t ref_i, size_t qry_i):
        # If repeat or heterogeneous region, return NINF
        for k in range(self.num_repeats):
            if ref_i == self.repeat_idxs[k]:
                return NINF
        for k in range(self.num_hets):
            if ref_i == self.het_idxs[k]:
                return NINF
        cdef size_t prev_idx = self.idx(ref_i - 1, qry_i)
        return max(self.m_arr[prev_idx] + self.gap_open + self.gap_incentives[ref_i],
                   self.d_arr[prev_idx] + self.gap_ext)

    cdef float r_value(self, Base *query, size_t k, size_t rep_i, size_t qry_i):
        # Match/mismatch score
        cdef float score
        if query[qry_i - 1] == self.repeats[k][rep_i]:
            score = self.match
        elif query[qry_i - 1] & self.repeats[k][rep_i]:
            score = self.degen_match
        else:
            score = self.mismatch
        # Compute scores
        cdef size_t prev_rep_i = rep_i - 1 if rep_i > 0 else self.repeat_lens[k] - 1
        cdef float r_m_score = self.r_arrs[k][self.idx(prev_rep_i, qry_i - 1)] + 0.1 * score
        cdef float r_i_score = self.r_arrs[k][self.idx(rep_i,      qry_i - 1)] + 0.1 * self.gap_rep
        cdef float r_d_score = self.r_arrs[k][self.idx(prev_rep_i, qry_i    )] + 0.1 * self.gap_rep
        cdef float best_score = max(r_m_score, r_i_score, r_d_score)
        # Best score/direction out of repeats
        cdef size_t prev_idx
        cdef float m_score
        cdef Direction m_dir
        if rep_i == 0:
            prev_idx = self.idx(self.repeat_idxs[k] - 1, qry_i - 1)
            m_score = self.m_arr[prev_idx] + 0.1 * score
            m_dir = LEFT
            if m_score >= best_score:
                best_score = m_score
                m_dir = DIAG
            self.tb_arr[self.idx(self.repeat_idxs[k], qry_i)] = m_dir
        # Best direction within repeats
        cdef Direction r_dir = DIAG
        if r_i_score == best_score:
            r_dir = LEFT
        elif r_d_score == best_score:
            r_dir = UP
        return best_score

    cdef void fill_buffers(self, Base *query, size_t q_len):
        # Initialize edges
        cdef size_t idx = self.idx(0, 0)
        self.m_arr[idx] = 0
        self.i_arr[idx] = NINF
        self.d_arr[idx] = NINF
        self.tb_arr[idx] = <Direction>0
        for ref_i in range(1, self.ref_len + 1):
            idx = self.idx(ref_i, 0)
            self.m_arr[idx] = NINF
            self.i_arr[idx] = NINF
            self.d_arr[idx] = ref_i * self.gap_ext
            self.tb_arr[idx] = UP
        for qry_i in range(1, q_len + 1):
            idx = self.idx(0, qry_i)
            self.m_arr[idx] = NINF
            self.i_arr[idx] = qry_i * self.gap_ext
            self.d_arr[idx] = NINF
            self.tb_arr[idx] = LEFT
        for k in range(self.num_repeats):
            for rep_i in range(self.repeat_lens[k]):
                for qry_i in range(q_len + 1):
                    idx = self.idx(rep_i, qry_i)
                    self.r_arrs[k][idx] = NINF
        # Fill in matrix
        cdef float max_score
        cdef Direction max_direction
        cdef size_t prev_idx
        cdef size_t done = 0
        for ref_i in range(1, self.ref_len + 1):
            for qry_i in range(1, q_len + 1):
                idx = self.idx(ref_i, qry_i)
                self.m_arr[idx] = self.m_value(query, ref_i, qry_i)
                self.i_arr[idx] = self.i_value(ref_i, qry_i)
                self.d_arr[idx] = self.d_value(ref_i, qry_i)
                done = 0
                # Repeat regions handle direction arrays while computing r_value
                for k in range(self.num_repeats):
                    if ref_i == self.repeat_idxs[k]:
                        for rep_i in range(self.repeat_lens[k]):
                            self.r_arrs[k][self.idx(rep_i, qry_i)] = \
                                self.r_value(query, k, rep_i, qry_i)
                        done = 1
                        break
                if done:
                    continue
                # Heterogeneous regions handle direction arrays while computing m_value
                for k in range(self.num_hets):
                    if ref_i == self.het_idxs[k]:
                        done = 1
                        break
                if done:
                    continue
                # For constant regions, determine direction based on array values
                max_score = max(self.m_arr[idx], self.i_arr[idx], self.d_arr[idx])
                max_direction = <Direction>0
                if is_close(self.m_arr[idx], max_score):
                    max_direction = DIAG
                elif is_close(self.i_arr[idx], max_score):
                    max_direction = LEFT
                elif is_close(self.d_arr[idx], max_score):
                    max_direction = UP
                self.tb_arr[idx] = max_direction

    @staticmethod
    def _new_size(query_len):
        size = 1
        while size <= query_len:
            size *= 2
        return size - 1

    def alignment_str(self, query: str):
        max_score = NINF
        max_idx = -1
        q_len = len(query)
        idx = self.idx(self.ref_len, q_len)
        for ref_i in range(self.ref_len + 1):
            # Check if ref_i is a repeat
            for k in range(self.num_repeats):
                if ref_i == self.repeat_idxs[k]:
                    max_idx_score = NINF
                    for rep_i in range(self.repeat_lens[k]):
                        max_idx_score = max(max_idx_score, self.r_arrs[k][self.idx(rep_i, q_len)])
                    break
            else:
                idx = self.idx(ref_i, q_len)
                max_idx_score = max(self.m_arr[idx], self.i_arr[idx], self.d_arr[idx])
            # If new max, set max_score and max_idx
            if max_idx_score > max_score:
                max_score = max_idx_score
                max_idx = ref_i
        # Assemble alignment string
        rev_ref = []
        rev_qry = []
        for ref_i in range(self.ref_len, max_idx, -1):
            rev_ref.append(self.ref_str[ref_i - 1])
            rev_qry.append('-')
        ref_i = max_idx
        qry_i = q_len
        direction = self.tb_arr[self.idx(max_idx, q_len)]
        while direction:
            if direction & DIAG:
                rev_ref.append(self.ref_str[ref_i - 1])
                rev_qry.append(query[qry_i - 1])
                ref_i -= 1
                qry_i -= 1
            elif direction & LEFT:
                if '0' <= self.ref_str[ref_i - 1] <= '9' or self.ref_str[ref_i - 1] == '*':
                    rev_ref.append(self.ref_str[ref_i - 1])
                else:
                    rev_ref.append('-')
                rev_qry.append(query[qry_i - 1])
                qry_i -= 1
            else:
                rev_ref.append(self.ref_str[ref_i - 1])
                rev_qry.append('-')
                ref_i -= 1
            direction = self.tb_arr[self.idx(ref_i, qry_i)]
        return ''.join(rev_ref)[::-1] + '\n' + ''.join(rev_qry)[::-1]

    @staticmethod
    def score_alignment(str alignment) -> float:
        ref_align, qry_align = alignment.split('\n')
        ref_trunc = ref_align[:len(qry_align.rstrip('-'))]
        qry_len = (qry_align.count('A') + qry_align.count('C') +
                   qry_align.count('G') + qry_align.count('T'))
        ref_len = len(ref_trunc) - ref_trunc.count('-')
        raw_score = 0.0
        for r, q in zip(ref_align, qry_align):
            if r in 'ACGT' and r == q:
                raw_score += 1
            elif r in '0123456789' and q != '-':
                raw_score += 1
        qry_frac = raw_score / qry_len
        ref_frac = raw_score / ref_len
        # Harmonic mean of qry_frac and ref_frac
        return 200 / (1 / qry_frac + 1 / ref_frac)

    cdef align(self, str query):
        if <size_t>len(query) > self.max_query_len:
            self.make_buffers(PowtnrkaAligner._new_size(len(query)))
        cdef Base *q_bases = map_bases(query)
        self.fill_buffers(q_bases, len(query))
        PyMem_Free(q_bases)

    def __call__(self, str query) -> Tuple[str, float]:
        self.align(query)
        alignment = self.alignment_str(query)
        return alignment, PowtnrkaAligner.score_alignment(alignment)

    def print_directions(self, str query):
        dir_map = ['0', 'U', 'L', '3', 'D', '5', '6', '7']
        print('  ' + query)
        print_str = [' ']
        for q_i in range(len(query) + 1):
            print_str.append(dir_map[int(self.tb_arr[self.idx(0, q_i)])])
        print_str = ''.join(print_str)
        print(print_str)
        for r_i in range(1, self.ref_len + 1):
            print_str = [self.ref_str[r_i - 1]]
            for q_i in range(len(query) + 1):
                print_str.append(dir_map[int(self.tb_arr[self.idx(r_i, q_i)])])
            print_str = ''.join(print_str)
            print(print_str)


cdef class RepeatAligner:
    # Scoring information
    cdef float match, degen_match, mismatch, gap_penalty
    # Reference information
    cdef readonly str ref_str
    cdef readonly str aux_str
    cdef size_t ref_len
    cdef Base *reference
    # Scoring/traceback arrays
    cdef size_t max_query_len
    cdef float *score_data
    cdef size_t *tb_arr
    # Subalign matrix
    cdef float *sub_score_data
    # Score cache
    cdef dict score_cache

    def __init__(self, str reference, float match = 5.0, float degen_match = 1.0,
                 float mismatch = -2.0, float gap_penalty = -5.0, size_t init_query_len = 64,
                 str aux_str: Optional[str] = None):
        # Reference information
        self.ref_str = reference
        self.aux_str = aux_str
        self.ref_len = len(reference)
        self.reference = map_bases(reference)
        # Scoring information
        self.match, self.degen_match, self.mismatch = match, degen_match, mismatch
        self.gap_penalty = gap_penalty
        # Initialize buffers
        self.max_query_len = 0
        self.score_data = NULL
        self.tb_arr = NULL
        self.make_buffers(init_query_len)
        self.sub_score_data = <float *>PyMem_Malloc((self.ref_len + 1) *
                                                    (2 * self.ref_len) * sizeof(float))
        self.score_cache = {}

    def __dealloc__(self):
        PyMem_Free(self.reference)
        PyMem_Free(self.score_data)
        PyMem_Free(self.tb_arr)
        PyMem_Free(self.sub_score_data)

    @staticmethod
    def _new_size(query_len):
        size = 1
        while size < query_len:
            size *= 2
        return size

    cdef void make_buffers(self, size_t max_query_len):
        self.score_data = <float *>PyMem_Realloc(self.score_data, max_query_len * sizeof(float))
        self.tb_arr = <size_t *>PyMem_Realloc(self.tb_arr, max_query_len * sizeof(size_t))
        self.max_query_len = max_query_len

    cdef inline size_t idx(self, size_t r_idx, size_t q_idx):
        return r_idx * (2 * self.ref_len) + q_idx

    cdef float sub_score(self, str substr):
        if substr in self.score_cache:
            return self.score_cache[substr]
        cdef Base *query = map_bases(substr)
        cdef size_t qry_len = len(substr)
        # Initialize subscore matrix
        self.sub_score_data[self.idx(0, 0)] = 0
        for r_idx in range(1, self.ref_len + 1):
            self.sub_score_data[self.idx(r_idx, 0)] = r_idx * self.gap_penalty
        for q_idx in range(1, qry_len + 1):
            self.sub_score_data[self.idx(0, q_idx)] = q_idx * self.gap_penalty
        cdef float m_score, i_score, d_score
        for r_idx in range(1, self.ref_len + 1):
            for q_idx in range(1, qry_len + 1):
                if self.reference[r_idx - 1] == query[q_idx - 1]:
                    m_score = self.match
                elif self.reference[r_idx - 1] & query[q_idx - 1]:
                    m_score = self.degen_match
                else:
                    m_score = self.mismatch
                m_score = self.sub_score_data[self.idx(r_idx - 1, q_idx - 1)] + m_score
                i_score = self.sub_score_data[self.idx(r_idx, q_idx - 1)] + self.gap_penalty
                d_score = self.sub_score_data[self.idx(r_idx - 1, q_idx)] + self.gap_penalty
                self.sub_score_data[self.idx(r_idx, q_idx)] = max(m_score, i_score, d_score)
        cdef float max_score = NINF
        for r_idx in range(self.ref_len + 1):
            max_score = max(max_score, self.sub_score_data[self.idx(r_idx, qry_len)])
        if qry_len != self.ref_len:
            max_score = max_score + self.gap_penalty
        self.score_cache[substr] = max_score
        return max_score

    cdef void fill_buffers(self, str query):
        cdef float score
        cdef float best_score
        cdef size_t best_size
        cdef size_t start
        for i in range(len(query)):
            best_score = NINF
            best_size = 0
            start = max(0, i + 2 - 2 * self.ref_len)
            for j in range(start, i + 1):
                score = self.sub_score(query[j:i + 1])
                if j > 0:
                    score = score + self.score_data[j - 1]
                if score > best_score:
                    best_score = score
                    best_size = i - j + 1
            self.score_data[i] = best_score
            self.tb_arr[i] = best_size

    def traceback(self, str query):
        repeats = []
        while query:
            last_size = self.tb_arr[len(query) - 1]
            query, repeat = query[:-last_size], query[-last_size:]
            repeats.append(repeat)
        repeats.reverse()
        return repeats

    def __call__(self, str query):
        if <size_t>len(query) > self.max_query_len:
            self.make_buffers(RepeatAligner._new_size(len(query)))
        self.fill_buffers(query)
        return self.traceback(query)
