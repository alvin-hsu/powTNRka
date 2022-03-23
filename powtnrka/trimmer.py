# Copyright (C) 2024   Alvin Hsu
from typing import Tuple

import numpy as np


class PowtnrkaTrimmer(object):
    def __init__(self, window_size: int, min_q: float):
        self.window_size = window_size
        self.min_q = min_q
        if window_size % 2 == 0:
            offset_size = window_size // 2
            self.l_offset = offset_size
            self.r_offset = offset_size
        else:
            offset_size = (window_size - 1) // 2
            self.l_offset = offset_size
            self.r_offset = offset_size + 1

    def __call__(self, read: str, qual: np.ndarray) -> Tuple[str, np.ndarray]:
        i = 0
        for i in range(self.l_offset, len(read) - self.r_offset):
            mean_q = qual[i - self.l_offset : i + self.r_offset].mean()
            if mean_q < self.min_q:
                break
        return read[:i], qual[:i]


class NullTrimmer(PowtnrkaTrimmer):
    def __init__(self):
        super().__init__(1, 0)

    def __call__(self, read: str, qual: np.ndarray) -> Tuple[str, np.ndarray]:
        return read, qual
