#!/usr/bin/env python

import numpy as np
import os

class Win:
    def __init__(self, prefix):
        file_win = prefix + ".win"
        assert os.path.exists(file_win)

        self._read_win_file(file_win)

    def _read_win_file(self, file_win):
        with open (file_win) as fp:
            self.lines = fp.readlines()

        self.num_bands = self._get_param_keyword("num_bands", 0, dtype=int)
        self.num_wann = self._get_param_keyword("num_wann", 0, dtype=int)

        if self.num_bands > self.num_wann:
            self.dis_num_iter = self._get_param_keyword("dis_num_iter", 200, dtype=int)
        else:
            self.dis_num_iter = 0
        self.num_iter = self._get_param_keyword("num_iter", 200, dtype=int)
        self.dis_froz_min = self._get_param_keyword("dis_froz_min", -100000, dtype=float)
        self.dis_froz_max = self._get_param_keyword("dis_froz_max", +100000, dtype=float)
        self.dis_win_min = self._get_param_keyword("dis_win_min", -100000, dtype=float)
        self.dis_win_max = self._get_param_keyword("dis_win_max", +100000, dtype=float)
        self.dis_mix_ratio = self._get_param_keyword("dis_mix_ratio", 0.5, dtype=float)

        p = self._get_param_keyword("mp_grid", dtype = str)
        if p is not None:
            self.mp_grid = np.array([ int(x) for x in p.split() ])
        else:
            raise Exception('mp_grid is not defined')

    def _get_param_keyword(self, keyword, default_value = None, dtype = int):
        data = None
        for line in self.lines:
            if line.startswith(keyword):
                assert data == None, keyword + " is defined more than once"
                if len(line.split("=")) > 1:
                    data = line.split("=")[1]
                elif len(line.split(":")) > 1:
                    data = line.split(":")[1]
        if data == None:
            data = default_value
        if data == None:
            return None
        return dtype(data)
