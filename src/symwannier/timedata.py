#!/usr/bin/env python

import time
import datetime

class TimeData:
    def __init__(self):
        self.start_time = time.time()
        self.start_datetime = datetime.datetime.now()
        self.clock_count = {}
        self.clock_start = {}
        self.clock_total = {}

    def start_clock(self, clock_name):
        self.clock_count[clock_name] = self.clock_count.get(clock_name, 0) + 1
        self.clock_start[clock_name] = time.time()

    def stop_clock(self, clock_name):
        self.clock_total[clock_name] = time.time() - self.clock_start.get(clock_name, self.start_time) + self.clock_total.get(clock_name, 0)

    def show_clock(self, clock_name):
        total = self.clock_total[clock_name]
        count = self.clock_count[clock_name]
        print(" {:25s}: {:12.2f}m  ({:12.2f}s for each call) ".format(clock_name, total / 60, total / count))

    def show_time(self):
        tot_time = time.time() - self.start_time
        print(" total time               : {:12.2f}m".format(tot_time / 60))

    def get_time(self):
        return time.time() - self.start_time

    def show_all(self):
        tot_time = time.time() - self.start_time
        print("")
        print(" total time               : {:12.2f}m".format(tot_time / 60))
        for key in self.clock_count.keys():
            key_tot = self.clock_total.get(key, 0)
            key_count = self.clock_count.get(key, 1)
            print(" {:25s}: {:12.2f}m  ({:12.2f}s for each call) ".format(key, key_tot / 60, key_tot / key_count))

