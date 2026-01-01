#!/usr/bin/env python

import time
import datetime
import logging

class TimeData:
    """Timing and profiling utilities for tracking execution time."""
    def __init__(self, log=None):
        """Initialize timing data structures.

        Parameters
        ----------
        log : logging.Logger, optional
            Logger instance.
        """
        self.log = log or logging.getLogger(__name__)
        if not self.log.handlers:
            logging.basicConfig(level=logging.INFO, format="%(message)s")
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
        """Display timing for a specific clock.

        Parameters
        ----------
        clock_name : str
            Name of the clock.
        """
        total = self.clock_total[clock_name]
        count = self.clock_count[clock_name]
        self.log.info(f"{clock_name:25s}: {total/60:12.2f}m  ({total/count:12.2f}s for each call)")

    def show_time(self):
        """Display total elapsed time."""
        tot_time = time.time() - self.start_time
        self.log.info(f"total time               : {tot_time/60:12.2f}m")

    def get_time(self):
        return time.time() - self.start_time

    def show_all(self):
        """Display all timing statistics."""
        tot_time = time.time() - self.start_time
        self.log.info("")
        self.log.info(f"total time               : {tot_time/60:12.2f}m")
        for key in self.clock_count.keys():
            key_tot = self.clock_total.get(key, 0)
            key_count = self.clock_count.get(key, 1)
            self.log.info(f"{key:25s}: {key_tot/60:12.2f}m  ({key_tot/key_count:12.2f}s for each call)")

