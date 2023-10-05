# -*- coding: utf-8 -*-
# SPDX-FileCopyrightText: : 2020-2023 The PyPSA-Eur Authors
#
# SPDX-License-Identifier: MIT
"""

"""

from __future__ import absolute_import, print_function

import logging
import os
import sys
import time

logger = logging.getLogger(__name__)

# TODO: provide alternative when multiprocessing is not available
try:
    from multiprocessing import Pipe, Process
except ImportError:
    from multiprocessing.dummy import Process, Pipe

from memory_profiler import _get_memory, choose_backend


# The memory logging facilities have been adapted from memory_profiler
class MemTimer(Process):
    """
    Write memory consumption over a time interval to file until signaled to
    stop on the pipe.
    """

    def __init__(
        self, monitor_pid, interval, pipe, filename, max_usage, backend, *args, **kw
    ):
        self.monitor_pid = monitor_pid
        self.interval = interval
        self.pipe = pipe
        self.filename = filename
        self.max_usage = max_usage
        self.backend = backend

        self.timestamps = kw.pop("timestamps", True)
        self.include_children = kw.pop("include_children", True)

        super(MemTimer, self).__init__(*args, **kw)

    def run(self):
        # get baseline memory usage
        cur_mem = _get_memory(
            self.monitor_pid,
            self.backend,
            timestamps=self.timestamps,
            include_children=self.include_children,
        )

        n_measurements = 1
        mem_usage = cur_mem if self.max_usage else [cur_mem]

        if self.filename is not None:
            stream = open(self.filename, "w")
            stream.write("MEM {0:.6f} {1:.4f}\n".format(*cur_mem))
            stream.flush()
        else:
            stream = None

        self.pipe.send(0)  # we're ready
        stop = False
        while True:
            cur_mem = _get_memory(
                self.monitor_pid,
                self.backend,
                timestamps=self.timestamps,
                include_children=self.include_children,
            )

            if stream is not None:
                stream.write("MEM {0:.6f} {1:.4f}\n".format(*cur_mem))
                stream.flush()

            n_measurements += 1
            if not self.max_usage:
                mem_usage.append(cur_mem)
            else:
                mem_usage = max(cur_mem, mem_usage)

            if stop:
                break
            stop = self.pipe.poll(self.interval)
            # do one more iteration

        if stream is not None:
            stream.close()

        self.pipe.send(mem_usage)
        self.pipe.send(n_measurements)


class memory_logger(object):
    """
    Context manager for taking and reporting memory measurements at fixed
    intervals from a separate process, for the duration of a context.

    Parameters
    ----------
    filename : None|str
        Name of the text file to log memory measurements, if None no log is
        created (defaults to None)
    interval : float
        Interval between measurements (defaults to 1.)
    max_usage : bool
        If True, only store and report the maximum value (defaults to True)
    timestamps : bool
        Whether to record tuples of memory usage and timestamps; if logging to
        a file timestamps are always kept (defaults to True)
    include_children : bool
        Whether the memory of subprocesses is to be included (default: True)

    Arguments
    ---------
    n_measurements : int
        Number of measurements that have been taken
    mem_usage : (float, float)|[(float, float)]
        All memory measurements and timestamps (if timestamps was True) or only
        the maximum memory usage and its timestamp

    Note
    ----
    The arguments are only set after all the measurements, i.e. outside of the
    with statement.

    Example
    -------
    with memory_logger(filename="memory.log", max_usage=True) as mem:
        # Do a lot of long running memory intensive stuff
        hard_memory_bound_stuff()

    max_mem, timestamp = mem.mem_usage
    """

    def __init__(
        self,
        filename=None,
        interval=1.0,
        max_usage=True,
        timestamps=True,
        include_children=True,
    ):
        if filename is not None:
            timestamps = True

        self.filename = filename
        self.interval = interval
        self.max_usage = max_usage
        self.timestamps = timestamps
        self.include_children = include_children

    def __enter__(self):
        backend = choose_backend()

        self.child_conn, self.parent_conn = Pipe()  # this will store MemTimer's results
        self.p = MemTimer(
            os.getpid(),
            self.interval,
            self.child_conn,
            self.filename,
            backend=backend,
            timestamps=self.timestamps,
            max_usage=self.max_usage,
            include_children=self.include_children,
        )
        self.p.start()
        self.parent_conn.recv()  # wait until memory logging in subprocess is ready

        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        if exc_type is None:
            self.parent_conn.send(0)  # finish timing

            self.mem_usage = self.parent_conn.recv()
            self.n_measurements = self.parent_conn.recv()
        else:
            self.p.terminate()

        return False


class timer(object):
    level = 0
    opened = False

    def __init__(self, name="", verbose=True):
        self.name = name
        self.verbose = verbose

    def __enter__(self):
        if self.verbose:
            if self.opened:
                sys.stdout.write("\n")

            if len(self.name) > 0:
                sys.stdout.write((".. " * self.level) + self.name + ": ")
            sys.stdout.flush()

            self.__class__.opened = True

        self.__class__.level += 1

        self.start = time.time()
        return self

    def print_usec(self, usec):
        if usec < 1000:
            print("%.1f usec" % usec)
        else:
            msec = usec / 1000
            if msec < 1000:
                print("%.1f msec" % msec)
            else:
                sec = msec / 1000
                print("%.1f sec" % sec)

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.opened and self.verbose:
            sys.stdout.write(".. " * self.level)

        if exc_type is None:
            stop = time.time()
            self.usec = usec = (stop - self.start) * 1e6
            if self.verbose:
                self.print_usec(usec)
        elif self.verbose:
            print("failed")
        sys.stdout.flush()

        self.__class__.level -= 1
        if self.verbose:
            self.__class__.opened = False
        return False


class optional(object):
    def __init__(self, variable, contextman):
        self.variable = variable
        self.contextman = contextman

    def __enter__(self):
        if self.variable:
            return self.contextman.__enter__()

    def __exit__(self, exc_type, exc_val, exc_tb):
        if self.variable:
            return self.contextman.__exit__(exc_type, exc_val, exc_tb)
        return False
