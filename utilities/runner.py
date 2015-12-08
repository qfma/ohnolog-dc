#!/usr/bin/env python2.7
''' Run different programs

This file provided various python functions that use the os and sys
modules in order to list files, folders etc.

'''
import sys
import os
import time
from subprocess import Popen, list2cmdline
import subprocess
import multiprocessing
import functools
import logging
import pprint

log = logging.getLogger('pipeline')


def cpu_count():
    '''Returns the number of CPUs in the system'''
    num = 1
    if sys.platform == 'win32':
        try:
            num = int(os.environ['NUMBER_OF_PROCESSORS'])
        except (ValueError, KeyError):
            pass
    elif sys.platform == 'darwin':
        try:
            num = int(os.popen('sysctl -n hw.ncpu').read())
        except ValueError:
            pass
    else:
        try:
            num = os.sysconf('SC_NPROCESSORS_ONLN')
        except (ValueError, OSError, AttributeError):
            pass
    return num


def exec_commands(cmds, max_task=cpu_count()):
    '''
    Exec commands in parallel in multiple process. By default it uses
    cpu_count() to determine the number of processors and uses all of them.
    The user can specify a lower number manually using the max_task parameter.
    '''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)
    # If the number of processors set is higher than available CPUs
    # use all.
    if max_task > cpu_count():
        max_task = cpu_count()
    processes = []
    # While loop submits commands to precessors.
    while True:
        while cmds and len(processes) < max_task:
            task = cmds.pop()
            print list2cmdline(task)
            processes.append(Popen(task))

        for p in processes:
            if done(p):
                if success(p):
                    processes.remove(p)
                else:
                    fail()

        if not processes and not cmds:
            break
        else:
            time.sleep(0.05)


def exec_in_row(cmds):
    '''
    Exec commands one after the other until finished. This is helpful
    if a program is already parallelized, but we have to submit many jobs
    after each other.
    '''
    if not cmds:
        return  # empty list

    def done(p):
        return p.poll() is not None

    def success(p):
        return p.returncode == 0

    def fail():
        sys.exit(1)

    for task in cmds:
        print task
        p = Popen(task)
        p.wait()

    if done(p):
            if success(p):
                print "done!"
            else:
                fail()


def sub_call(command, stdout=None, shell=False):
    '''
    Subprocess manager: will log commands and raise an error
    if return code does not equal 0.
    '''
    pretty_command = pprint.pformat(command, indent=4)
    log.debug('running command:\n{}'.format(pretty_command))
    if stdout:
        log.debug('writing to: ' + stdout.name)
    subprocess.check_call(command, stdout=stdout, shell=shell)


def multiprocess(iterator):
    '''
    Returns the multiprocess decorator. This function's only purpose
    is to pass its iterator argument along.
    '''

    def multiprocess_decorator(func):
        '''
        Returns the multiprocess wrapper
        '''

        @functools.wraps(func)
        def multiprocess_wrapper():
            '''
            Multiprocesses a given function with arguments specified
            by the iterator.  This will create a separate process for each
            item in the iterator.
            '''
            # log.info('#' * 45 + ' -> start')
            # log.info('Calling multi_processor with target: ' + func.__name__)
            jobs = []
            for it in iterator:
                if isinstance(it, tuple):
                    j = multiprocessing.Process(target=func,
                                                name=(func.__name__ + ' ' + str(it)),
                                                args=it)
                else:
                    j = multiprocessing.Process(target=func,
                                                name=(func.__name__ + ' ' + it),
                                                args=(it,))
                jobs.append(j)
                j.start()

            for j in jobs:
                j.join()
            for j in jobs:
                if j.exitcode != 0:
                    raise AssertionError('multi process returned with non-0 '
                                         'exit code')
            # log.info('#' * 45 + ' end <-')
        return multiprocess_wrapper

    return multiprocess_decorator
