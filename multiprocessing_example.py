#!/usr/bin/python


import os
from Queue import Queue
from multiprocessing import Process
import signal

class KeyboardInterruptError(Exception): pass

def manip_data(data):
    try:
        for datum in data:
            test = datum * data
    except KeyboardInterrupt:
        raise KeyboardInterruptError

    return 'no'

def manip_data_args(data, arg1=0):
    try:
        for datum in data:
            test = datum * data
    except KeyboardInterrupt:
        raise KeyboardInterruptError

    return arg1

def init_worker():
    signal.signal(signal.SIGINT, signal.SIG_IGN)

def main():

    import numpy as np

    # http://stackoverflow.com/questions/10797998/is-it-possible-to-multiprocesspython-a-function-that-returns-something

    import multiprocessing
    import time

    data = (
        ['a', '2'], ['b', '4'], ['c', '6'], ['d', '8'],
        ['e', '1'], ['f', '3'], ['g', '5'], ['h', '7']
    )

    data = np.random.random((1000, 10000))

    p = multiprocessing.Pool()
    args = np.arange(0, 20, 1)
    results = np.zeros(args.shape)
    for i, arg in enumerate(args):
        results[i] = p.apply_async(manip_data_args, args=(data, arg))
    results.get()
    p.close()
    print 'pool map complete'

    '''
    try:
        p = multiprocessing.Pool()
        test = p.map(manip_data, data)
        p.close()
        print 'pool map complete'
    except KeyboardInterrupt:
        print 'got ^C while pool mapping, terminating the pool'
        p.terminate()
        p.wait()
        print 'pool is terminated'
    finally:
        print 'joining pool processes'
        p.join()
        print 'join complete'
    print 'the end'
    '''

if __name__ == '__main__':
	main()

