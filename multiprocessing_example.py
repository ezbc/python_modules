#!/usr/bin/python


import os
from Queue import Queue
from multiprocessing import Process

def returning_wrapper(func, *args, **kwargs):
    queue = kwargs.get("multiprocess_returnable")
    del kwargs["multiprocess_returnable"]
    queue.put(func(*args, **kwargs))

class Multiprocess(object):
    """Cute decorator to run a function in multiple processes."""
    def __init__(self, func):
        self.func = func
        self.processes = []

    def __call__(self, *args, **kwargs):
        num_processes = kwargs.get("multiprocess_num_processes", 12) # default to two processes.
        return_obj = kwargs.get("multiprocess_returnable", Queue()) # default to stdlib Queue
        kwargs["multiprocess_returnable"] = return_obj
        for i in xrange(num_processes):
            pro = Process(target=returning_wrapper, args=tuple([self.func] + list(args)), kwargs=kwargs)
            self.processes.append(pro)
            pro.start()
        return return_obj

@Multiprocess
def run():
    import numpy as np
    n = np.random.randint(0, 1, size=(100, 1000))
    for i in xrange(n.shape[0]):
        for j in xrange(n.shape[1]):
        	a = n * i * j

    return a

def for_loop(n):
    for i in xrange(n.shape[0]):
        for j in xrange(n.shape[1]):
        	a = n * i * j

    return a

def main():

    import numpy as np

    # http://stackoverflow.com/questions/10797998/is-it-possible-to-multiprocesspython-a-function-that-returns-something
    '''
    items = (10, 10, 10, 10)

    x = split_processing(items)

    print x

    info('main line')
    p = Process(target=f, args=('bob',))
    p.start()
    p.join()

    data = info()
    print 'data', data
    data = run()
    #print data.get(False)

    print 'stage 1 done'

    @Multiprocess
    def info():
        print 'module name:', __name__
        print 'parent process:', os.getppid()
        print 'process id:', os.getpid()
        return 4 * 22

    data = info()
    print data.get(False)

    print 'stage 2 done'

    data = Multiprocess(for_loop)
    #data = info(for_loop)

    array = np.random.randint(0, 1, size=(100,100))

    output = data(*(array,))

    print 'before'
    #print help(output)
    print output.empty()
    print output.get(False)
    #print(help(data))

    '''

    from multiprocessing import Pool

    array = np.random.randint(0, 1, size=(100,100))

    pool = Pool()
    pool.map(for_loop, (array,))
    pool.get(timeout=10)


if __name__ == '__main__':
	main()

