import pickle

def check_file(filename, clobber=False):

    ''' Determines if a file exists.
    '''

    import os

    exists = False

    if os.path.isfile(filename) or os.path.isdir(filename):
        exists = True
        if clobber:
            os.system('rm -rf {:s}'.format(filename))
            exists = False

    return exists

def save_pickle(filename, data):
    with open(filename, 'wb') as f:
        pickle.dump(data, f)
    f.close()

def load_pickle(filename):
    with open(filename, 'rb') as f:
        data = pickle.load(f)
    f.close()

    return data

