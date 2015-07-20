

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
