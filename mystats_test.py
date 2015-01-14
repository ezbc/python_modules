
def main():

    import numpy as np
    from mystats import calc_symmetric_error

    x = np.linspace(0,10,100)
    y = np.exp(-(x-5)**2 / (2 * 2**2)) + np.exp(-(x-8)**2 / (2 * 1**2))

    median, high_error, low_error = calc_symmetric_error(x, y, alpha=1.0 - 0.68)

    print 'center', median, '+', high_error, '-', low_error

if __name__ == '__main__':
    main()

