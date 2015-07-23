
def test_calc_symmetric_error():

    import numpy as np
    from mystats import calc_symmetric_error
    import matplotlib.pyplot as plt

    x = np.linspace(-5,15,1000)
    y = np.exp(-(x-3)**2 / (2 * 1**2)) + \
        3 * np.exp(-(x-9)**2 / (2 * 0.5**2)) + \
        0.1 * np.exp(-(x+0.5)**2 / (2 * 0.5**2))

    x_clip = x[x < 0]
    y[x < 0] = 0.1 * x_clip + 2.5 / 5

    median, high_error, low_error = calc_symmetric_error(x, y, alpha=1.0 - 0.68)

    plt.plot(x, y)
    plt.axvspan(median - low_error, median + high_error, alpha=0.5)
    plt.axvline(median, linestyle='--')
    plt.savefig('test_plots/test_calc_symmetric_error.png')





