#!/usr/bin/python

def test_calc_T():

    from myscience import wolfire03 as w03

    n = 10 # cm^-3
    T = w03.calc_T(n)

    assert T < 188 and T > 168

def test_calc_T_asarray():

    from myscience import wolfire03 as w03

    n = (10,10) # cm^-3
    T = w03.calc_T(n)

    assert all(T < 188) and all(T > 168)

def test_w03Temp_vs_pressureTemp():

    import matplotlib.pyplot as plt
    from myscience import wolfire03 as w03
    from myscience import calc_temperature
    import numpy as np

    # range of densities to compare
    n = np.logspace(0,2,100)

    # get Wolfire+03 temperatures for several rad fields
    G_0_list = [0.6,  0.8, 0.64,]
    cloud_list = ['California', 'Perseus', 'Taurus']
    temp_w03_list = []
    for G_0 in G_0_list:
        temp_w03_list.append(w03.calc_T(n, G_0=G_0))

    # get temperatures assumes regular pressure equilibrium
    temp_pressure = calc_temperature(n_H=n, calc_error=False)

    # Plot the figure
    # -----------------------------------------------------------------------------
    figure = plt.figure(figsize=(5,5))
    for i in xrange(len(G_0_list)):
        G_0 = G_0_list[i]
        cloud = cloud_list[i]
        temp_w03 = temp_w03_list[i]
        plt.plot(n, temp_w03,
                 linestyle='--',
                 linewidth='2',
                 label='W+03, $G_0$={0:.2f}, {1:s}'.format(G_0, cloud),
                 )

    plt.plot(n, temp_pressure,
             linestyle='-',
             linewidth=2,
             label='Ideal Gas Law',
             )

    plt.xlabel(r'$n$ [cm$^{-3}$]')
    plt.ylabel(r'$T$ [K]')
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='best')
    plt.savefig('test_plots/myscience_test.test_w03Temp_vs_pressureTemp.png')

