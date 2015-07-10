

def test_calc_logL():

    import numpy as np
    from mystats import calc_logL

    data = np.random.normal(1, size=10)
    model = np.random.normal(1, size=10)
    error = 1
    weights = np.linspace(0,100,10)

    calc_logL(model, data, error, weights=weights)




