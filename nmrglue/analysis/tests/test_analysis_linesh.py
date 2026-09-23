import numpy as np
from nmrglue.analysis.linesh import fit_spectrum, sim_NDregion

def test_fit_spectrum():
    lineshapes = ['g']
    params = [[(13797.0, 2.2495075273313034)],
              [(38979.0, 5.8705185693227664)],
              [(39066.0, 5.7125954296137103)],
              [(39153.0, 5.7791485451283791)],
              [(41649.0, 4.260242375400459)],
              [(49007.0, 4.2683625950679964)],
              [(54774.0, 3.2907139764685569)]]
    amps = [35083.008667, 32493.824402, 32716.156556, 33310.711914, 82682.928405,
            82876.544313, 85355.658142]
    bounds = [[[(None, None), (0, None)]], [[(None, None), (0, None)]],
              [[(None, None), (0, None)]], [[(None, None), (0, None)]],
              [[(None, None), (0, None)]], [[(None, None), (0, None)]],
              [[(None, None), (0, None)]]]
    ampbounds = [None, None, None, None, None, None, None]
    centers = [(13797.0,), (38979.0,), (39066.0,), (39153.0,), (41649.0,),
               (49007.0,), (54774.0,)]
    rIDs = [1, 2, 3, 4, 5, 6, 7]
    box_width = (5,)

    generated_data = sim_NDregion(shape=(65536,), lineshapes=['g'], params=params, amps=amps)

    params_best, amp_best, iers = fit_spectrum(
        generated_data, lineshapes, params, amps, bounds, ampbounds, centers,
        rIDs, box_width, error_flag=False, verb=False)
    
    assert np.allclose(params_best, params, rtol=1e-4)
    assert np.allclose(amp_best, amps, rtol=1e-4)
