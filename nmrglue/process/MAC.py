from numpy import *

def bruk_ranceY(dic, data, rotate_phase=True):
    """
    Create a unit conversion object which accepts/returns NMRPipe indices.

    Parameters
    ----------
    dic : dict
        Dictionary of NMRPipe parameters
    data : ndarray
        Array of NMR data.
    rotate_phase : bool, optional
        Remove the requirement for a 90 degree zero-order phase correction

     Returns
    -------
    ndic : dict
        Dictionary of updated NMRPipe parameters.
    ndata : ndarray
        Array of NMR data which has been reshuffled according to the Rance-Kay scheme.

    """
    shuffled_data = np.empty(data.shape, data.dtype)

    for i in range(0, data.shape[0], 2):
        shuffled_data[i] = (1.*(data[i].real - data[i+1].real) +
                            1.*(data[i].imag - data[i+1].imag)*1j)
        if rotate_phase == True:
            shuffled_data[i+1] = (-1.*(data[i].imag + data[i+1].imag) +
                                   1.*(data[i].real + data[i+1].real)*1j)
        else:
            shuffled_data[i+1] = (1.*(data[i].real + data[i+1].real) +
                                  1.*(data[i].imag + data[i+1].imag)*1j)

    return dic, shuffled_data

