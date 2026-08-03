"""
Denoising using non-ML and ML algorithms.
Using wrapper functions for existing data within the NMR denoiser toolbox, 
as well as in-house denoising models.
"""

import numpy as np
import scipy.optimize
import nmrdenoise as nd
from nmrglue.fileio.bruker import scale_pdata


def denoise(nmr_file=None, data=None, model="Desperate"):
    """
    Computes denoised data across all models (ML and non-ML).
    
    Example use: denoise(nmr_file=None, data=data, model="Desperate")

    Parameters:
    -----------
    data : array_like (2D array real + imaginary)
        The input NMR data after FFT has been applied by nmrglue.
    model : str, optional
        The algorithm to use for denoising. Defaults to "Desperate".
    nmr_file : str, optional
        Path to the NMR file (e.g., Bruker format), if applicable.
        Currently only bruker format is supported for converting x-axis to parts-per-million (ppm) extraction.
        
    Returns:
    --------
    Tuple or ndarray
        Denoised data, and frequency axis if nmr_file is provided.
    """
    T = np.size(data)  # Total duration in seconds

    multiple = 128
    new_length = (len(data) // multiple) * multiple

    cropped_data = data[:new_length]
    og_data = data[:new_length]

    # LD-net requires data normalization as it's not in-built
    data_min = cropped_data.min()
    data_max = cropped_data.max()
    normalized_data = (cropped_data - data_min) / (data_max - data_min)

    fn = len(normalized_data)
    cropped_data_reform = normalized_data.reshape(1, fn)

    if model == "AEdenoise":
        """
        @OG Repository Author: Kinkini Monaragala, Built into nmr denoise as a piece of original work
        Paper DOI: tba (coming soon)
        """
        ae = nd.AE_denoise()
        autoencoder_path = ae.load_denoising_checkpoint(
            model_name="autoencoder_20000_bruker_ver_135_percent_noise_with_noshuffle_v4_rm_spikelets_with_lb_5.pt"
        )
        denoised_data, confidence = ae.denoise_component(cropped_data, autoencoder_path)

    elif model == "Desperate":
        """
        @OG Repository Author: Haolin Zhan
        Link to repository: https://github.com/rschurko/DESPERATE/tree/master
        DOI: https://doi.org/10.1016/j.jmr.2022.107320
        """
        dn = nd.Desperate()
        denoised_data, coeffin, coeffs = dn.wavelet_denoise(
            7, cropped_data, 0, wave='bior2.2', threshold='mod', alpha=0
        )

    elif model == "LDNet":
        """
        @OG Repository Authors: Harris Mason and Adam Altenhof
        Link to repository: https://github.com/HaolinZhan/LD-Net-for-NMR-spectroscopy-denoising
        """
        ld = nd.LD_net()
        denoised_data = ld.test(cropped_data_reform, fn)

    # Denormalize output
    # Can make this entirely optional if preferred 
    #this converts data from 0 to 1 scale to original scale
    denoised_data = denoised_data * (data_max - data_min) + data_min

    min_dat_length = min(len(og_data), len(denoised_data))

    # if bruker file is present (can extend this to other formats as well)
    if nmr_file:
        try:
            ae_reader = nd.AE_denoise()

            # converts the data to ppm (matches software such as mesternova following it's metadata exactly!)
            frq, _, _ = ae_reader.read_bruker(nmr_file)
            return frq[:min_dat_length], denoised_data[:min_dat_length]
        except Exception as e:
            print(f"Error reading Bruker file: {e}")

    return denoised_data[:min_dat_length]

