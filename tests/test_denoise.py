import nmrglue as ng
import numpy as np



def create_noise(data, snr: int, seed=None):

    """
    Injects white noise to create desired snr based on zhangs peak snr calcualtion    

    Parameters:
    -----------
    data : 
        The input NMR spectrum, it has to be real part.
    snr : int
        The target SNR to create artificially

    Returns:
    --------
       The data with noise injected to meet snr criterion
    """

    signal_max = np.max(np.real(data))
    target_noise_std = signal_max / snr
    
    rng = np.random.default_rng(seed=seed)
    noise = rng.normal(0, target_noise_std, data.shape)
    
    noise = noise - np.mean(noise)
    
    data_with_noise = data + noise
    return data_with_noise

def normalize_spectrum(spec):
    """
    Normalizes the spectrum using max peak and baseline median    

    Parameters:
    -----------
    spec : array
        The input NMR spectrum.

    Returns:
    --------
    Normalized NMR spectrum : array
    """
    baseline = np.median(spec)
    shifted = spec - baseline
    peak_max = np.max(shifted)
    
    if peak_max <= 0:
        return shifted
        
    return shifted / peak_max

def generate_synthetic_nmr_fid(n_points=2048, sw=5000.0, freqs=None, amps=None,t2=None, seed=None):
    """
    Generates a simulated 1D NMR FID using a sum of decaying peaks.

    Parameters:
    -----------
    n_points : int, optional
        The number of data points in time-domain FID. 
    sw : float, optional
        The spectral width in Hz. 
    freqs : list of float, optional
        Resonance frequencies (in Hz). 
    amps : list of float, optional
        Initial amplitudes for each peak. 
    t2 : list of float, optional
        Spin-to-spin relaxation time constants controlling peak widths. 
    seed : int, optional
        Random seed for reproducibility.

    Returns:
    --------
    fid : array
        Complex valued 'synthetic' FID.
    """

    #  defaults if none provided
    if freqs is None: freqs = [500.0, 1200.0, -800.0]
    if amps is None: amps = [0.5, 0.4, 0.3]
    if t2 is None: t2 = [0.5, 0.4, 0.6]
    
    t = np.linspace(0, n_points / sw, n_points, endpoint=False)
    fid = np.zeros(n_points, dtype=complex)

    # equation inspired and taken from https://github.com/DavideCandoli/PULSEE/blob/master/Code/Simulation.py
    # lines 675 to 677
    #    for t in times:
    #    dm_t = dm.free_evolution(h_unperturbed, t)
    #    FID.append((dm_t*I_plus_rotated*np.exp(-t/T2)*np.exp(-1j*2*math.pi*reference_frequency*t)).trace())

    # looping through all peaks in array
    for f, a, t2_val in zip(freqs, amps, t2):
        # this is creating a complex valued sum of decaying peaks using exponentional
        # adding new peak to to our fid with added dampening
        fid += a * np.exp(-t / t2_val) * np.exp(1j * 2 * np.pi * f * t)
        
    return fid

def denoise_synthetic_data(noise_snr=20, model="Desperate", seed=None, sw=None, amps=None, freqs=None, t2=None):
    """
    Generates synthetic NMR data, injects noise to synthetic data and applies a chosen denoising model for comparative evaluation.

    Parameters:
    -----------
    noise_snr : float, optional
        The target SNR Ratio for the noise injection. 
    model : str, optional
        The name of the denoising model to apply (e.g., "Desperate", "AEdenoise", "LDNet"). Default is "Desperate".
    seed : int, optional
        Random seed. Defaults to seed (508948) if None.
    sw : float, optional
        The spectral width in Hz for FID generator.
    amps : list of float, optional
        Peak amplitudes for FID generator.
    freqs : list of float, optional
        Peak frequencies for FID generator.
    t2 : list of float, optional
        $T_2$ relaxation for FID generator.

    Returns:
    --------
    spec_real_synthetic_norm : array
        The normalized synthetic spectrum.
    spec_real_noisy_norm : array
        The normalized noisy spectrum befoe to denoising.
    spec_denoised_norm : array
        The normalized spectrum after denoising.
    """
    # This seed shows an adequete noise floor to show denoising ability across models
    if seed is None:
        seed = 508948
        
    # Create a synthetic time-domain data
    fid_synthetic = generate_synthetic_nmr_fid()         

    # add noise with desired noise snr level
    fid_noisy = create_noise(fid_synthetic, noise_snr, seed)

    # Please avoid aggressive apodization; the models do not respond well when this setting is applied heavily.
    # do note however, that this is syntehtically generated data hence not realistic to standrad lb values set in real nmr data
    fid_apodized_noisy = ng.proc_base.em(fid_noisy, lb=0.0002)
    # standard fid processing with nmr glue applied to our synthetic data
    spec_complex_noisy = ng.proc_base.fft(fid_apodized_noisy)
    # for our example demonstration we will have imaginary component disgarded       
    spec_real_noisy = ng.proc_base.di(spec_complex_noisy)           
    

    # apply same steps to synthetic data prior to noise injection to have a point of comparison
    fid_apodized_synthetic = ng.proc_base.em(fid_synthetic, lb=0.0002)
    spec_complex_synthetic = ng.proc_base.fft(fid_apodized_synthetic)
    spec_real_synthetic = ng.proc_base.di(spec_complex_synthetic)  

    # Apply denoising with desired model
    spec_denoised = ng.proc_denoise.denoise(nmr_file=None, data=spec_real_noisy, model=model)
    
    spec_real_synthetic_norm = normalize_spectrum(spec_real_synthetic)
    spec_real_noisy_norm = normalize_spectrum(spec_real_noisy)
    spec_denoised_norm = normalize_spectrum(spec_denoised)

    return spec_real_synthetic_norm, spec_real_noisy_norm, spec_denoised_norm

#example usage: spec_real_synthetic, spec_real_noisy, spec_denoised = denoise_synthetic_data(noise_snr=0.7, model="Desperate", seed=None)
