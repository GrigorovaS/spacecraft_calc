import numpy as np
from scipy.special import erfc
from scipy.constants import c, pi

# Estimate the Free Space Losses

sma = 7028.0e3  # Semi-Major Axis, m
freq_hz = 8400e6  # Frequency, Hz 
earth_radius_m = 6371.0e3  # Earth's radius, m

# Calculate slant range
slr_m = sma + earth_radius_m

# Free Space Loss formula
fsl_db = 20 * np.log10(slr_m) + 20 * np.log10(freq_hz) + 20 * np.log10(4 * pi / c)

print(f"Free Space Loss: {fsl_db:.2f} dB")

# Assumption based on ITU-R recommendations for X-band
# and typical values for Northern Latitudes
atmospheric_loss_db = 3.0 

print(f"Atmospheric Losses: {atmospheric_loss_db} dB")

# Link parameters
eirp_dbw = 3  # EIRP, dBW
gt_dbk = 27  # G/T of the ground station, dB/K
link_bw_hz = 30e6  # Link Bandwidth, Hz

# Calculate C/N ratio from the link budget
cn_ratio_db = eirp_dbw + gt_dbk - fsl_db - atmospheric_loss_db - 10 * np.log10(link_bw_hz)

# Convert C/N ratio from dB to linear scale for BER calculation
cn_ratio_linear = 10 ** (cn_ratio_db / 10)

# Calculate Eb/N0 for QPSK
# Assuming the spectral efficiency of QPSK is 2 bits/s/Hz, then Eb/N0 = (C/N) / 2
eb_n0_linear = cn_ratio_linear / 2

# Estimate the Bit Error Rate for QPSK
ber = 0.5 * erfc(np.sqrt(eb_n0_linear))

print(f"BER (QPSK): {ber:.2e}")

