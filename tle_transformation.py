import numpy as np 
import math
from datetime import datetime, timedelta

G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
M = 5.97219e24   # mass of the Earth, kg


def extract_tle_time(tle_line1):
    # Extract epoch year and day
    epoch_year_short = int(tle_line1[18:20])
    epoch_day = float(tle_line1[20:32])

    # Convert short year to full year
    epoch_year = 2000 + epoch_year_short

    # Convert day of the year to month and day
    date = datetime(epoch_year, 1, 1) + timedelta(days=int(epoch_day) - 1)
    month = date.month
    day_of_month = date.day

    # Convert fractional day to hours, minutes, and seconds
    fractional_day = epoch_day - int(epoch_day)
    hours = int(fractional_day * 24)
    minutes = int((fractional_day * 24 - hours) * 60)
    seconds = int((((fractional_day * 24 - hours) * 60) - minutes) * 60)

    return epoch_year, month, day_of_month, hours, minutes, seconds


def tle_to_keplerian(tle_line1, tle_line2):
    # Extract parameters from TLE lines
    inc = float(tle_line2[8:16])  # Inclination [degrees]
    raan = float(tle_line2[17:25])  # Right Ascension of Ascending Node [degrees]
    ecc = float("0." + tle_line2[26:33])  # Eccentricity
    argp = float(tle_line2[34:42])  # Argument of Perigee [degrees]
    ma = float(tle_line2[43:51])  # Mean Anomaly [degrees]
    n = float(tle_line2[52:63])  # Mean Motion [revs per day]

    # Convert Mean Motion from revs per day to radians per second
    n_rad_per_sec = n * 2 * math.pi / (24 * 3600)

    # Calculate semi-major axis from mean motion
    a = ( G * M / (n_rad_per_sec ** 2)) ** (1. / 3.)  # semi-major axis, m

    # Convert angles from degrees to radians
    inc_rad = math.radians(inc)
    raan_rad = math.radians(raan)
    argp_rad = math.radians(argp)
    M_rad = math.radians(ma)

    return (a, ecc, inc_rad, raan_rad, argp_rad, M_rad)




def solve_keplers_equation(M_rad, e, tolerance=1e-6):
    E = M_rad  # Initial guess for Eccentric Anomaly
    while True:
        delta = E - e * np.sin(E) - M_rad
        if abs(delta) < tolerance:
            break
        E -= delta / (1 - e * np.cos(E))
    return E

def keplerian_to_cartesian(a, e, inc_rad, raan_rad, argp_rad, M_rad):

    # Solve for Eccentric Anomaly
    E = solve_keplers_equation(M_rad, e)

    # True Anomaly
    ν = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2), np.sqrt(1 - e) * np.cos(E / 2))

    # Distance
    r = a * (1 - e * np.cos(E))

    # Position in the orbital plane
    x_op = r * np.cos(ν)
    y_op = r * np.sin(ν)

    # Velocity in the orbital plane
    vx_op = -np.sqrt( G * M / a) * np.sin(E) / (1 - e * np.cos(E))
    vy_op = np.sqrt( G * M / a) * (np.cos(E) - e) / (1 - e * np.cos(E))

    # Rotation matrices
    cos_argp = np.cos(argp_rad)
    sin_argp = np.sin(argp_rad)
    cos_inc = np.cos(inc_rad)
    sin_inc = np.sin(inc_rad)
    cos_raan = np.cos(raan_rad)
    sin_raan = np.sin(raan_rad)

    # Position in ECI coordinates
    x = x_op * (cos_argp * cos_raan - sin_argp * cos_inc * sin_raan) - y_op * (sin_argp * cos_raan + cos_argp * cos_inc * sin_raan)
    y = x_op * (cos_argp * sin_raan + sin_argp * cos_inc * cos_raan) + y_op * (cos_argp * cos_inc * cos_raan - sin_argp * sin_raan)
    z = x_op * sin_argp * sin_inc + y_op * cos_argp * sin_inc

    # Velocity in ECI coordinates
    vx = vx_op * (cos_argp * cos_raan - sin_argp * cos_inc * sin_raan) - vy_op * (sin_argp * cos_raan + cos_argp * cos_inc * sin_raan)
    vy = vx_op * (cos_argp * sin_raan + sin_argp * cos_inc * cos_raan) + vy_op * (cos_argp * cos_inc * cos_raan - sin_argp * sin_raan)
    vz = vx_op * sin_argp * sin_inc + vy_op * cos_argp * sin_inc

    position = np.array([x, y, z])  # Position, m
    velocity = np.array([vx, vy, vz])  # Velocity, m/s

    return position, velocity 


#tle_line1 = "1 25544U 98067A   21053.17708333  .00001295  00000-0  29636-4 0  9993"
#tle_line2 = "2 25544  51.6445  47.1467 0002776  45.4628  68.8661 15.48915128202379"
#
#
## Parse TLE and calculate Keplerian elements
#keplerian_elements = tle_to_keplerian(tle_line1, tle_line2)
#
#print("Keplerian Elements:")
#print("Semi-major axis (km):", keplerian_elements[0])
#print("Eccentricity:", keplerian_elements[1])
#print("Inclination (rad):", keplerian_elements[2])
#print("RAAN (rad):", keplerian_elements[3])
#print("Argument of Perigee (rad):", keplerian_elements[4])
#print("Mean Anomaly (rad):", keplerian_elements[5])
#
## Convert to Cartesian coordinates
#position, velocity = keplerian_to_cartesian(*keplerian_elements)
#
#print("Position Vector (m):", position)
#print("Velocity Vector (m/s):", velocity)

