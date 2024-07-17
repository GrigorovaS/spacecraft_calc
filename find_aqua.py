from sgp4.api import Satrec, jday, SatrecArray
from fetch_tle import tle_line1, tle_line2

# TLE lines for Aqua 
tle_line1_old = "1 27424U 02022A   21033.23334491  .00000140  00000-0  10267-3 0  9998"
tle_line2_old = "2 27424  98.2046 128.9474 0002851  97.9550 262.1850 14.57113657462354"

#tle_line1_new = tle_line1_old 
#tle_line2_new = tle_line2_old 

tle_line1_new = tle_line1
tle_line2_new = tle_line2

print(tle_line1_new)
print(tle_line2_new)

# Initialize record from the TLE lines
aqua_old = Satrec.twoline2rv(tle_line1_old, tle_line2_old)
aqua_new = Satrec.twoline2rv(tle_line1_new, tle_line2_new)

# Specify the date and time we would like to know the location of Aqua
year = 2023
month = 2
day = 1
hour = 16
minute = 23
second = 5

# Convert time to Julian date
jd, fr = jday(year, month, day, hour, minute, second)


# Find the satellite's position and velocity
e_old, position_old, velocity_old = aqua_old.sgp4(jd, fr)
e_new, position_new, velocity_new = aqua_new.sgp4(jd, fr)

#e, r, v = sat_array.sgp4(jd, fr)

#print(f"Position using old TLE: {position_old}km")
#print(f"Position using new TLE: {position_new}km")
#print(r)


if __name__ == "__main__":
    e_old, position_old, velocity_old = aqua_old.sgp4(jd, fr)
    e_new, position_new, velocity_new = aqua_new.sgp4(jd, fr)
    print(f"Position using old TLE: {position_old}km")
    print(f"Position using new TLE: {position_new}km")

