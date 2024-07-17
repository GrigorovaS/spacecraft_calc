import requests

# Download data from CelesTrak
print("Downloading data from Celestrak")
text = requests.get('https://celestrak.org/NORAD/elements/gp.php?GROUP=active').text

# Split by line
lines = text.split('\r\n')

# Find AQUA by name
aqua_start_line = list(filter(lambda line: 'AQUA' in line, lines))[0]

# Get index
aqua_start_index = lines.index(aqua_start_line)
tle_line1 = lines[aqua_start_index + 1]
tle_line2 = lines[aqua_start_index + 2]

print(tle_line1)
print(tle_line2)
