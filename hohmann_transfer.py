import numpy as np 

mass = 10 # kg
force_of_engine = 0.01 # N
force_of_engine2 = 500 # N

def calculate_hohmann_transfer(final_position, final_velocity, mass, force_of_engine):
    G = 6.67430e-11  # gravitational constant, m^3 kg^-1 s^-2
    M = 5.97219e24   # mass of the Earth, kg
    r1 = np.linalg.norm(final_position)  
    print(f"Initial Orbit at Radius = {r1} m")
    r1_analytical = np.linalg.norm([-5529203, -2217254, 3399353]) # initial orbit radius, m
    v_initial_analytical = np.sqrt(G * M / r1_analytical) # initial velocity, m/s 
    print(f"Initial Analytical Orbit Radius = {r1_analytical} m")
    r2 = 6950.0 * 1000  # Final orbit radius in meters

    # Calculate the initial and final orbital speeds from the Vis-viva equation
    v_initial = final_velocity # initial velocity calculated by the RK4 propagator, m/s
    v_final = np.sqrt(G * M / r2)
    
    # Calculate delta Vs 
    delta_v1 = v_initial * (np.sqrt(2 * r2 / (r1 + r2)) - 1) 
    delta_v2 = v_final * (1 - np.sqrt(2 * r1 / (r1 + r2)))
    
    total_delta_v = abs(delta_v1) + abs(delta_v2)
    
    # Calculate the total burn duration
    burn_duration = total_delta_v * mass / force_of_engine
    
    # Calculate delta Vs analytically
    delta_v1_analytical = v_initial_analytical * (np.sqrt(2 * r2  / (r1_analytical + r2)) - 1) 
    delta_v2_analytical = v_final * (1 - np.sqrt(2 * r1_analytical / (r1_analytical + r2)))
    
    total_delta_v_analytical = abs(delta_v1_analytical) + abs(delta_v2_analytical)
    
    # Calculate the total burn duration
    burn_duration_analytical = total_delta_v_analytical * mass / force_of_engine


    return burn_duration, total_delta_v, total_delta_v_analytical


# Convert seconds to HH:MM:SS
def convert_seconds_to_hms(seconds):
    days = seconds // (3600*24)
    hours = (seconds % (3600*24)) // 3600
    minutes = (seconds % 3600) // 60
    remaining_seconds = seconds % 60
    return days, hours, minutes, remaining_seconds


