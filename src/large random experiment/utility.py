import random
import math
import time
import os
from scipy.spatial import KDTree

random.seed(time.time())

def haversine(lat1, lon1,  lat2, lon2):
    """
    from http://stackoverflow.com/questions/4913349/
    Calculate the great circle distance in kilometers between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(math.radians, [lon1, lat1, lon2, lat2])

    # haversine formula 
    dlon = lon2 - lon1 
    dlat = lat2 - lat1 
    a = math.sin(dlat/2)**2 + math.cos(lat1) * math.cos(lat2) * math.sin(dlon/2)**2
    c = 2 * math.asin(math.sqrt(a)) 
    R = 3956 # Radius of earth in miles. Use 6371 for Km. Determines return value units.
    return c * R

# Change this function as required
def get_distance(p_1, p_2):
    dist = haversine(p_1.latitude, p_1.longitude, p_2.latitude, p_2.longitude)
    return dist

def get_travel_energy(p_1, p_2, m):
    """
    Energy requirement for travel in kWh
    """
    dist = get_distance(p_1, p_2)
    energy = dist / m
    return energy


def get_travel_time(p_1, p_2, v):
    """
    Travel time in min
    """
    dist = get_distance(p_1, p_2)
    energy = dist / v
    return energy



def create_loc(central_lat, central_lon, radius):
    """
    Creates a random location (lat,lon) in the radius of radius
    from a central location

    """
    R = 3956 # Radius of earth in miles. Use 6371 for Km. Determines return value units.
    distance = random.uniform(0, radius)  # Random distance between 0 and radius_km
    bearing = random.uniform(0, 2 * math.pi)  # Random bearing (angle in radians) from 0 to 2pi

    # Calculate the new latitude and longitude
    lat1 = math.radians(central_lat)
    lon1 = math.radians(central_lon)
    lat2 = math.asin(math.sin(lat1) * math.cos(distance / R) +
                     math.cos(lat1) * math.sin(distance / R) * math.cos(bearing))
    lon2 = lon1 + math.atan2(math.sin(bearing) * math.sin(distance / R) * math.cos(lat1),
                             math.cos(distance / R) - math.sin(lat1) * math.sin(lat2))
    
    # Convert back to degrees
    lat2 = math.degrees(lat2)
    lon2 = math.degrees(lon2)
    
    return lat2, lon2


def create_kd_tree(cp, C, theta, min_lat, max_lat, min_lon, max_lon, min_theta, max_theta, denom_lat, denom_lon, denom_theta):
    """
    Create a KD-Tree from a list of Fog objects
    """
    # lats = [c.latitude for c in cp[1:]] # first cp is the depot
    # lons = [c.longitude for c in cp[1:]]
    # min_lat = min(lats)
    # min_lon = min(lons)
    # min_theta = min(theta[1:])
    # max_lat = max(lats)
    # max_lon = max(lons)
    # max_theta = max(theta[1:])
    # denom_lat = max_lat - min_lat
    # denom_lon = max_lon - min_lon
    # denom_theta = max_theta - min_theta
    points = []
    for i in range(1, len(C)):
        norm_lat = (cp[i].latitude - min_lat) / denom_lat
        norm_lon = (cp[i].longitude - min_lon) / denom_lon
        norm_theta = (theta[i] - min_theta) / denom_theta
        points.append((norm_lat, norm_lon, norm_theta))
    kd_tree = KDTree(points)
    return kd_tree

def find_nearest_static_fog(kd_tree, loc_lat, loc_lon, min_lat, max_lat, min_lon, max_lon, min_theta, max_theta, denom_lat, denom_lon, denom_theta):
    """
    Function to find the nearest fog of a waypoint using the KD-Tree.

    Parameters
    ----------
    kd_tree : TYPE
    static_fogs : TYPE
        A list of static fogs on which the KD-tree is built.
    loc_lat : TYPE
        Latitude of the waypoint.
    loc_lon : TYPE
        Longitube of the waypoint.

    Returns
    -------
    nearest_fog : TYPE
        Nearest Fog of the waypoint.
    haversine_distance : TYPE
        Distance from the waypoint to the fog in Km.

    """
    # Query the KD-Tree for the nearest neighbor to the waypoint
    norm_lat = (loc_lat - min_lat) / denom_lat
    norm_lon = (loc_lon - min_lon) / denom_lon
    norm_theta = 0 # (min_theta - min_theta) / denom_theta # we try to choose the cp with lowest theta
    point = (norm_lat, norm_lon, norm_theta)
    _, index = kd_tree.query(point)

    return index+1 # As while creating the KD tree we don't take the first cp, i.e. the depot