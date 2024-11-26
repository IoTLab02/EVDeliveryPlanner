
import pandas as pd
import math
import os
import time 
import random
from st_dbscan import ST_DBSCAN
from collections import defaultdict
from utility import create_kd_tree, find_nearest_static_fog
from delivery import Delivery
from cp import CP
from ev import EV


# random.seed(time.time())

DEBUG = False
def printDebug(string, val):
    if DEBUG == True:
        print(string, val)


# Read the CSV file
csv_file = 'combined_data.csv'  # Replace with your CSV file path
data = pd.read_csv(csv_file)
# data = data1.iloc[::-1].reset_index(drop=True)


#user input
row_max = 230 # how many rows to take from the csv
time_slot_max = 60 # max deadline for all deliveries

# Number of  (We can control it as user input)
nW = 1  # Warehouse depots
# nD = len(deliveries)  # Delivery points
nD = 25  # max value 20
#nC = len(cp) - 1 # Charging points 
nC = 2 # max total cp including warehouse should be 10. nW + nC <=0
nE = 2  # Electric delivery trucks
nS = math.ceil(nD / nE)  # max Subtrips per EV

# Constants
mj = 2 # mileage in miles per kWh
# vj = 0.5 # EV velocity in miles per min. 30 mph -> 0.5 miles per min
beta_f = 30 # EV full capacity 60 kWh
M = 99999 # very large M for Big M constraint



# objects of Delivery and CP classes
deliveries = []
cp = []
EVs = []

# Maps to track local ID to global ID and global ID to local ID
delivery_map_local_to_global = {}
delivery_map_global_to_local = {}
cp_map_local_to_global = {}
cp_map_global_to_local = {}

# Iterate over the DataFrame and create objects based on the type column
for index, row in data.iterrows():
    # if index == row_max:
    #     break
    location = eval(row['location'])  # Convert the string representation of the dict to an actual dict
    lat = location['latitude']
    long = location['longitude']
    global_id = row['index']  # Assuming 'index' column is the global ID

    # Check the type column and create the appropriate object
    if row['type'] == 'restaurant':
        local_id = len(deliveries)  # Local ID is the index in the deliveries array
        deliveries.append(Delivery(lat, long, local_id, global_id))
        delivery_map_local_to_global[local_id] = global_id
        delivery_map_global_to_local[global_id] = local_id
    elif row['type'] == 'gas_station':
        local_id = len(cp)  # Local ID is the index in the cp array
        cp.append(CP(lat, long, local_id, global_id))
        cp_map_local_to_global[local_id] = global_id
        cp_map_global_to_local[global_id] = local_id


# List of entities
W = {0} # Warehouse
D = set(range(0, nD)) 
C = set(range(1, nC + 1)) # CP starts from 1 as W is also considered as a CP and with id 0
C1 = W.union(C)
E = set(range(0, nE))
S = set(range(0, nS))

# We trim and use the list as per user input
deliveries = deliveries[:nD]
cp = cp[:(nW + nC)]


# create EV object
for j in E:
    EVs.append(EV(j, cp[0], beta_f))


# Initialize window time for every delivery point
random.seed(123)
tau_start = {}
tau_end = {}
for y in D:
    tau_start[y] = random.randint(10, time_slot_max) # 10 min to 480 min
    tau_end[y] = tau_start[y] + 30


# Per kWh energy cost at a charging station
random.seed(123)
theta = [round(random.uniform(0.4, 0.6),2) for _ in C1]
theta[0] = 0.36

dist_DD = {} # Distance between any two delivery points
dist_DC = {} # Distance between a delivery point and a charging point
psi_DD = {} # Energy requirement to travel between two delivery points
psi_DC = {} # Energy requirement to travel between a delivery point and a charging point
gamma_DD = {} # Time to reach from one delivery point to another
gamma_DC = {} # Time to reach from one delivery point to charging point or opposite


file_path_diatance = 'distance_matrix.csv'
matrix = pd.read_csv(file_path_diatance, index_col=0)
dist = {}
for i in matrix.index[:]:
    for j in matrix.columns[:]:
        dist[(int(i), int(j))] = matrix.loc[i, j] * 1.5 # we consider the distance is 1.5 times
        
        
file_path_duration = 'complete_duration_matrix2.csv'
matrix_time = pd.read_csv(file_path_duration, index_col=0)
time1 = {}
for i in matrix_time.index[:]:
    for j in matrix_time.columns[:]:
        time1[(int(i), int(j))] = round(matrix_time.loc[i, j] / 60, 2) # in min



# Store energy and time requirement for traveling each edge
for d1 in D:
    for d2 in D:
        d1_global = delivery_map_local_to_global[d1]
        d2_global = delivery_map_local_to_global[d2]
        # dist12 = random.randint(1, 2 * radius_mile) # in miles
        dist12 = dist[(d1_global, d2_global)]
        # dist_DD[(d1,d2)] = dist12
        # dist_DD[(d2,d1)] = dist12
        energy12 = round(dist12 /mj, 2) # in kWh
        psi_DD[(d1,d2)] = energy12
        # time12 = round(dist12 /vj, 2) # in min
        time12 = time1[(d1_global, d2_global)]
        gamma_DD[(d1,d2)] = time12
    for c in C1:
        d1_global = delivery_map_local_to_global[d1]
        c_global = cp_map_local_to_global[c]
        dist1c = dist[(d1_global, c_global)]
        # dist1c = random.randint(1, 2 * radius_mile) # in miles
        # dist_DC[(d1,c)] = dist1c
        # dist_DC[(c,d1)] = dist1c
        energy1c = round(dist1c /mj, 2) # in kWh
        psi_DC[(d1,c)] = energy1c
        # time12 = round(dist1c /vj, 2) # in min
        time12 = time1[(d1_global, c_global)]
        gamma_DC[(d1,c)] = time12
        
rateE = {} # charge acceptance rate at EV
random.seed(123)

for j in E:
    rateE[j] = random.randint(90, 350) # in kW per min
    
rateC = {} # charging rate at CP
for y in C:
    rateC[y] = random.randint(50, 350) # in kW per min



start_time = time.time() # Algorithm start time

# Create the KD-Tree
lats = [c.latitude for c in cp[1:]] + [d.latitude for d in deliveries] # first cp is the depot
lons = [c.longitude for c in cp[1:]] + [d.longitude for d in deliveries]
min_lat = min(lats)
min_lon = min(lons)
min_theta = min(theta[1:])
max_lat = max(lats)
max_lon = max(lons)
max_theta = max(theta[1:])
denom_lat = max_lat - min_lat
denom_lon = max_lon - min_lon
denom_theta = max_theta - min_theta
kd_tree = create_kd_tree(cp, C, theta, min_lat, max_lat, min_lon, max_lon, min_theta, max_theta, denom_lat, denom_lon, denom_theta)

# Find the nearest fog for each delivery and 
# prepare a list of spatio-temporal parameters of the deliveries for ST-DBSCAN 
deliveries_spatiotemporal_params = []
min_deadline = min([tau_end[y] for y in D])
max_deadline = max([tau_end[y] for y in D])
denom_deadline = max_deadline - min_deadline

for d in deliveries:
    nearest_fog_index = find_nearest_static_fog(kd_tree, d.latitude, d.longitude, min_lat, max_lat, min_lon, max_lon, min_theta, max_theta, denom_lat, denom_lon, denom_theta)
    d.nearest_CP_local_id = nearest_fog_index
    local_id_d = d.local_id
    # create a list of waypoints with the deadlines (spatio-temporal parameters)
    param = [(tau_end[local_id_d]-min_deadline)/denom_deadline, (d.latitude-min_lat)/denom_lat, (d.longitude-min_lon)/denom_lon]
    deliveries_spatiotemporal_params.append(param)
    

# Cluster the deliveries depending on their spatio-temporal parameters
# Hence eps1, that is the spatial density threshold (maximum spatial distance) between 
# 
# eps2 is the temporal threshold (maximum temporal distance) between two 
# points to be considered related.     
## Avg Euclidean dist among the fog lat lon values: 0.0037424493333129533

# tune eps1 and eps2 depending on the values of  deliveries_spatiotemporal_params
st_dbscan = ST_DBSCAN(eps1 = 0.2, eps2 = 0.2, min_samples = 1) 
st_dbscan.fit(deliveries_spatiotemporal_params) 
labels = st_dbscan.labels
clusters = defaultdict(list)
clusters_min_deadline = {}
for i in D:
    delivery = deliveries[i]
    label = labels[i]
    clusters[label].append(delivery)

# Check the output and tune parameters eps1, eps2, and min_samples
for cluster_label, incluster_deliveries in clusters.items():
    clusters_min_deadline[cluster_label] = min([tau_end[d.local_id] for d in incluster_deliveries])
    # print(f"Cluster {cluster_label}:")
    # for delivery in incluster_deliveries:
    #     print(f"  {delivery.local_id}")
        
# Sort the cluster ids depending on their min deadlines
sorted_cluster_ids = sorted(clusters_min_deadline, key=lambda k: clusters_min_deadline[k])


# Assign deliveries to the EVs
max_subtrip_per_EV = nS
for c_id in sorted_cluster_ids:
    # sort the deliveries in each cluster in the descending order of their deadline
    sorted_deliveries = sorted(clusters[c_id], key=lambda d: tau_end[d.local_id])
    
    for delivery in sorted_deliveries:
        d_id = delivery.local_id
        #print("cluster:", c_id, " delivery: ", d_id, "delivery start time: ", tau_start[d_id])
        
        energy_cost_min = 9999999 # cost is measured in terms of energy used
        assigned_EV = None
        time_req_for_delivery = -1
        for e in EVs:
            if e.totalSubtrip >= max_subtrip_per_EV:
                printDebug("Max delivery reached. EV: ", e.local_id)
                continue
            last_loc = e.A[-1]
            duration_lastLoc2delivery = 0
            energy_req = 0 # energy required to reach the current delivery point from the last point
            energy_req_d2W = 0 # energy required to reach warehouse from the current delivery point
            if last_loc.type == "delivery":
                duration_lastLoc2delivery = gamma_DD[(last_loc.local_id, d_id)]
                energy_req = psi_DD[(last_loc.local_id, d_id)]
                
            else:
                duration_lastLoc2delivery = gamma_DC[(d_id, last_loc.local_id)]
                energy_req = psi_DC[(d_id, last_loc.local_id)]
            
            energy_req_d2W = psi_DC[(d_id, 0)]
            nearest_cp = cp[delivery.nearest_CP_local_id]
            
            if energy_req < energy_cost_min:
                # Check the battery constraints 
                if energy_req + energy_req_d2W > e.B_res:
                    if last_loc.type == "delivery": # if already at a CP then charge is already full
                        # Send it to the nearest CP for charging
                        nearest_cp_id_from_last_point = last_loc.nearest_CP_local_id
                        e_id = e.local_id
                        e.A.append(cp[nearest_cp_id_from_last_point])
                        req_charge = beta_f - e.B_res - psi_DC[(last_loc.local_id, nearest_cp_id_from_last_point)]
                        e.T.append(e.t_ed)
                        e.t_ed = gamma_DC[(last_loc.local_id, nearest_cp_id_from_last_point)] +  req_charge/min(rateE[e_id], rateC[nearest_cp_id_from_last_point]) # reaching time + charging time
                        e.B_res = beta_f # reset battery
                        e.B_trace.append(beta_f)
                    continue
                
                reaching_time = e.t_ed + duration_lastLoc2delivery
                if reaching_time <= tau_end[d_id]:
                    assigned_EV = e
                    energy_cost_min = energy_req
                    time_req_for_delivery = duration_lastLoc2delivery
        if assigned_EV != None:
            delivery.assigned_EV = assigned_EV
            assigned_EV.assigneddeliveries.append(delivery.local_id)
            # last_point = assigned_EV.A[-1]
            # time_flying = getFlyingTime(last_point, delivery, p.horizontalSpeed) # time to fly from last loc to delivery waypoint
            t_ed = assigned_EV.t_ed # earliest departure time from the last location
            #t_dep = max((delivery.latest_sensing_start_timeinstance - time_flying), t_ed) # actual departure time from the last location
            t_dep = t_ed # As there is no sensing start time
            assigned_EV.A.append(delivery)
            assigned_EV.T.append(t_dep)
            #assigned_EV.t_ed = max((t_dep +  time_req_for_delivery),  tau_end[d_id]) # old version
            assigned_EV.t_ed = min(max((t_dep +  time_req_for_delivery), tau_start[d_id]), tau_end[d_id])# earliest departure time at the delivery waypoint ## AK: New change
            #assigned_EV.t_ed = (t_dep +  time_req_for_delivery) # if only deadline exists. no slot start time
            # assigned_EV.totalEnergyUsed = assigned_EV.totalEnergyUsed + energy_req_for_delivery
            B_res = assigned_EV.B_res - energy_cost_min
            assigned_EV.B_res = B_res
            assigned_EV.B_trace.append(B_res)
            
            
            printDebug("energy_cost_min: ", energy_cost_min)


####            
unassigned_deliveries = [t for t in deliveries if t.assigned_EV==None]











for e in EVs:
    last_loc = e.A[-1]
    if last_loc.type == "delivery":
        e.A.append(cp[0]) # return to the warehouse
    if last_loc.type == "cp":
        if last_loc.local_id != 0:
            A1 = e.A[:-1] # remove the last CP
            A1.append(cp[0]) # Add the depot as last charging point
            e.A = A1



end_time = time.time()
elapsed_time = (end_time - start_time) * 1000


# result variables    
total_deliveries_completed = 0
total_cost = 0
total_energy = 0
delivery_distribution = {}
for e in EVs:
    A = e.A
    # total_cost_e = 0
    energy_per_trip = 0
    delivery_done = 0
    last_loc = A[0]
    # print("Trace::")
    # print(last_loc.local_id)
    for loc in A[1:]:
        energy_req = 0
        if loc.type == "delivery":
            delivery_done += 1
            if last_loc.type == "delivery":
                energy_req = psi_DD[(last_loc.local_id, loc.local_id)]
            else:
                energy_req = psi_DC[(loc.local_id, last_loc.local_id)]
            energy_per_trip += energy_req
            last_loc = loc
            # print("->",loc.local_id, " energy: ", energy_req)
        else:
            if last_loc.type == "delivery":
                energy_req = psi_DC[(last_loc.local_id, loc.local_id)]
            else:
                print("ERROR: coming to a cp from a cp")
            # print("->",loc.local_id, " energy: ", energy_req)
            energy_per_trip += energy_req
            total_energy += energy_per_trip
            total_cost += energy_per_trip * theta[loc.local_id]
            last_loc = loc
            energy_per_trip = 0
    total_deliveries_completed += delivery_done
    delivery_distribution[e.local_id] = delivery_done

print("****Heuristic****")          
print("Total delivery: ", nD, " Total CP: ", nC+1, " Total EV: ", nE, " time_slot_max:", time_slot_max)          
print(f"Elapsed time: {elapsed_time:.2f} ms")
print(f"Total cost: {total_cost}")
print(f"Total successful delivery: {total_deliveries_completed}")
print("::Delivery Distribution::")
for j in E:
    print(f"delivery_distribution[{j}]:", delivery_distribution.get(j, 0))
            
    



        
