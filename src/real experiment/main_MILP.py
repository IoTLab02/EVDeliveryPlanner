import numpy as np
import gurobipy as gp
from gurobipy import GRB
import math

#from point import Point
#import math
#from cp import CP
#from ev import EV
#from utility import generate_cps, generate_deliveries, create_kd_tree, find_nearest_cp, findAvailableEVs, get_travel_energy, sendForCharging, get_travel_time
#import array
import time 
import random
import os
import pandas as pd
from delivery import Delivery
from cp import CP
from ev import EV



seed_val = 123
# random.seed(time.time())
random.seed(seed_val)

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y



# random.seed(time.time())

DEBUG = False
def printDebug(string, val):
    if DEBUG == True:
        print(string, val)


# Read the CSV file
csv_file = 'combined_data.csv'  # Replace with your CSV file path
data = pd.read_csv(csv_file)


#user input
row_max = 230 # how many rows to take from the csv
time_slot_max = 60 # max deadline for all deliveries

# Number of  (We can control it as user input)
nW = 1  # Warehouse depots
# nD = len(deliveries)  # Delivery points
nD = 5  # max value 20
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
    if index == row_max:
        break
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
random.seed(seed_val)
tau_start = {}
tau_end = {}
for y in D:
    tau_start[y] = random.randint(10, time_slot_max) # 10 min to 480 min
    tau_end[y] = tau_start[y] + 30


# Per kWh energy cost at a charging station
random.seed(seed_val)
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
for i in matrix.index[:row_max]:
    for j in matrix.columns[:row_max]:
        dist[(int(i), int(j))] = matrix.loc[i, j] * 1.5 # we consider the distance is 1.5 times
        
        
file_path_duration = 'complete_duration_matrix2.csv'
matrix_time = pd.read_csv(file_path_duration, index_col=0)
time1 = {}
for i in matrix_time.index[:row_max]:
    for j in matrix_time.columns[:row_max]:
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
random.seed(seed_val)

for j in E:
    rateE[j] = random.randint(90, 350) # in kW per min
    
rateC = {} # charging rate at CP
for y in C:
    rateC[y] = random.randint(50, 350) # in kW per min


# Tunable parameters
alpha1 = 25 # parameter for tuning number of successful deliveries
alpha2 = 5 # parameter for tuning the energy cost

total_deliveries_completed = 0
total_cost = 0
delivery_distribution = {}
elapsed_time = 0


# Create the Gurobi model
model = gp.Model("TripPlanner")





options = {
# Gurobi WLS license file
# Your credentials are private and should not be shared or copied to public repositories.
# Visit https://license.gurobi.com/manager/doc/overview for more information.
"WLSACCESSID" : "47a60875-8894-48f5-866e-942b043fb06f",
"WLSSECRET" : "0cef4017-a933-4dbd-8bec-d54378908634",
"LICENSEID" : 2558567,
}



with gp.Env(params=options) as env, gp.Model(env=env) as model:
    start_time = time.time()
    chiCD = model.addVars(C1, D, E, S, vtype=GRB.BINARY, name="chiCD") # Decision var for edge b/w C and D
    chiDC = model.addVars(D, C1, E, S, vtype=GRB.BINARY, name="chiDC")
    chiDD = model.addVars(D, D, E, S, vtype=GRB.BINARY, name="chiDD")

    TarrC = model.addVars(C1, E, S, vtype=GRB.CONTINUOUS, lb=0, ub= GRB.INFINITY, name="arrival_time_at_C")
    TarrD = model.addVars(D, vtype=GRB.CONTINUOUS, lb=0, ub= GRB.INFINITY, name="arrival_time_at_D")
    TdepC = model.addVars(C1, E, S, vtype=GRB.CONTINUOUS, lb=0, ub= GRB.INFINITY, name="departure_time_at_C")
    TdepD = model.addVars(D, vtype=GRB.CONTINUOUS, lb=0, ub= GRB.INFINITY, name="departure_time_at_D")
    
    reqB = model.addVars(E, S, vtype=GRB.CONTINUOUS, lb=0, ub= GRB.INFINITY, name="reqB")       # battery required  
    
    
    u = model.addVars(D, E, S, vtype=GRB.CONTINUOUS, lb=1, ub=len(D), name="u")

    


    #obj = model.addVar(lb=0, name="obj")  # lb=0 ensures Z is non-negative
    #objective = gp.LinExpr()
    objective = gp.QuadExpr()
    
    ## Objective 2
    for j in E:
        for l in S:
            for x in C1:
                for y in D:
                    objective += chiCD[x, y, j, l] * alpha1
            for x in D:
                for y in D:
                    if x != y:
                        objective += chiDD[x, y, j, l] * alpha1
    for j in E:
        for l in S:
            for x in D:
                for y in C1:
                    term0_component = chiDC[x, y, j, l] * theta[y]  # This is symbolic during model building
                    objective -= term0_component * reqB[j, l] * alpha2
                    
    model.setObjective(objective, GRB.MAXIMIZE)
    
    
    
    
    
    # ## Objective 1
    # for j in E:
    #     for l in S:
    #         for x in D:
    #             for y in C1:
    #                 term0_component = chiDC[x, y, j, l] * theta[y]  # This is symbolic during model building
    #                 objective += term0_component * reqB[j, l]

    # objective += gp.quicksum(eta for x in D)
    
    # for j in E:
    #     for l in S:
    #         for x in C1:
    #             for y in D:
    #                 objective -= chiCD[x, y, j, l] * eta
    #         for x in D:
    #             for y in D:
    #                 if x != y:
    #                     objective -= chiDD[x, y, j, l] * eta
                                          
    # #model.addConstr(obj == objective, name="objective_constraint")
    # # Set the objective in the model (minimization)
    # #model.setObjective(obj, GRB.MINIMIZE)
    # model.setObjective(objective, GRB.MINIMIZE)
    
    
    
    
    #### Constraints ####
    # Additional Constraint: All self edges to be 0
    for j in E:
        for l in S:
            for x in D:
                model.addConstr(chiDD[x, x, j, l] == 0, name=f"chiDD_{j}_{l}_{x}_{x}")
                

    # Miller-Tucker-Zemlin (MTZ) Formulation: sub-tour elimination constraints
    n = len(D)
    for j in E:
        for l in S:
            for x in D:
                for y in D:
                    if x != y:
                        model.addConstr(u[x, j, l] - u[y, j, l] + n * chiDD[x, y, j, l] <= n - 1, f"mtz_{j}_{l}_{x}_{y}")


    # Every subtrip starts with a CS or D0 and end at a CS or D0  (Eq 5,6,7,8)
    # Create a linear expression
    # Initialize z_notIsEmptySubRoute variable: 1 when subtrip is not empty, 0 otherwise
    z_notIsEmptySubRoute = model.addVars(
        E, S, vtype=GRB.BINARY, name="z_notisemptysubroute")
    
    for j in E:
        for l in S:
            
            # Initialize emptySubRouteCheck as a linear expression
            emptySubRouteCheck = gp.LinExpr()

            for x in D:
                for y in D:
                    emptySubRouteCheck += chiDD[x, y, j, l]
                for y in C1:
                    emptySubRouteCheck += chiCD[y, x, j, l]
    
            # Big M Implementation for the constraints
            M = len(D) + 10000
            emptySubRouteCheck -= M * z_notIsEmptySubRoute[j, l]
            
            # First constraint: (sum of all edges to delivery points in a subroute) - M * z <= 0
            model.addConstr(emptySubRouteCheck <= 0, name=f"empty_sub_route_check1_{j}_{l}")
    
            # Second constraint: (sum of all edges to delivery points in a subroute) + M - M * z >= 1e-6
            emptySubRouteCheck += M  # Adding the M constant
            epsilon = 1e-6
            model.addConstr(emptySubRouteCheck >= epsilon, name=f"empty_sub_route_check2_{j}_{l}")
            
            # Constraints 1 and 2 for CP nodes at the start and end of subroutes
            constraint1 = gp.LinExpr()
            constraint2 = gp.LinExpr()
    
            for x in C1:
                for y in D:
                    constraint1 += chiCD[x, y, j, l]
    
            for x in D:
                for y in C1:
                    constraint2 += chiDC[x, y, j, l]
    
            # Adding constraints 5 and 6 to the model
            model.addConstr(constraint1 == z_notIsEmptySubRoute[j, l], name=f"c1_{j}_{l}")
            model.addConstr(constraint2 == z_notIsEmptySubRoute[j, l], name=f"c2_{j}_{l}")
    
    # Constraint 3: A delivery point should be visited at most once (Eq. 15)
    for x in D:
        constraint3 = gp.LinExpr()
        
        for j in E:
            for l in S:
                for y in D:
                    constraint3 += chiDD[y, x, j, l]
                
                for y in C1:
                    constraint3 += chiCD[y, x, j, l]
    
        # Add the constraint to the model: constraint3 <= 1
        model.addConstr(constraint3 <= 1, name=f"c3_{x}")
    
    
    # Constraint 4: An EV visiting a delivery point 'x' in trip 'l' should also come out from there (Eq. 14) # problem
    for j in E:
        for l in S:
            for x in D:
                constraint4A = gp.LinExpr()  # Initialize a linear expression
                constraint4B = gp.LinExpr()  # Initialize a linear expression
                
                for y in D:
                    constraint4A += chiDD[y, x, j, l]
                for y in D:
                    constraint4B += chiDD[x, y, j, l]
                    
                for y in C1:
                    constraint4A += chiCD[y, x, j, l]
                for y in C1:
                    constraint4B += chiDC[x, y, j, l]
    
                # Add the constraint to the model: constraint3 == 0
                model.addConstr(constraint4A - constraint4B == 0, name=f"c4_{j}_{l}_{x}")
                
    # Constraint 5, 6, 7: An EV starts from a depot and return at the end only once (Eq. 12)
    for j in E:
        constraint5 = gp.LinExpr()  # Initialize a linear expression
        for l in S:
            for x in D:
                constraint5 += chiCD[0, x, j, l]   
        model.addConstr(constraint5 <= 1, name=f"c5_{j}_{l}")
        
        constraint6 = gp.LinExpr()  # Initialize a linear expression
        for l in S:
            for x in D:
                constraint6 += chiDC[x, 0, j, l]
        model.addConstr(constraint6 <= 1, name=f"c5_{j}_{l}")

            # Add the constraint to the model: constraint3 == 0
        model.addConstr(constraint5 - constraint6 == 0, name=f"c5_{j}_{l}")
        
    
    # Constraint 8: An EV visiting a CP in trip 'l' should also come out from there at trip 'l+1'(Eq. 13)
    for j in E:
        for l in range(1, len(S)):
            for x in C:
                constraint8A = gp.LinExpr()  # Initialize a linear expression
                constraint8B = gp.LinExpr()  # Initialize a linear expression
                for y in D:   
                    constraint8A += chiDC[y, x, j, (l-1)]
                    constraint8B += chiCD[x, y, j, l]
        
                # Add the constraint to the model: constraint3 == 0
                model.addConstr(constraint8A == constraint8B, name=f"c8_{j}_{l}")
            
    # # Additional Constraint: If a subtrip is non-empty, its previous subtrip should also be non-empty
    # for j in E:
    #     for l in range(1, len(S)):  # Start from 1 to access the previous subtrip (l-1)
    #         constraint9 = z_notIsEmptySubRoute[j, l] - z_notIsEmptySubRoute[j, l-1]
            
    #         # Add the constraint to the model: z[k, l] - z[k, l-1] <= 0 # kept less than to avoid floating point error
    #         model.addConstr(constraint9 <= 0, name=f"c9_{j}_{l}")
    
    
    # Depot Constraint 1: A EV must start its first subtrip from the depot (fog 0)
    for j in E:
        depotConstraint1 = gp.LinExpr()
    
        # Sum over all tasks for the first subtrip (l=0)
        for y in D:
            depotConstraint1 += chiCD[0, y, j, 0]
    
        # Add the constraint to the model: depotConstraint1 == 1
        model.addConstr(depotConstraint1 == 1, name=f"depot_constraint1_{j}")
        
        
        
    # Energy Constraint
    for j in E:
        for l in S:
            constraint10 = gp.LinExpr()  # Initialize a linear expression
            for x in C1:
                for y in D:
                    constraint10 += chiCD[x, y, j, l] * psi_DC[(y, x)]
                
                for y in D:
                    constraint10 += chiDC[y, x, j, l] * psi_DC[(y, x)]
            for x in D:
                for y in D:
                    if x != y:
                        constraint10 += chiDD[x, y, j, l] * psi_DD[(x, y)]
                        
            # Add the constraint to the model: all energy ensumption for a subtrip <= beta_f
            model.addConstr(reqB[j, l] == constraint10, name=f"energy_constraint1_{j}_{l}")
            model.addConstr(constraint10 <= beta_f , name=f"energy_constraint_{j}_{l}")
            
            
    # Time Constraints
    for j in E:
        model.addConstr(TdepC[0, j, 0] >= 0, name=f"time_constraint1_{j}")
        
    for j in E:
        for l in S:
            for y in D:
                  constraint11 = gp.QuadExpr()  # Initialize a linear expression  
                  for x in D:
                      constraint11 += chiDD[x, y, j, l] * (TdepD[x] + gamma_DD[(x,y)])
                  for x in C1:
                      constraint11 += chiCD[x, y, j, l] * (TdepC[x, j, l] + gamma_DC[(y, x)])
                  model.addConstr(TarrD[y] >= constraint11, name=f"time_constraint2_{j}_{l}_{y}")
    
    
    for j in E:
        for l in S:
            for y in C:
                constraint12 = gp.QuadExpr()  # Initialize a linear expression  
                for x in D:
                    constraint12 += chiDC[x, y, j, l] * (TdepD[x] + gamma_DD[(x,y)])
                model.addConstr(TarrC[y, j, l] >= constraint12, name=f"time_constraint3_{j}_{l}_{y}")
    
    
    for y in D:
        # # If we have just the deadline of delivery
        # model.addConstr(TarrD[y] <= tau_end[y], name=f"time_constraint4C_{y}") 
        # model.addConstr(TdepD[y] >= TarrD[y], name=f"time_constraint4A_{y}")
        
        # If we have time slot (both upper and lower time limit) for delivery
        model.addConstr(TdepD[y] >= TarrD[y], name=f"time_constraint4A_{y}")
        model.addConstr(TdepD[y] >= tau_start[y], name=f"time_constraint4B_{y}")
        model.addConstr(TdepD[y] <= tau_end[y], name=f"time_constraint4C_{y}")
    
    for j in E:
        for l in range(1, len(S)):
            for y in C:
                constraint13 = gp.LinExpr()  # Initialize a linear expression
                chargingTime = reqB[j, l] * 60 / min(rateE[j], rateC[y]) # *60 for converting reqB from kWh to kW
                constraint13 = TarrC[y, j, (l-1)] + chargingTime
                model.addConstr(TdepC[y, j, l] >= constraint13, name=f"time_constraint5_{j}_{l}_{y}")
                    

    # Optimize the model
    model.optimize()
    
    end_time = time.time()
    elapsed_time = (end_time - start_time) * 1000
    
    # Print term1 after optimization if an optimal solution is found
    if model.status == GRB.OPTIMAL:
        for j in E:
            for l in S:
                for xp in D:
                    for yp in C1:
                        if chiCD[yp, xp, j, l].x == 1:
                            print(f"chiCD[{yp}, {xp}, {j}, {l}] = {chiCD[yp, xp, j, l].x:.2f}")
                        if chiDC[xp, yp, j, l].x == 1:
                            print(f"chiDC[{xp}, {yp}, {j}, {l}] = {chiDC[xp, yp, j, l].x:.2f}")
                    for yp in D:
                        if chiDD[xp, yp, j, l].x == 1:
                            print(f"chiDD[{xp}, {yp}, {j}, {l}] = {chiDD[xp, yp, j, l].x:.2f}")
                            
                print(f"reqB[{j},{l}] = {reqB[j,l].x:.2f}")
                
        for j in E:
            delivery_by_j = 0
            for l in S:
                for x1 in C1:
                    for y in D:
                        if chiCD[x1, y, j, l].x > 0.5:
                            total_deliveries_completed  += 1
                            delivery_by_j += 1
                for x1 in D:
                    for y in D:
                        if x1 != y and chiDD[x1, y, j, l].x > 0.5:
                            total_deliveries_completed  += 1
                            delivery_by_j += 1

            delivery_distribution[j] = delivery_by_j
            
        for j in E:
            for l in S:
                for x in D:
                    for y in C1:
                        if chiDC[x, y, j, l].x > 0.5:
                            total_cost += theta[y] * reqB[j, l].x  # This is symbolic during model building
                            
    
    # if model.status == GRB.OPTIMAL:
    #     for j in E:
    #         for l in S:
    #             for x in D:
    #                 for y in C1:
    #                     # Evaluate term1 using the optimized values of chiCD, chiDD, and chiDC
    #                     term1_value = (
    #                         beta_f -
    #                         sum(chiCD[xp, yp, j, l].x * psi_DC[(xp, yp)] for xp in C1 for yp in D) -
    #                         sum(chiDD[xp, yp, j, l].x * psi_DD[(xp, yp)] for xp in D for yp in D) -
    #                         sum(chiDC[xp, yp, j, l].x * psi_DC[(xp, yp)] for xp in D for yp in C1)
    #                     )
                        
    #                     print(f"Evaluated term1 for j={j}, l={l}, x={x}, y={y}: {term1_value:.2f}")
    # else:
    #     print("No optimal solution found.")


print("****Optimal****") 
print("Total delivery: ", nD, " Total CP: ", nC+1, " Total EV: ", nE)
print(f"Elapsed time: {elapsed_time:.2f} ms")
print(f"Total cost: {total_cost}")
print(f"Total successful delivery: {total_deliveries_completed}")
print("::Delivery Distribution::")
for j in E:
    print(f"delivery_distribution[{j}]:", delivery_distribution.get(j, 0))
    







    
# #    # Use FeasRelax to relax constraints and find a feasible solution
# #    #feas_relax_result = model.feasRelax(relaxobjtype=0, minrelax=False)

#     # Optimize the model
#     model.optimize()

#     # Infeasibility Diagnosis: Use IIS (Irreducible Inconsistent Subsystem)
#     #model.computeIIS()
#     #model.write("infeasible_model.ilp") # This will write the IIS to a file for review

#     #model.write("model.lp")     

    
# #    # Print the results
# #    if model.status == GRB.OPTIMAL:
# #        print(f"Optimal tour cost: {model.objVal}")
# #        solution = model.getAttr('chi', chi)
# #        
# #        # Extract the tour from the solution
# #        tour = []
# #        current_city = 0
# #        while len(tour) < nD:
# #            tour.append(current_city)
# #            for j in range(nD):
# #                if solution[current_city, j] > 0.5:
# #                    current_city = j
# #                    break
# #        tour.append(tour[0])  # To complete the tour by returning to the starting city
# #        
# #        print(f"Optimal tour: {tour}")
# #    
# #    if model.status == GRB.INFEASIBLE:
# #        vars = model.getVars()
# #        ubpen = [1.0]*model.numVars
# #        model.feasRelax(1, False, vars, None, ubpen, None, None)
# #        model.optimize()
