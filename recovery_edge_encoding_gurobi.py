"""Recovery optimization with Gurobi: minimum cumulative lack of resilience.

ONE-TIME MAC SETUP (in Terminal):
    python3 -m pip install gurobipy
    Activate your Academic Named-User license using grbgetkey
    (from the official Gurobi installation). The license must already be
    configured; this script does not need the license key in its source.

RUN (from the directory containing this script and your input file):
    python3 recovery_gurobi_post.py 43_equipment_case2.txt

The input file uses the same format as recovery_dwave.py.
Resilience loss is measured AFTER each recovery action:
    sum(60 - sink_flow[t] for t = 1, ..., number_of_broken_nodes)

Writes <input_filename>_results.csv next to the input file.
"""

import ast
import csv
from pathlib import Path
import sys
import gurobipy as gp
from gurobipy import GRB

MAX_CAPACITY = 60  # Fully recovered system flow (per recovery step)

# some infos are hard-coded for the given example
output_edges = [(19,25),(20,26),(21,27),(22,28),(23,29),(24,30)]
edge_to_equipment = {
    (1, 4): 1,
    (2, 5): 2,
    (3, 6): 3,
    (4, 7): 4,
    (4, 8): 5,
    (5, 7): 6,
    (5, 8): 7,
    (6, 7): 8,
    (6, 8): 9,
    (7, 8): 10,
    (8, 7): 10,
    (7, 9): 11,
    (8, 10): 12,
    (9, 11): 13,
    (9, 12): 14,
    (10, 11): 15,
    (10, 12): 16,
    (11, 13): 17,
    (12, 14): 18,
    (13, 15): 19,
    (13, 16): 20,
    (14, 15): 21,
    (14, 16): 22,
    (15, 16): 23,
    (16, 15): 23,
    (15, 17): 24,
    (16, 18): 25,
    (17, 19): 26,
    (17, 20): 27,
    (17, 21): 28,
    (17, 22): 29,
    (17, 23): 30,
    (17, 24): 31,
    (18, 19): 32,
    (18, 20): 33,
    (18, 21): 34,
    (18, 22): 35,
    (18, 23): 36,
    (18, 24): 37,
    (19, 25): 38,
    (20, 26): 39,
    (21, 27): 40,
    (22, 28): 41,
    (23, 29): 42,
    (24, 30): 43
}

def read_input():
    input_file = open(sys.argv[1], 'r')
    n = int(input_file.readline().strip())
    
    edge_list = []
    edge_weights = {}
    source = 0
    sink = n-1
    
    for u in range(n):
        line = input_file.readline().split()
        for v in range(n):
            weight = int(line[v].strip())

            # if there is an edge from u to v
            if not(weight == 0):
                edge_list.append((u,v))
                edge_weights[(u,v)] = weight

    broken_edges = eval(input_file.readline().strip())


    input_file.close()

    return n, source, sink, edge_list, edge_weights, broken_edges
#        weight_matrix = []
#        edge_list = []
#        for i in range(n):
#            weights = list(map(int, input_file.readline().split()))
#            weight_matrix.append(weights)
#            for j, weight in enumerate(weights):
#                if weight != 0:
#                    edge_list.append((i, j))
#        broken_nodes = ast.literal_eval(input_file.readline().strip())
#        node_capacity = ast.literal_eval(input_file.readline().strip())
#    return n, weight_matrix, broken_nodes, edge_list, node_capacity


def build_model(n, source, sink, edge_list, edge_weights, broken_edges):
    model = gp.Model('recovery')


    # hard-code 2 special cases (E10 and E23), only keep 1 of the 
    # bidirectional edge since both directions are repaired at the same time
    if (7,8) in broken_edges and (8,7) in broken_edges:
        broken_edges.remove((8,7))
    if (15,16) in broken_edges and (16,15) in broken_edges:
        broken_edges.remove((16,15))
    
    k = len(broken_edges)
    # define recovery variables
    recovery_vars = {(u,v,t): model.addVar(vtype=GRB.BINARY, name=f'x_{u}_{v}_{t}')
            for (u,v) in broken_edges for t in range(1,k+1)
            }
    
    # define flow variables
    y_vars = {(u,v,t): model.addVar(vtype=GRB.INTEGER,lb=0,
                                    ub=max(0, edge_weights[u,v]),
                                    name=f'y_{u}_{v}_{t}')
            for (u,v) in edge_list for t in range(1,k+1)
            }


    #print(y_vars)

    # Exactly t broken vertices have been recovered at step t.
    for t in range(1, k + 1):
        model.addConstr(gp.quicksum(recovery_vars[u,v,t] for (u,v) in broken_edges) == t,
                        name=f'number_recovered_{t}')

    # Recovered vertices remain functional.
    for t in range(1,k):
        for (u,v) in broken_edges:
            model.addConstr(recovery_vars[u,v,t+1] >= recovery_vars[u,v,t],
                            name=f'remains_functional_{u}_{v}_{t}')

    # Flow-conservation constraints
    # 1. y_e_t <= W(e)x_e_t i.e., flow can only pass through the edge if it is fixed for
    # broken edges

    for (u,v) in broken_edges:
        for t in range(1,k+1):
            model.addConstr(y_vars[u,v,t] <= edge_weights[u,v]*recovery_vars[u,v,t],
                    name=f'flow_value_bound_{u}_{v}_{t}')
    
    # two special cases for (8,7) and (16,15)
    if (7,8) in broken_edges:
        for t in range(1,k+1):
            model.addConstr(y_vars[8,7,t] <= edge_weights[8,7]*recovery_vars[7,8,t],
                    name=f'flow_value_bound_8_7_{t}')
    if (15,16) in broken_edges:
        for t in range(1,k+1):
            model.addConstr(y_vars[16,15,t] <= edge_weights[16,15]*recovery_vars[15,16,t],
                    name=f'flow_value_bound_16_15_{t}')
                   
    # 2. y_e_t <= W(e) for functional edges, do not need to encode since the upper-bound is 
    # set explicitly when defining the flow variables

    # 3. sum_u_v_t y_u_v_t = sum_v_w_t y_v_w_t for nodes other than source and sink 
    # i.e., incoming flow is equal to outgoing flow
    model.update()
    for v in range(4,24):
        for t in range(1,k+1):
            #print(v,t)
            LHS=gp.quicksum(y_vars[u,v,t] for u in range(n) if (u,v) in edge_list)
            #print(str(LHS))
            RHS=gp.quicksum(y_vars[v,w,t] for w in range(n) if (v,w) in edge_list)
            #print(str(RHS))
            model.addConstr(LHS == RHS, name=f'flow_conserv_{v}_{t}')

    # Set objective, note that the super sink is not used in order to save some
    # variables, we directly sum over the edges in the last layer of the network
    # last layer of edges (output_edges) is hard-coded at the top of the script
    objective = 0
    for t in range(1,k+1):
        objective += (MAX_CAPACITY - sum(y_vars[u,v,t] for (u,v) in output_edges))
    
    model.setObjective(objective,GRB.MINIMIZE)

    return model,recovery_vars,y_vars

def process_solution(model,n,broken_edges,recovery_vars,y_vars):

    if model.SolCount == 0:
        print('No feasible solution found. Status:', model.Status)
        model.computeIIS()
        model.write('iis.ilp')
        return
    if model.Status == 2:
        print('Solution is optimal')
    k = len(broken_edges)
    
    #for t in range(1,k+1):
    #    for (u,v) in broken_edges:
    #        print(recovery_vars[u,v,t].X, end=' ')
    #    print()



    recovery_seq = []
    # frist recovery is calculated separately (since we do not have
    # variables for step 0)
    for (u,v) in broken_edges:
        if recovery_vars[u,v,1].X > 0.5:
            recovered_edge_index = edge_to_equipment[(u,v)]
            recovery_seq.append(recovered_edge_index)
    # Exactly one edge is recovered at teach step
    for t in range(2,k+1):
        recovered = [(u,v) for (u,v) in broken_edges
                    if recovery_vars[u,v,t].X > 0.5
                    and recovery_vars[u,v,t-1].X < 0.5
                ]
        if len(recovered) != 1:
            raise ValueError(f'Expected one recovered edge at step {t}; found {recovered}')
        recovered_edge_index = edge_to_equipment[recovered[0]]
        recovery_seq.append(recovered_edge_index)
    for t in range(1,k+1):
        flow = 0
        for (u,v) in output_edges:
            flow += y_vars[u,v,t].X
    
        print(f'flow at iteration {t} = ',flow)

    print('objective = ', model.ObjVal)

    print('recovery_seq = ',recovery_seq)

    #for (u,v) in edge_list:
    #    if (u,v) in edge_to_equipment:
    #        e = edge_to_equipment[(u,v)]
    #        print(f'E{e} = ', y_vars[u,v,1].X)
n, source, sink, edge_list, edge_weights, broken_edges = read_input()
model,recovery_vars,y_vars = build_model(n, source, sink, edge_list, edge_weights, broken_edges)
model.Params.TimeLimit = 30
model.Params.OutputFlag = False
model.optimize()

process_solution(model,n,broken_edges,recovery_vars,y_vars)
#print(n)
#print(source)
#print(sink)
#print(edge_list)
#print(edge_weights)
#print(broken_edges)
#if __name__ == '__main__':
#    if len(sys.argv) < 2:
#        sys.exit('Usage: python3 recovery_gurobi_post.py input_file.txt')
#    n, weight_matrix, broken_nodes, edge_list, node_capacity = read_input()
#    model, var_maps = build_model(n, weight_matrix, broken_nodes, edge_list, node_capacity)
#    model.Params.TimeLimit = 40
#    try:
#        model.optimize()
#    except gp.GurobiError as exc:
#        sys.exit(f'Gurobi error: {exc}. Check that your academic license is activated.')
#    process_solution(model, n, broken_nodes, var_maps)
