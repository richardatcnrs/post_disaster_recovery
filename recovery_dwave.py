import dimod
import dwave.embedding
from dwave.system import LeapHybridCQMSampler 

import sys
from itertools import groupby


def read_input():
    input_file = open(sys.argv[1], 'r')

    n = int(input_file.readline().strip())
    weight_matrix = [[None for i in range(n)] for y in range(n)]
    edge_list = []

    for i in range(n):
        weights = input_file.readline().split()
        #print(weights)
        for j in range(n):
            #print(i,j)
            temp = int(weights[j].strip()) 
            weight_matrix[i][j] = temp
            if not(temp == 0):
                edge_list.append((i,j))

    broken_nodes = eval(input_file.readline().strip())
    node_capacity = eval(input_file.readline().strip())
    input_file.close()
    k = len(broken_nodes)
    return n,weight_matrix,broken_nodes,edge_list,node_capacity

def print_constraint(model):
    for const in model.constraints:
        
        print(model.constraints[const].to_polystring())

def build_model(n, weight_matrix, broken_nodes, edge_list, node_capacity):
    model = dimod.ConstrainedQuadraticModel()
    k = len(broken_nodes)
     
    # setting up the variables to encode the flow values
    # lower and upper bounds set at initialization
    fv_vars = {}
    for v in range(n):
        for t in range(k+1):
            var_label = 'fv_' + str(v) + '_' + str(t)
            fv_vars[v,t] = var_label
            #print(node_capacity[v])
            #node_flow_vars[v,t] = 
            #append(var_label)
            model.add_variable('INTEGER',var_label,lower_bound=0,upper_bound=node_capacity[v])
    #print(var_list)
    #model.add_variables('INTEGER', node_flow_vars)
    
    #print('ffffffffffff')
          
    fe_vars = {}
    for u in range(n):
        for v in range(n):
            for t in range(k+1):
                var_label = 'fe_' + str(u) + '_' + str(v) + '_' + str(t)
                fe_vars[u,v,t] = var_label
                if (u,v) in edge_list:
                    model.add_variable('INTEGER',var_label,lower_bound=0,upper_bound=weight_matrix[u][v])
                else:
                    model.add_variable('INTEGER',var_label,lower_bound=0,upper_bound=0)
    
    # setting up the variables for the recovery sequence

    recovery_vars = {}
    for u in broken_nodes:
        for t in range(k+1):
            var_label = 'x_' + str(u) + '_' + str(t)
            recovery_vars[u,t] = var_label
            model.add_variable('BINARY', var_label)

    #print(model.variables)

  

    # constraints for the recovery matrix
    # sum x_u,t = t for each t
    # there are t functional vertices at time t
    for t in range(1, k+1):
        terms = []
        for u in broken_nodes:
            var_label = recovery_vars[u,t]
            terms.append([var_label,1])
        model.add_constraint_from_iterable(terms, '==', rhs = t, label = 'Number of functional vertices at time ' + str(t))
    
    # x_u,t <= x_u,t+1 for each broken u and 2 <= t <= k
    # fixed vertices must remain functional
    # convert to x_u,t+1 - x_u,t >= 0 for the API

    for t in range(k):
        for u in broken_nodes:
            var_1 = recovery_vars[u,t+1]
            var_2 = recovery_vars[u,t]
            terms = [[var_1,1],[var_2,-1]]
            model.add_constraint_from_iterable(terms,'>=',rhs = 0, label = 'Vertex ' + str(u) + ' remains functional at time ' + str(t))



    
    # constraints for the correct flow values
    
    # f_v,t <= W(v)*x_v,t for all v and t
    # flow can only reach vertex v at step t if v is functional at step t
    # convert to W(v)*x_v,t - f_v,t >= 0
    for t in range(1,k):
        for v in broken_nodes: # source is always at full capacity
            var_1 = recovery_vars[v,t]

            var_2 = fv_vars[v,t]
            terms = [[var_1,node_capacity[v]],[var_2,-1]]
            model.add_constraint_from_iterable(terms,'>=',rhs = 0, label = 'FLow reach vertex ' + str(v) + ' at time ' + str(t) + ' if fixed') 

    
    # f_v,t = sum_u->v f_(u,v),t for all v and t
    # flow value reaching vertex v is equal to the sum of flow values of 
    # all edges reaching v
    # convert to f_v,t - sum_u->v f_(u,v),t = 0

    for t in range(1,k):
        for v in range(1,n): # don't need to consider folw reaching the source
            terms = [[fv_vars[v,t],1]]
            for u in range(n):
                if (u,v) in edge_list:
                    terms.append([fe_vars[u,v,t],-1])

            model.add_constraint_from_iterable(terms,'==',rhs = 0, label = 'Flow reaching ' + str(v) + ' at time ' + str(t))


    # f_v,t >= sum_v->u f_(v,u),t for all v and t
    # flow leaving vertex v cannot exceed flow reaching v
    # convert f_v,t - sum_v->u f_(v,u),t >= 0 
    for t in range(1,k):
        for v in range(n-1): # don't need to consider flow leaving the sink
            terms = [[fv_vars[v,t],1]]
            for u in range(n):
                if (v,u) in edge_list:
                    terms.append([fe_vars[v,u,t],-1])
            model.add_constraint_from_iterable(terms,'>=',rhs = 0, label = 'Flow leaving ' + str(v) + ' at time ' + str(t))
    

    obj = []
    
    for t in range(1,k+1):
       obj.append([fv_vars[n-1,t],-1])
    model.set_objective(obj)
     

    # setting constant values for vars

    # soure always have max flow value
    for t in range(1,k+1):
        model.fix_variable(fv_vars[0,t], node_capacity[0])

    # source and sink always functional
    #for t in range(k+1):
    #    model.fix_variable()
    
    # broken vertices start as broken
    for u in broken_nodes:
        model.fix_variable(recovery_vars[u,0],0)
    
    var_maps = {}
    var_maps['recovery_vars'] = recovery_vars
    var_maps['fv_vars'] = fv_vars
    var_maps['fe_vars'] = fe_vars


    #print('kkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkkk')
    return model,var_maps

def run_cqm_and_collect_solutions(model, sampler,time_limit):
    #time_limit = sampler.min_time_limit(model)
    print('solve time =', 35)
    #print('solve time =', time_limit)
    sampleset = sampler.sample_cqm(model,time_limit=35)
    #sampleset = sampler.sample_cqm(model)
    return sampleset

def process_solutions(sampleset,n,broken_nodes,var_maps,edge_list):
    k = len(broken_nodes)
    feasible_solutions = []
    recovery_vars = var_maps['recovery_vars']
    fv_vars = var_maps['fv_vars']
    fe_vars = var_maps['fe_vars']
    sorted_sampleset = sampleset.truncate(100, sorted_by='energy')

    for solution in sampleset:
        #print(solution)
        if(model.check_feasible(solution, atol = 0)):
        #if check_feasibility(solution, n, k) == True:
            feasible_solutions.append(solution)
    #print(sampleset.variables) 
    for record in sorted_sampleset.record:
        if record.is_feasible:
            sample = record.sample
            energy = record.energy
            recovery_seq = []
            for t in range(1,k+1):
                for u in broken_nodes:
                    # fetch variable index in the sample
                    var_label = recovery_vars[u,t]
                    index = sampleset.variables.index(var_label)
                    if sample[index] == 1 and not(u in recovery_seq):
                        recovery_seq.append(u)
            print(recovery_seq)
            print('Energy = ', record.energy)
        
            # print flow value for sink at each iteration

            for t in range(1,k+1):
                var_label = fv_vars[n-1,t]
                index = sampleset.variables.index(var_label) 
                print('fv_' + str(n-1) + '_' + str(t) + ' = ', str(sample[index]))

            # print all vertex flow values
            #for t in range(1,k+1):
            #    for u in range(1,n):
            #        var_label = fv_vars[u,t]
            #        index = sampleset.variables.index(var_label)
            #        print('fv_' + str(u) + '_' + str(t) + ' = ' + str(sample[index]))
            #print('ssssssssssssssssssssssssssss')

            #for t in range(1,k+1):
            #    for (u,v) in edge_list:
            #        var_label = fe_vars[u,v,t]
            #        index = sampleset.variables.index(var_label)
            #        print('fe_' + str(u) + '_' + str(v) + '_' + str(t) + ' = ' + str(sample[index]))
def check_feasibility(solution, n, m, T, weight_matrix):
   
    #print(solution)
    

    # constraint 1 - all vertices are visited
    vertices_visited = []
    for v in range(1, n):
        for i in range(m):
            for t in range(1, T+1):
                var_label = 'x_' + str(i) + '_' + str(v) + '_' + str(t)
                # check if v is visited
                if solution[var_label] == 1:
                    if v not in vertices_visited:
                        vertices_visited.append(v)
                    else:
                        print(v, ' is visited more than once')
                        print(solution)
                        sys.exit()

# set time limit
time_limit = 40

token_file = open('/home/richard/Desktop/data/dwave_token','r')
token = token_file.readline()
#print(token)
token_file.close()


n, weight_matrix, broken_nodes,edge_list, node_capacity = read_input()
sampler = LeapHybridCQMSampler()
model,var_maps = build_model(n, weight_matrix, broken_nodes, edge_list, node_capacity)
sampleset = run_cqm_and_collect_solutions(model, sampler,time_limit)

process_solutions(sampleset,n,broken_nodes,var_maps,edge_list)

#    print(best_makespan)
#    for route in routes:
#        print(route)

#print('Optimal cost =', min_cost)
#print('Optimal path =', min_path)





