from amplify import VariableGenerator
from amplify import equal_to, less_equal, greater_equal,Model
from amplify import FixstarsClient, solve
from amplify import sum as amplifysum
from datetime import timedelta
import networkx as nx
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
    return n,weight_matrix,broken_nodes,edge_list,node_capacity

def decode_solution(values,recovery_var,broken_nodes,k):
    
    # output the recovery matrix
    print('Recovery matrix =')
    for v in broken_nodes:
        for t in range(1,k+1):
            print(values[recovery_var[v,t]], end=' ')
        print()


    recovery_sequence = []
    for u in broken_nodes:
        if values[recovery_var[u,1]] == 1:
            recovery_sequence.append(u)
    for t in range(2,k+1):
        for u in broken_nodes:
            if values[recovery_var[u,t-1]] == 0 and values[recovery_var[u,t]] == 1:
                recovery_sequence.append(u)

    print('Recovery sequence =',recovery_sequence)
    print('Flow value after each iteration')
    for t in range(1,k+1):
        print(values[node_flow_var[n-1,t]])

# solver time in milli seconds

time = 5000
if '-t' in sys.argv:
    # convert input time to seconds
    time = int(sys.argv[sys.argv.index('-t')+1].strip())*1000

# replace token with your access token from Fixstars Amplify
token_file = open('/home/richard/Desktop/data/fs_token','r')
token = token_file.readline().strip()
#print(token)
token_file.close()


n, weight_matrix, broken_nodes,edge_list,node_capacity = read_input()
#print(n) 
#print(weight_matrix)
#print(broken_nodes)




# setting up the variables for the recovery matrix
# note that I'm assuming the broken vertices are labelled 
# from 1 to k
k = len(broken_nodes)
gen = VariableGenerator()

recovery_var = gen.array('Binary', shape=(n,k+1))
#print(recovery_var)

# setting up variables to encode the flow values
# note that not all variables defined are used
# upper bounds for the flow values can be set here or as constraints
node_flow_var = gen.array('Integer',shape=(n,k+1),bounds=(0,None))
for v in range(n):
    for t in range(k+1):    
        node_flow_var[v,t].lower_bound = 0
        node_flow_var[v,t].upper_bound = node_capacity[v]


edge_flow_var = gen.array('Integer',shape=(n,n,k+1),bounds=(0,None))
for u in range(n):
    for v in range(n):
        for t in range(k+1):
            edge_flow_var[u,v,t].lower_bound = 0
            if (u,v) in edge_list:
                edge_flow_var[u,v,t].upper_bound = weight_matrix[u][v]
            else:

                edge_flow_var[u,v,t].upper_bound = 0






# setting up the constraints
# source always have the max flow value
node_flow_var[0,:] = node_capacity[0]

# source and sink are always functional
recovery_var[0,:] = 1
recovery_var[n-1:] = 1

# broken nodes start as broken
for u in broken_nodes:
    recovery_var[u,0] = 0
#print(recovery_var)

# constraints for the recovery matrix
# sum 1 <= i <= k x_i,t = t for each 1 <= t <= k there are 
# t functional vertices at time t.
constraint_set_1 = []
for t in range(1,k+1):
    constraint = equal_to(sum(recovery_var[u,t] for u in broken_nodes),t,label='there are ' + str(t) + ' fixed nodes at step ' + str(t))
    constraint_set_1.append(constraint) 

# x_i,t <= x_i,t+1 for each 1 <= i <= k and 2 <= t <= k fixed vertex 
# must remain functional
# need to convert this to x_i,t+1 - x_i,t >= 0 for the API

constraint_set_2 = []
for t in range(k):
    for u in broken_nodes:
        constraint = greater_equal(recovery_var[u,t+1]-recovery_var[u,t],0,label='if vertex ' + str(u) + ' is fixed at step ' + str(t) + ' it reamins fixed')
        constraint_set_2.append(constraint)




# constraints for the correct flow values
# f_v,t <= W(v)*x_v,t for all v and 1 <= t <= k
# flow can only reach vertex v at step t if v is functional at step t
constraint_set_3 = []
for t in range(1,k):
    for u in range(n):
        constraint = greater_equal(node_capacity[u]*recovery_var[u,t]-node_flow_var[u,t],0,label='Node flow constraint for vertex ' + str(u) + ' at step ' + str(t))
        constraint_set_3.append(constraint)

# f_v,t = sum_u->v f_(u,v),i for all v and 1 <= t <= k
# flow value reaching vertex v is equal to the sum of flow values of all 
# edges reaching v
constraint_set_4 = []
for t in range(1,k):
    for v in range(1,n): # we don't need to consider flow reaching the source
        expression = node_flow_var[v,t]
        for u in range(n): 
            if (u,v) in edge_list:
                expression = expression - edge_flow_var[u,v,t]
        constraint = equal_to(expression,0,label='flow reaching vertex' + str(v) + ' at step ' + str(t))
        constraint_set_4.append(constraint) 
                
# f_v,t >= sum_v->u f_(v,u),i for all v and 1 <= t <= k
# flow leaving vertex v cannot exceed flow reaching v
constraint_set_5 = []
for t in range(1,k):
    for v in range(n-1): # we don't need to consider flow leaving the sink
        expression = node_flow_var[v,t]
        #print(expression)
        for u in range(n):
            if (v,u) in edge_list:

                #print(v, 'ccccccccccc',(v,u))
                expression = expression - edge_flow_var[v,u,t]
                #print(expression)
        constraint = greater_equal(expression,0,label='flow leaving vertex ' + str(v) + ' at step ' + str(t))
        constraint_set_5.append(constraint)





# set up the objective
# min sum -f_sink,t   
# minimizing the flow reaching the sink at each iteration

objective = sum(node_flow_var[n-1,t] for t in range(1,k+1))

model = Model(-1*objective)

# adjusting penalty weights before adding to model, current values are arbitrary
a = 500
b = 700
c = 200
d = 200
e = 200

for constraint in constraint_set_1:
    constraint.weight = a

for constraint in constraint_set_2:
    constraint.weight = b

for constraint in constraint_set_3:
    constraint.weight = c

for constraint in constraint_set_4:
    constraint.weight = d

for constraint in constraint_set_5:
    constraint.weight = e


model += amplifysum(constraint for constraint in constraint_set_1)
model += amplifysum(constraint for constraint in constraint_set_2)
model += amplifysum(constraint for constraint in constraint_set_3)
#print(model)
model += amplifysum(constraint for constraint in constraint_set_4)
model += amplifysum(constraint for constraint in constraint_set_5)
#print(model)

client = FixstarsClient()
client.token = token.strip()
#client.parameters.timeout = timedelta(milliseconds=time)
client.parameters.timeout = timedelta(milliseconds=10000)


result = solve(model,client,filter_solution=False)
#print(result)
#print(result.best.objective)
#print(result.best.feasible)
#values = (result.best.values)
#print(len(result.solutions))

for solution in result.solutions:
    # output solution feasibility and unsatisfied constraints (if any)
    print('all constraints satisfied =', solution.feasible)
    print('unsatisfied constraints =')
    infeasible_constraints = []
    for constraint in constraint_set_1:
        if not(constraint.is_satisfied(solution.values)):
            infeasible_constraints.append(constraint)
    print('Set 1 ------------------------')
    #print(infeasible_constraints)
    for constraint in infeasible_constraints:
        print(constraint)


    infeasible_constraints = []
    for constraint in constraint_set_2:
        if not(constraint.is_satisfied(solution.values)):
            infeasible_constraints.append(constraint)
    print('Set 2 ------------------------')
    for constraint in infeasible_constraints:
        print(constraint)


    infeasible_constraints = []
    for constraint in constraint_set_3:
        if not(constraint.is_satisfied(solution.values)):
            infeasible_constraints.append(constraint)
    print('Set 3 ------------------------')
    for constraint in infeasible_constraints:
        print(constraint)


    infeasible_constraints = []
    for constraint in constraint_set_4:
        if not(constraint.is_satisfied(solution.values)):
            infeasible_constraints.append(constraint)
    print('Set 4 ------------------------')
    for constraint in infeasible_constraints:
        print(constraint)
   

    infeasible_constraints = []
    for constraint in constraint_set_5:
        if not(constraint.is_satisfied(solution.values)):
            infeasible_constraints.append(constraint)
    print('Set 5 ------------------------')
    for constraint in infeasible_constraints:
        print(constraint)
   
    # output the recovery matrix if solution is feasible
    if solution.feasible:
            decode_solution(solution.values,recovery_var,broken_nodes,k)
    print('Solution objective value =', solution.objective)


