"""Recovery optimization with Gurobi MIQP and explicit QUBO penalties.

SETUP:
    python3 -m pip install gurobipy
    Configure your Gurobi license before running this script.

RUN:
    python3 recovery_edge_encoding_qubo_miqp.py input_file -t 30

The input format and post-repair resilience-loss objective match the original.
All decision variables are binary. Flow and slack integers use bounded binary
encodings. Repair counts, persistence, damaged-edge capacity and conservation
are embedded directly in the quadratic objective.

Output includes each penalty, the QUBO energy, the original resilience loss,
the throughput after each repair and the equipment recovery sequence.
"""

import ast
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
    with open(sys.argv[1], 'r') as input_file:
        n = int(input_file.readline().strip())
        if n < 31:
            raise ValueError('The hard-coded output edges require at least 31 nodes.')

        edge_list = []
        edge_weights = {}
        source = 0
        sink = n - 1

        for u in range(n):
            line = input_file.readline().split()
            if len(line) != n:
                raise ValueError(f'Expected {n} capacities in row {u}.')
            for v in range(n):
                weight = int(line[v])
                if weight < 0:
                    raise ValueError('Edge capacities must be nonnegative integers.')
                if weight != 0:
                    edge_list.append((u, v))
                    edge_weights[u, v] = weight

        damaged = ast.literal_eval(input_file.readline().strip())

    if not isinstance(damaged, (list, tuple, set)):
        raise ValueError('Expected a collection of broken edges.')

    broken_edges = []
    for edge in damaged:
        if (not isinstance(edge, (list, tuple)) or len(edge) != 2
                or any(type(v) is not int for v in edge)):
            raise ValueError(f'Invalid broken edge: {edge!r}')
        edge = tuple(edge)
        if edge not in edge_weights or edge not in edge_to_equipment:
            raise ValueError(f'Unknown or absent equipment edge: {edge}')
        if edge not in broken_edges:
            broken_edges.append(edge)

    for edge in output_edges:
        if edge not in edge_weights:
            raise ValueError(f'Missing output edge: {edge}')

    return n, source, sink, edge_list, edge_weights, broken_edges


def binary_integer(model, upper_bound, name):
    # The weights sum to upper_bound and represent every integer from 0 to it.
    # Example: upper_bound = 10 gives weights [1, 2, 4, 3].
    if upper_bound == 0:
        return gp.LinExpr(0)

    bit_count = upper_bound.bit_length()
    weights = [2 ** j for j in range(bit_count - 1)]
    weights.append(upper_bound - sum(weights))

    bits = [
        model.addVar(vtype=GRB.BINARY, name=f'{name}_b{j}')
        for j in range(bit_count)
    ]
    return gp.quicksum(weight * bit for weight, bit in zip(weights, bits))


def build_model(n, source, sink, edge_list, edge_weights, broken_edges):
    model = gp.Model('recovery_qubo_miqp')

    # E10 and E23 each have two directed edges but only one repair action.
    # Canonicalize either listed direction to the forward edge.
    for forward, reverse in [((7, 8), (8, 7)), ((15, 16), (16, 15))]:
        if forward in broken_edges or reverse in broken_edges:
            if forward not in edge_weights or reverse not in edge_weights:
                raise ValueError(f'Shared equipment requires both {forward} and {reverse}.')
            broken_edges[:] = list(dict.fromkeys(
                forward if edge == reverse else edge
                for edge in broken_edges
            ))

    k = len(broken_edges)

    recovery_vars = {
        (u, v, t): model.addVar(vtype=GRB.BINARY, name=f'x_{u}_{v}_{t}')
        for u, v in broken_edges for t in range(1, k + 1)
    }

    # Integer flows are linear expressions in binary variables.
    # Their encodings enforce 0 <= y[u,v,t] <= edge_weights[u,v].
    y_vars = {
        (u, v, t): binary_integer(model, edge_weights[u, v], f'y_{u}_{v}_{t}')
        for u, v in edge_list for t in range(1, k + 1)
    }

    count_penalty = gp.QuadExpr()
    persistence_penalty = gp.QuadExpr()
    capacity_penalty = gp.QuadExpr()
    conservation_penalty = gp.QuadExpr()

    # Exactly t pieces of equipment have been repaired at step t.
    for t in range(1, k + 1):
        residual = gp.quicksum(
            recovery_vars[u, v, t] for u, v in broken_edges
        ) - t
        count_penalty += residual * residual

    # Penalize a repaired edge becoming broken again: x_t * (1 - x_{t+1}).
    for u, v in broken_edges:
        for t in range(1, k):
            persistence_penalty += (
                recovery_vars[u, v, t] * (1 - recovery_vars[u, v, t + 1])
            )

    repair_state = {edge: edge for edge in broken_edges}
    if (7, 8) in repair_state:
        repair_state[8, 7] = (7, 8)
    if (15, 16) in repair_state:
        repair_state[16, 15] = (15, 16)

    # y <= W*x becomes y + slack = W*x, with 0 <= slack <= W.
    for (u, v), (repair_u, repair_v) in repair_state.items():
        capacity = edge_weights[u, v]
        for t in range(1, k + 1):
            slack = binary_integer(model, capacity, f'slack_{u}_{v}_{t}')
            residual = (
                y_vars[u, v, t] + slack
                - capacity * recovery_vars[repair_u, repair_v, t]
            )
            capacity_penalty += residual * residual

    # Match the attached script: flow conservation at nodes 4 through 24.
    for v in range(4, 25):
        for t in range(1, k + 1):
            LHS = gp.quicksum(
                y_vars[u, v, t] for u in range(n) if (u, v) in edge_weights
            )
            RHS = gp.quicksum(
                y_vars[v, w, t] for w in range(n) if (v, w) in edge_weights
            )
            residual = LHS - RHS
            conservation_penalty += residual * residual

    # Resilience loss is measured AFTER each recovery action.
    objective = gp.LinExpr(0)
    for t in range(1, k + 1):
        objective += MAX_CAPACITY - gp.quicksum(
            y_vars[u, v, t] for u, v in output_edges
        )

    # This exceeds the full range of the unpenalized objective.
    # Every nonzero penalty is an integer >= 1, so a global minimizer
    # satisfies all constraints. A time-limited incumbent may violate them.
    penalty_weight = k * sum(edge_weights[edge] for edge in output_edges) + 1
    total_penalty = (
        count_penalty + persistence_penalty
        + capacity_penalty + conservation_penalty
    )
    model.setObjective(objective + penalty_weight * total_penalty, GRB.MINIMIZE)

    model._original_objective = objective
    model._penalty_weight = penalty_weight
    model._penalties = {
        'repair count': count_penalty,
        'repair persistence': persistence_penalty,
        'damaged-edge capacity': capacity_penalty,
        'flow conservation': conservation_penalty,
    }

    return model, recovery_vars, y_vars


def process_solution(model, n, broken_edges, recovery_vars, y_vars):
    if model.SolCount == 0:
        print('No QUBO solution found. Status:', model.Status)
        return

    if model.Status == GRB.OPTIMAL:
        print('QUBO solution is optimal')
    else:
        print('Gurobi status =', model.Status)

    print('penalty_weight =', model._penalty_weight)
    print('QUBO objective =', model.ObjVal)

    # Check feasibility before interpreting the repair sequence.
    penalty_values = {
        name: expression.getValue()
        for name, expression in model._penalties.items()
    }
    for name, value in penalty_values.items():
        print(f'{name} penalty =', value)

    if any(abs(value) > 1e-6 for value in penalty_values.values()):
        print('The returned QUBO solution violates the original constraints.')
        print('resilience loss before penalties =', model._original_objective.getValue())
        return

    k = len(broken_edges)
    recovery_seq = []
    for t in range(1, k + 1):
        recovered = [
            (u, v) for u, v in broken_edges
            if recovery_vars[u, v, t].X > 0.5
            and (t == 1 or recovery_vars[u, v, t - 1].X < 0.5)
        ]
        if len(recovered) != 1:
            raise ValueError(f'Expected one recovered edge at step {t}; found {recovered}')
        recovery_seq.append(edge_to_equipment[recovered[0]])

    for t in range(1, k + 1):
        flow = sum(y_vars[u, v, t].getValue() for u, v in output_edges)
        print(f'flow at iteration {t} = ', flow)

    print('objective = ', model._original_objective.getValue())
    print('recovery_seq = ', recovery_seq)


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('Usage: python3 recovery_edge_encoding_qubo_miqp.py input_file [-t seconds]')

    try:
        time_limit_s = 30
        if '-t' in sys.argv:
            time_limit_s = int(sys.argv[sys.argv.index('-t') + 1].strip())
        if time_limit_s <= 0:
            raise ValueError('Time limit must be positive.')

        n, source, sink, edge_list, edge_weights, broken_edges = read_input()
        model, recovery_vars, y_vars = build_model(
            n, source, sink, edge_list, edge_weights, broken_edges
        )

        model.Params.TimeLimit = time_limit_s
        model.Params.PreQLinearize = -1  # Automatic MIQP reformulation.
        model.Params.OutputFlag = False
        model.optimize()

        process_solution(model, n, broken_edges, recovery_vars, y_vars)
    except (OSError, ValueError, SyntaxError, IndexError, gp.GurobiError) as error:
        sys.exit(str(error))
