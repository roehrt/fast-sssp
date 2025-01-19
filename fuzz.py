from timeit import timeit

import networkx as nx
import random
import subprocess
import numpy

random.seed(42)
numpy.random.seed(42)
def gen(n, m):
    lo, hi = -10**9, 10**9
    while True:
        T = nx.dfs_tree(nx.random_labeled_tree(n), 0)
        for u, v in T.edges():
            T[u][v]['weight'] = random.randint(lo, hi)
        dist = [float('inf')] * n
        dist[0] = 0
        for u, v in nx.dfs_edges(T, 0):
            dist[v] = dist[u] + T[u][v]['weight']
        for _ in range(m - (n-1)):
            while True:
                u, v = random.sample(range(n), 2)
                if T.has_edge(u, v):
                    continue
                if not (lo <= dist[v] - dist[u] <= hi):
                    continue
                T.add_edge(u, v, weight=random.randint(max(lo, dist[v] - dist[u]), hi))
                break
        assert not nx.negative_edge_cycle(T)
        yield T, dist

for G, solution in gen(5*10**4, 10**5):
    with open('fuzz.in', 'w') as f:
        f.write(f'{G.number_of_nodes()} {G.number_of_edges()}\n')
        for u, v in G.edges():
            f.write(f'{u+1} {v+1} {G[u][v]["weight"]}\n')
        f.write(' '.join(map(str, solution)))
    print('running...')
    with open('fuzz.in', 'r') as f:
        subprocess.run('./cmake-build-release/sssp', check=True, stdin=f)
    print('ok')

