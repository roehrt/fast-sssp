import networkx as nx
import random
import subprocess

def gen():
    while True:
        G = nx.gnp_random_graph(100, 0.04, directed=True)
        for u, v in G.edges():
            G[u][v]['weight'] = random.randint(-3, 20)
        if nx.negative_edge_cycle(G):
            continue
        yield G, nx.single_source_bellman_ford_path_length(G, source=0)

for G, solution in gen():
    with open('fuzz.in', 'w') as f:
        f.write(f'{G.number_of_nodes()} {G.number_of_edges()}\n')
        for u, v in G.edges():
            f.write(f'{u+1} {v+1} {G[u][v]["weight"]}\n')
        f.write(' '.join(str(solution.get(u, (2**63-1) // 2)) for u in range(G.number_of_nodes())))
    print('running...')
    with open('fuzz.in', 'r') as f:
        subprocess.run('./cmake-build-release/sssp', check=True, stdin=f)
    print('ok')

