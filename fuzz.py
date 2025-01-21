import networkx as nx
import random
import subprocess
import numpy

random.seed(42)
numpy.random.seed(42)

def gen(n, m):
    lo, hi = -10**9, 10**9
    while True:
        G = nx.dfs_tree(nx.random_labeled_tree(n), 0)
        for u, v in G.edges():
            G[u][v]['weight'] = random.randint(lo, hi)
        dist = [float('inf')] * n
        dist[0] = 0
        for u, v in nx.dfs_edges(G, 0):
            dist[v] = dist[u] + G[u][v]['weight']
        for _ in range(m - (n-1)):
            u, v = random.sample(range(n), 2)
            while G.has_edge(u, v) or not (lo <= dist[v] - dist[u] <= hi):
                u, v = random.sample(range(n), 2)
            G.add_edge(u, v, weight=random.randint(max(lo, dist[v] - dist[u]), hi))
        assert not nx.negative_edge_cycle(G)
        yield G, dist

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('executable')
    parser.add_argument('-n', type=int, default=100)
    parser.add_argument('-m', type=int, default=500)
    args = parser.parse_args()
    for i, (G, solution) in enumerate(gen(args.n, args.m), start=1):
        with open('fuzz.in', 'w') as f:
            f.write(f'{G.number_of_nodes()} {G.number_of_edges()}\n')
            for u, v in G.edges():
                f.write(f'{u+1} {v+1} {G[u][v]["weight"]}\n')
            f.write(' '.join(map(str, solution)))
        print('running...', end=' ', flush=True)
        with open('fuzz.in', 'r') as f:
            subprocess.run(args.executable, check=True, stdin=f, stdout=subprocess.DEVNULL)
        print(f'ok ({i})')

if __name__ == '__main__':
    main()
