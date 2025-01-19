/**
 * Author: Simon Lindholm
 * Date: 2015-02-23
 * License: CC0
 * Source: http://en.wikipedia.org/wiki/Bellman-Ford_algorithm
 * Description: Calculates shortest paths from $s$ in a graph that might have negative edge weights.
 * Unreachable nodes get dist = inf; nodes reachable through negative-weight cycles get dist = -inf.
 * Assumes $V^2 \max |w_i| < \tilde{} 2^{63}$.
 * Time: O(VE)
 * Status: Tested on kattis:shortestpath3
 */

#include <bits/stdc++.h>
#define all(x) begin(x), end(x)
#define sz(x) int(x.size())
#define rep(i,a,b) for (int i = (a); i < (b); ++i)

using namespace std;
using ll = long long;

const ll inf = LLONG_MAX;
struct Ed { int a, b, w, s() { return a < b ? a : -a; }};
struct Node { ll dist = inf; int prev = -1; };

void bellmanFord(vector<Node>& nodes, vector<Ed>& eds, int s) {
	nodes[s].dist = 0;
	sort(all(eds), [](Ed a, Ed b) { return a.s() < b.s(); });

	int lim = sz(nodes) / 2 + 2; // /3+100 with shuffled vertices
	rep(i,0,lim) for (Ed ed : eds) {
		Node cur = nodes[ed.a], &dest = nodes[ed.b];
		if (abs(cur.dist) == inf) continue;
		ll d = cur.dist + ed.w;
		if (d < dest.dist) {
			dest.prev = ed.a;
			dest.dist = (i < lim-1 ? d : -inf);
		}
	}
	rep(i,0,lim) for (Ed e : eds) {
		if (nodes[e.a].dist == -inf)
			nodes[e.b].dist = -inf;
	}
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cin.exceptions(std::istream::failbit);

    int n, m; std::cin >> n >> m;
    vector<Node> nodes(n);
    vector<Ed> eds(m);
    for (int i = 0; i < m; ++i) {
        int u, v, w; std::cin >> u >> v >> w;
        --u, --v;
        eds[i] = {u, v, w};
    }

    auto start = std::chrono::high_resolution_clock::now();
    bellmanFord(nodes, eds, 0);
    auto end = std::chrono::high_resolution_clock::now();
    std::cerr << "Execution time: " << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms\n";
}

