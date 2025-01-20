#include <bits/extc++.h>

std::mt19937 rng(42);

using len = int64_t;
len operator ""_l(unsigned long long x) { return static_cast<len>(x); }
const len inf = std::numeric_limits<len>::max() / 2;
using graph = std::vector<std::vector<std::pair<int, len>>>;
using priority_queue = __gnu_pbds::priority_queue<std::pair<len, int>, std::greater<>>;

graph transpose(const graph &G) {
    graph Gt(std::ssize(G));
    for (int u = 0; u < std::ssize(G); ++u)
        for (const auto &[v, w] : G[u])
            Gt[v].emplace_back(u, w);
    return Gt;
}

void remove_edges(graph &G, std::vector<std::pair<int, int>> edges) {
    std::ranges::sort(edges);
    for (int u = 0; u < std::ssize(G); ++u) {
        auto &adj = G[u];
        adj.erase(std::remove_if(adj.begin(), adj.end(), [&](const auto &p) {
            return std::ranges::binary_search(edges, std::pair(u, p.first));
        }), adj.end());
    }
}

auto strongly_connected_components(const graph &G) {
    int T = 0, cn = 0;
    std::vector<int> c(G.size()), b(G.size()), st;
    std::function<int(int)> f = [&] (int u) {
        int t = b[u] = ++T;
        st.push_back(u);
        for (auto [v, _] : G[u]) if (!c[v])
            t = std::min(t, b[v] ?: f(v));
        if (t == b[u] && ++cn) while (!c[u])
            c[st.back()] = cn, st.pop_back();
        return b[u] = t;
    };
    for (int u = 0; u < std::ssize(G); ++u) if (!c[u]) f(u);
    return std::pair{cn+1, c};
}

len minimum_weight(const graph &G) {
    len W = inf;
    for (int u = 0; u < std::ssize(G); ++u)
        for (const auto &[v, w] : G[u])
            W = std::min(W, w);
    return W;
}

void reweight(graph &G, const std::vector<len> &phi) {
    for (int u = 0; u < std::ssize(G); ++u)
        for (auto &[v, w] : G[u])
            w += phi[u] - phi[v];
}

std::vector<double> approximate_ball_sizes(const graph& G, len r, double epsilon) {
    const double MAGIC_CONSTANT = 0.01; // Used to be 5;
    const int samples = ceil(MAGIC_CONSTANT * log(G.size()) / std::pow(epsilon, 2));
    assert(samples > 0);

    std::vector<int> count(G.size());
    std::uniform_int_distribution<int> vertices(0, std::ssize(G) - 1);

    std::vector dist(G.size(), inf);
    priority_queue pq;
    std::vector<priority_queue::point_iterator> it(G.size());
    for (int i = 0; i < samples; ++i) {
        int s = vertices(rng);
        std::vector<int> out_ball;

        pq.push({dist[s] = 0, s});
        while (!pq.empty()) {
            auto [_, u] = pq.top();
            pq.pop();
            for (const auto &[v, w] : G[u]) {
                if (dist[u] + w > r) continue;
                if (dist[v] == inf) {
                    it[v] = pq.push({dist[v] = dist[u] + w, v});
                    out_ball.push_back(v);
                } else if (dist[v] > dist[u] + w) {
                    pq.modify(it[v], {dist[v] = dist[u] + w, v});
                }
            }
        }

        for (int u : out_ball) ++count[u];
        for (int u : out_ball) dist[u] = inf;
    }

    std::vector<double> ball_size(G.size());
    for (int i = 0; i < std::ssize(G); ++i) ball_size[i] = count[i] / (double)samples;
    return ball_size;
}

std::vector<bool> in_light_vertices(const graph& G, len d) {
    double epsilon = 1.0 / 8.0;
    auto ball_size = approximate_ball_sizes(G, d / 4, epsilon);
    std::vector<bool> light(G.size());
    for (int i = 0; i < std::ssize(G); ++i)
        light[i] = ball_size[i] <= 5.0 / 8.0;
    return light;
}

std::vector<std::pair<int, int>> decompose(graph H, len d) {
    for (int u = 0; u < std::ssize(H); ++u)
        for (auto &[v, w] : H[u])
            w = std::max(w, 0_l);

    auto Ht = transpose(H);

    auto in_light = in_light_vertices(H, d);
    auto out_light = in_light_vertices(Ht, d);

    double p = std::min(1.0 - std::numeric_limits<double>::epsilon(), 20 * log(H.size()) / d);
    std::geometric_distribution<len> geom(p);

    std::vector<std::pair<int, int>> S;
    std::vector<bool> is_deleted(H.size());

    std::vector dist(H.size(), inf);
    priority_queue pq;
    std::vector<priority_queue::point_iterator> it(H.size());
    auto carve = [&](auto &G, auto &light) {
        for (int i = 0; i < std::ssize(G); ++i) {
            if (!light[i] || is_deleted[i]) continue;
            len r = geom(rng);
            pq.push({dist[i] = 0, i});
            while (!pq.empty()) {
                auto [_, u] = pq.top();
                pq.pop();
                is_deleted[u] = true;
                for (const auto &[v, w] : G[u]) {
                    if (is_deleted[v]) continue;
                    if (dist[u] + w > r) {
                        if (!is_deleted[v]) S.emplace_back(u, v);
                        continue;
                    }
                    if (dist[v] == inf)
                        it[v] = pq.push({dist[v] = dist[u] + w, v});
                    else if (dist[v] > dist[u] + w)
                        pq.modify(it[v], {dist[v] = dist[u] + w, v});
                }
            }
        }
    };
    carve(H, out_light);
    carve(Ht, in_light);
    return S;
}

void fix_dag_edges(const graph &G, std::vector<len> &phi) {
    auto [cn, c] = strongly_connected_components(G);
    std::vector<std::vector<int>> scc(cn);
    for (int i = 0; i < std::ssize(G); ++i) scc[c[i]].push_back(i);
    std::vector<len> phi_prime(cn);
    for (auto &s : scc) {
        for (int u : s) {
            for (const auto &[v, w] : G[u]) {
                if (c[u] == c[v]) continue;
                len w_prime = w + phi[u] - phi[v];
                phi_prime[c[u]] = std::min(phi_prime[c[u]], phi_prime[c[v]] + w_prime);
            }
        }
    }
    for (int u = 0; u < std::ssize(G); ++u) phi[u] += phi_prime[c[u]];
}

void bellman_ford_dijkstra(graph G, std::vector<len> &phi) {
    reweight(G, phi);
    graph G_neg(std::ssize(G));
    for (int u = 0; u < std::ssize(G); ++u) {
        for (auto &[v, w]: G[u]) if (w < 0) G_neg[u].emplace_back(v, w);
        std::erase_if(G[u], [](const auto &p) { return p.second < 0; });
    }

    std::vector dist(std::ssize(G), 0_l);
    priority_queue pq;
    std::vector<priority_queue::point_iterator> it(std::ssize(G));
    for (int i = 0; i < std::ssize(G); ++i) it[i] = pq.push({0, i});

    while (!pq.empty()) {
        // Dijkstra phase
        std::vector<int> A;
        while (!pq.empty()) {
            auto [_, u] = pq.top();
            pq.pop();
            it[u] = pq.end();
            A.push_back(u);
            for (const auto &[v, w]: G[u]) {
                if (dist[v] > dist[u] + w) {
                    dist[v] = dist[u] + w;
                    if (it[v] != pq.end()) pq.modify(it[v], {dist[v], v});
                    else pq.push({dist[v], v});
                }
            }
        }

        // Bellman-Ford phase
        for (int u: A) {
            for (const auto &[v, w]: G_neg[u]) {
                if (dist[v] > dist[u] + w) {
                    dist[v] = dist[u] + w;
                    if (it[v] != pq.end()) pq.modify(it[v], {dist[v], v});
                    else pq.push({dist[v], v});
                }
            }
        }
    }

    for (int u = 0; u < std::ssize(G); ++u) phi[u] += dist[u];
}

std::vector<len> scale(graph G) {
    auto W = -minimum_weight(G);

    for (int u = 0; u < std::ssize(G); ++u)
        for (auto &[v, w] : G[u])
            w += (W + 1) / 2;

    std::function<std::vector<len>(graph, len)> dfs = [&W, &dfs](graph H, len d) -> std::vector<len> {
        if (H.size() <= 1 || d <= W / 2) return std::vector(std::ssize(H), 0_l);
        auto S = decompose(H, d);
        auto G = H;
        remove_edges(H, S);

        auto [cn, c] = strongly_connected_components(H);
        std::vector<std::vector<int>> scc(cn);
        for (int i = 0; i < std::ssize(G); ++i) scc[c[i]].push_back(i);

        std::vector<int> lookup(std::ssize(G));
        std::vector phi(std::ssize(G), 0_l);
        for (int i = 0; i < cn; ++i) {
            graph U(scc[i].size());
            for (int j = 0; j < std::ssize(scc[i]); ++j) lookup[scc[i][j]] = j;
            for (int j = 0; j < std::ssize(scc[i]); ++j) {
                int u = scc[i][j];
                for (const auto &[v, w] : G[u]) {
                    if (c[u] == c[v]) U[j].emplace_back(lookup[v], w);
                }
            }
            auto di = U.size() <= 3.0 / 4.0 * G.size() ? d : d / 2;
            auto psi = dfs(U, di);
            for (int j = 0; j < std::ssize(U); ++j) phi[scc[i][j]] = psi[j];
        }

        fix_dag_edges(H, phi);
        bellman_ford_dijkstra(G, phi);
        return phi;
    };

    for (;;) {
        auto phi = dfs(G, (G.size() * W + 1) / 2);
        auto H = G;
        reweight(H, phi);
        if (minimum_weight(H) >= 0) return phi;
    }
}

std::pair<std::vector<len>, std::vector<int>> dijkstra(const graph& G, int s) {
    std::vector<len> dist(G.size(), inf);
    std::vector<int> par(G.size(), -1);
    priority_queue pq;
    std::vector<priority_queue::point_iterator> it(std::ssize(G));
    pq.push({dist[s] = 0, par[s] = s});
    while (!pq.empty()) {
        auto [_, u] = pq.top();
        pq.pop();
        for (const auto &[v, w] : G[u]) {
            if (dist[v] == inf)
                it[v] = pq.push({dist[v] = dist[u] + w, v}), par[v] = u;
            else if (dist[v] > dist[u] + w)
                pq.modify(it[v], {dist[v] = dist[u] + w, v}), par[v] = u;
        }
    }
    return std::pair{dist, par};
}

std::pair<std::vector<len>, std::vector<int>> single_source_shortest_path(const graph &G, int s) {
    graph H = G;
    for (int u = 0; u < std::ssize(H); ++u)
        for (auto &[v, w] : H[u])
            w *= std::ssize(H);
    while (minimum_weight(H) < -1) {
        auto psi = scale(H);
        reweight(H, psi);
    }
    for (int u = 0; u < std::ssize(H); ++u)
        for (auto &[v, w] : H[u])
            w += 1;
    auto [_, par] = dijkstra(H, s);

    auto Gt = transpose(G);
    std::vector dist(G.size(), inf);
    dist[s] = 0;
    std::function<len(int)> dfs = [&](int u) {
        if (dist[u] != inf || par[u] == -1) return dist[u];
        len W = inf;
        for (const auto &[v, w] : Gt[u]) if (v == par[u]) W = std::min(W, w);
        return dist[u] = dfs(par[u]) + W;
    };
    for (int u = 0; u < std::ssize(G); ++u) dfs(u);
    return {dist, par};
}

int main() {
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    std::cin.exceptions(std::istream::failbit);

    int n, m; std::cin >> n >> m;
    graph G(n);
    for (int i = 0; i < m; ++i) {
        int u, v; len w; std::cin >> u >> v >> w;
        --u, --v;
        G[u].emplace_back(v, w);
    }

    auto [dist, par] = single_source_shortest_path(G, 0);
    for (int j = 0; j < n; ++j) std::cout << dist[j] << " \n"[j == n - 1];
    for (int j = 0; j < n; ++j) std::cout << par[j] + 1 << " \n"[j == n - 1];
}