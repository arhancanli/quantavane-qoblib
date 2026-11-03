/*
 * qob_topo: low-diameter graphs for the QOBLIB topology class (order/degree problem: n nodes, maximum degree d,
 * minimise the diameter; github.com/ZIB-AOPT/QOBLIB, 10-topology).
 *
 * Method. Degrees are fixed at d (one node at d-1 when n*d is odd): the start is a circulant graph (offsets
 * 1..floor(d/2), plus the antipodal matching for odd d), scrambled by 20*|E| random degree-preserving switches (more
 * if the result is disconnected). Its diameter is the first incumbent. Then, for target D = incumbent - 1, simulated
 * annealing moves by switches: two edges (a,b),(c,e) become (a,c),(b,e) or (a,e),(b,c) when the new edges are new and
 * not loops, so degrees never change. The cost of a graph comes from the balls R_k(v) = nodes within distance k of v,
 * kept as bitsets and built level by level (R_{k+1}(v) = R_k(v) | OR over neighbours u of R_k(u)):
 *     cost = WF * far + near,  far = sum_v (n - |R_D(v)|)   (ordered pairs further apart than D),
 *                              near = sum_{k=1}^{D-1} sum_v (n - |R_k(v)|)   (rewards slack).
 * Each round cools geometrically from T0 to T1 over `sweep` moves, then reheats from the best graph of the round.
 * A graph with far = 0 has diameter <= D: it becomes the incumbent (diameter recomputed exactly) and the target drops
 * by one. Lower bound: the Moore bound. A graph of maximum degree d and diameter k has at most
 * 1 + d + d(d-1) + ... + d(d-1)^(k-1) nodes, so the diameter is at least LB = the smallest k whose bound reaches n
 * (LB = 1 only when d >= n - 1). The program computes LB, never targets a diameter below it, and stops as soon as the
 * incumbent equals LB: the graph is then proven optimal. Otherwise the run ends when its CPU time reaches `budget`
 * seconds, or after `stall` CPU seconds without a new incumbent. No known value is used.
 *
 * Output: OUT (DIMACS-like .gph as QOBLIB uses: "p edge n m", "e u v", 1-based), rewritten at every new incumbent;
 * TS, the run's objective time series in QOBLIB's format ([[{"Time": seconds, "Incumbent": diameter}, ...]]: the
 * start graph, every new incumbent, and the end of the run). Times are the process's own CPU time (user + system)
 * since the start, i.e. the run's time without any waiting for a core on a shared machine; the RESULT line on stdout
 * also gives wall-clock times.
 * Every result is re-checked by check.py (a separate Python code path: breadth-first search from every node).
 *
 * Usage: qob_topo n d budget stall seed OUT TS [sweep] [T0] [T1] [WF]   (defaults 10000000, 8, 0.5, 1)
 * Build: cc -O3 -o qob_topo qob_topo.c -lm
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, d, D, W;
static int *deg, *adj;           /* adj[v*d + i], i < deg[v] */
static uint64_t *R0, *R1;         /* n*W bitsets each */
static long long WFAR;

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline double urand(void) { return (rnd() >> 11) * (1.0 / 9007199254740992.0); }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double cpu(void) { return (double)clock() / CLOCKS_PER_SEC; }

static inline int has_edge(int a, int b) {
    for (int i = 0; i < deg[a]; i++) if (adj[a * d + i] == b) return 1;
    return 0;
}
static inline void replace_nb(int a, int old, int nw) {
    for (int i = 0; i < deg[a]; i++) if (adj[a * d + i] == old) { adj[a * d + i] = nw; return; }
    fprintf(stderr, "replace_nb: missing edge\n"); exit(3);
}
static void add_edge(int a, int b) { adj[a * d + deg[a]++] = b; adj[b * d + deg[b]++] = a; }

/* one level: nxt[v] = cur[v] | OR_{u in N(v)} cur[u]; returns the number of unreached ordered pairs */
static long long level(const uint64_t *cur, uint64_t *nxt) {
    long long unreached = 0;
    for (int v = 0; v < n; v++) {
        uint64_t *o = nxt + (size_t)v * W;
        const uint64_t *a = cur + (size_t)v * W;
        for (int w = 0; w < W; w++) o[w] = a[w];
        for (int i = 0; i < deg[v]; i++) {
            const uint64_t *b = cur + (size_t)adj[v * d + i] * W;
            for (int w = 0; w < W; w++) o[w] |= b[w];
        }
        int c = 0;
        for (int w = 0; w < W; w++) c += __builtin_popcountll(o[w]);
        unreached += n - c;
    }
    return unreached;
}
static void singletons(uint64_t *r) {
    memset(r, 0, sizeof(uint64_t) * (size_t)n * W);
    for (int v = 0; v < n; v++) r[(size_t)v * W + (v >> 6)] |= 1ULL << (v & 63);
}

/* cost for target D and the far count (pairs beyond D) */
static long long cost(long long *farp) {
    uint64_t *cur = R0, *nxt = R1;
    singletons(cur);
    long long nearsum = 0, far = 0;
    for (int k = 1; k <= D; k++) {
        long long unreached = level(cur, nxt);
        uint64_t *t = cur; cur = nxt; nxt = t;
        if (k < D) nearsum += unreached; else far = unreached;
    }
    *farp = far;
    return WFAR * far + nearsum;
}

/* Moore lower bound on the diameter: smallest k >= 1 with 1 + d * sum_{i<k} (d-1)^i >= n */
static int moore_lb(int nn, int dd) {
    double M = 1, term = dd;
    for (int k = 1; k < nn; k++) { M += term; if (M >= nn) return k; term *= (dd - 1); }
    return nn - 1;
}

/* exact diameter, or -1 when the graph is disconnected */
static int diameter(void) {
    uint64_t *cur = R0, *nxt = R1;
    singletons(cur);
    long long prev = -1;
    for (int k = 1; k <= n; k++) {
        long long unreached = level(cur, nxt);
        uint64_t *t = cur; cur = nxt; nxt = t;
        if (unreached == 0) return k;
        if (unreached == prev) return -1;
        prev = unreached;
    }
    return -1;
}

static int ne;                    /* number of edges */
static int *ea, *eb;              /* edge list */

static void build_edges(void) {
    ne = 0;
    for (int v = 0; v < n; v++) for (int i = 0; i < deg[v]; i++) { int u = adj[v * d + i]; if (v < u) { ea[ne] = v; eb[ne] = u; ne++; } }
}

/* try switch of edges x,y; orientation o. returns 1 if applied */
static int try_switch(int x, int y, int o, int *sa, int *sb, int *sc, int *se) {
    int a = ea[x], b = eb[x], c = ea[y], e = eb[y];
    if (o) { int t = c; c = e; e = t; }
    /* new edges (a,c), (b,e) */
    if (a == c || a == e || b == c || b == e) return 0;
    if (has_edge(a, c) || has_edge(b, e)) return 0;
    replace_nb(a, b, c); replace_nb(b, a, e); replace_nb(c, e, a); replace_nb(e, c, b);
    ea[x] = a; eb[x] = c; ea[y] = b; eb[y] = e;
    *sa = a; *sb = b; *sc = c; *se = e;
    return 1;
}
static void undo_switch(int x, int y, int a, int b, int c, int e) {
    /* current edges (a,c) at x and (b,e) at y; restore (a,b), (c,e) */
    replace_nb(a, c, b); replace_nb(c, a, e); replace_nb(b, e, a); replace_nb(e, b, c);
    ea[x] = a; eb[x] = b; ea[y] = c; eb[y] = e;
}

static void write_graph(const char *path, int diam, unsigned long long seed, double wall) {
    char tmp[4200]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w");
    if (!f) { perror(tmp); exit(2); }
    fprintf(f, "c Undirected graph with diameter %d (n = %d, max degree %d)\n", diam, n, d);
    fprintf(f, "c Found by qob_topo (simulated annealing over degree-preserving edge switches), seed %llu, %.1f s CPU\n", seed, wall);
    fprintf(f, "p edge %d %d\n", n, ne);
    for (int i = 0; i < ne; i++) fprintf(f, "e %d %d\n", ea[i] + 1, eb[i] + 1);
    fclose(f); rename(tmp, path);
}

#define TSMAX 4096
static double ts_t[TSMAX]; static int ts_v[TSMAX]; static int ts_n;
static void ts_add(double t, int v) { if (ts_n < TSMAX) { ts_t[ts_n] = t; ts_v[ts_n] = v; ts_n++; } }
static void ts_write(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); exit(2); }
    fprintf(f, "[[");
    for (int i = 0; i < ts_n; i++) fprintf(f, "%s{\"Time\": %.3f, \"Incumbent\": %d}", i ? ", " : "", ts_t[i], ts_v[i]);
    fprintf(f, "]]\n");
    fclose(f);
}

int main(int argc, char **argv) {
    if (argc < 8) { fprintf(stderr, "usage: %s n d budget stall seed OUT TS [sweep] [T0] [T1] [WF]\n", argv[0]); return 2; }
    n = atoi(argv[1]); d = atoi(argv[2]);
    double budget = atof(argv[3]), stall = atof(argv[4]);
    unsigned long long seed = strtoull(argv[5], 0, 10);
    const char *out = argv[6], *tspath = argv[7];
    long long sweep = argc > 8 ? atoll(argv[8]) : 10000000LL;
    double T0 = argc > 9 ? atof(argv[9]) : 8.0, T1 = argc > 10 ? atof(argv[10]) : 0.5;
    if (n < 3 || d < 2 || d >= n || sweep < 1) { fprintf(stderr, "bad n d sweep\n"); return 2; }
    rs = seed * 0x9E3779B97F4A7C15ULL + 12345; if (!rs) rs = 1;
    W = (n + 63) / 64;
    WFAR = argc > 11 ? atoll(argv[11]) : 1;
    deg = calloc(n, sizeof(int)); adj = malloc(sizeof(int) * (size_t)n * d);
    R0 = malloc(sizeof(uint64_t) * (size_t)n * W); R1 = malloc(sizeof(uint64_t) * (size_t)n * W);
    ea = malloc(sizeof(int) * (size_t)n * d); eb = malloc(sizeof(int) * (size_t)n * d);
    int *bestadj = malloc(sizeof(int) * (size_t)n * d), *bestdeg = malloc(sizeof(int) * n);
    if (!R0 || !R1 || !adj || !bestadj) { fprintf(stderr, "out of memory\n"); return 2; }

    /* circulant start: offsets 1..floor(d/2); odd d adds the antipodal matching (one node left at d-1 if n odd) */
    for (int v = 0; v < n; v++) for (int o = 1; o <= d / 2; o++) { int u = (v + o) % n; if (!has_edge(v, u)) add_edge(v, u); }
    if (d & 1) {
        int h = n / 2;
        for (int v = 0; v < h; v++) { int u = v + h; if (u < n && !has_edge(v, u) && deg[v] < d && deg[u] < d) add_edge(v, u); }
    }
    build_edges();
    double t0 = now(), c0 = cpu();
    int sa, sb, sc, se;
    int inc = -1;
    for (int tries = 0; tries < 50 && inc < 0; tries++) {      /* scramble until connected */
        for (long long i = 0; i < 20LL * ne; i++) {
            int x = rnd() % ne, y = rnd() % ne;
            if (x != y) try_switch(x, y, rnd() & 1, &sa, &sb, &sc, &se);
        }
        inc = diameter();
    }
    if (inc < 0) { fprintf(stderr, "start graph disconnected\n"); return 3; }
    int lb = moore_lb(n, d);
    write_graph(out, inc, seed, cpu() - c0);
    ts_add(cpu() - c0, inc);
    double clast = cpu() - c0, wbest = now() - t0, cbest = clast;
    int start = inc;
    long long moves = 0;
    fprintf(stderr, "qob_topo n=%d d=%d edges=%d seed=%llu start diameter %d, Moore lower bound %d\n", n, d, ne, seed, inc, lb);
    for (D = inc - 1; D >= lb; D = inc - 1) {
        long long far, cur = cost(&far), best = cur;
        memcpy(bestadj, adj, sizeof(int) * (size_t)n * d); memcpy(bestdeg, deg, sizeof(int) * n);
        int found = 0, stop = 0;
        while (!found && !stop) {
            double lr = log(T1 / T0) / (double)sweep;
            for (long long it = 0; it < sweep; it++) {
                double T = T0 * exp(lr * it);
                int x = rnd() % ne, y = rnd() % ne;
                if (x == y) continue;
                if (!try_switch(x, y, rnd() & 1, &sa, &sb, &sc, &se)) continue;
                long long f2, c2 = cost(&f2);
                moves++;
                long long dc = c2 - cur;
                if (dc <= 0 || urand() < exp(-(double)dc / T)) {
                    cur = c2; far = f2;
                    if (cur < best) { best = cur; memcpy(bestadj, adj, sizeof(int) * (size_t)n * d); memcpy(bestdeg, deg, sizeof(int) * n); }
                    if (far == 0) { found = 1; break; }
                } else undo_switch(x, y, sa, sb, sc, se);
                if ((moves & 15) == 0) { double c = cpu() - c0; if (c > budget || c - clast > stall) { stop = 1; break; } }
            }
            if (found || stop) break;
            memcpy(adj, bestadj, sizeof(int) * (size_t)n * d); memcpy(deg, bestdeg, sizeof(int) * n);  /* reheat */
            build_edges(); cur = cost(&far);
            double c = cpu() - c0; if (c > budget || c - clast > stall) stop = 1;
        }
        if (!found) break;
        int dd = diameter();
        if (dd < 0 || dd > D) { fprintf(stderr, "internal: diameter %d above target %d\n", dd, D); return 3; }
        inc = dd;
        wbest = now() - t0; cbest = cpu() - c0; clast = cbest;
        write_graph(out, inc, seed, cbest);
        ts_add(cbest, inc);
        fprintf(stderr, "diameter %d at %.2f s wall, %.2f s cpu, %lld moves\n", inc, wbest, cbest, moves);
        if (cbest > budget) break;
    }
    double wend = now() - t0, cend = cpu() - c0;
    ts_add(cend, inc);
    ts_write(tspath);
    printf("RESULT n=%d d=%d diameter=%d start=%d lb=%d proven=%d wall=%.2f cpu=%.2f best_wall=%.2f best_cpu=%.2f moves=%lld seed=%llu\n",
           n, d, inc, start, lb, inc == lb, wend, cend, wbest, cbest, moves, seed);
    return 0;
}
