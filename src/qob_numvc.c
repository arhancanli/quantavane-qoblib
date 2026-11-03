/*
 * qob_numvc: maximum independent set through minimum vertex cover (the complement of a cover is a stable set), for
 * the QOBLIB class 07-independentset (github.com/ZIB-AOPT/QOBLIB).
 *
 * Method: NuMVC (Cai, Su, Luo and Sattar 2013, "NuMVC: an efficient local search algorithm for minimum vertex
 * cover"): edge weights with forgetting, two-stage exchange and configuration checking.
 *   - Whenever the current set C covers every edge it is recorded and the vertex of C with the highest dscore is
 *     dropped, so the search then looks for a cover one smaller.
 *   - A step removes u in C with the highest dscore (oldest on ties; not the vertex added in the step before), then
 *     picks a random uncovered edge and adds its endpoint with confChange = 1 and the higher dscore (oldest on ties).
 *     Removing u clears confChange[u]; any change at a vertex sets confChange of its neighbours.
 *   - dscore(v) = weight of uncovered edges v would cover (v not in C), minus the weight that only v covers (v in C).
 *     Each step adds 1 to the weight of every uncovered edge; when the mean weight reaches gamma (default n/2) every
 *     weight becomes floor(rho * w) (rho = 0.3) and the scores are rebuilt.
 * Every reported set is re-checked here (each edge has an endpoint in the cover) and again by check.py
 * (a separate Python code path) before it is used.
 *
 * Output: "BEST size t=<wall s> cpu=<process CPU s> ..." lines (size = stable set = n - cover) on stdout; whenever the size exceeds `target`, OUT
 * is (re)written as a QOBLIB stable-set solution (1-based node indices, '#' comments).
 *
 * Usage: qob_numvc GRAPH.gph seconds seed target OUT [gamma] [rho] [stop] [TS]
 *   (`seconds` is a wall-clock limit; gamma, rho <= 0 keep the defaults; stop <= 0 means no early stop, otherwise
 *   the run ends once the stable set reaches `stop`; TS, when given, receives the run's objective time series in
 *   QOBLIB's format: [[{"Time": wall seconds since start, "Incumbent": stable set size}, ...]] with the first set,
 *   every improvement and the end of the run)
 * Build: cc -O3 -o qob_numvc qob_numvc.c
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

static int n, m;
static int *eu, *ev, *w;               /* edges and weights */
static int *xadj, *adjn, *adje;        /* CSR: neighbour and edge id */
static char *inC, *conf;
static long long *dscore, *stamp;
static int *unc, *upos, nunc;          /* uncovered edge list */
static int *cl, *cpos, ncl;            /* vertices in C (list for the scan) */

static uint64_t rs;
static inline uint64_t rnd(void) { rs ^= rs << 13; rs ^= rs >> 7; rs ^= rs << 17; return rs; }
static inline int rndint(int k) { return (int)((rnd() >> 33) % (uint64_t)k); }
static double now(void) { struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + 1e-9 * t.tv_nsec; }
static double cpu(void) { return (double)clock() / CLOCKS_PER_SEC; }

static inline void unc_add(int e) { upos[e] = nunc; unc[nunc++] = e; }
static inline void unc_del(int e) { int p = upos[e], last = unc[--nunc]; unc[p] = last; upos[last] = p; }
static inline void cl_add(int v) { cpos[v] = ncl; cl[ncl++] = v; }
static inline void cl_del(int v) { int p = cpos[v], last = cl[--ncl]; cl[p] = last; cpos[last] = p; }

static void add_v(int v, long long step) {
    inC[v] = 1; dscore[v] = -dscore[v]; stamp[v] = step; cl_add(v);
    for (int i = xadj[v]; i < xadj[v + 1]; i++) {
        int z = adjn[i], e = adje[i];
        if (!inC[z]) { dscore[z] -= w[e]; unc_del(e); }
        else dscore[z] += w[e];
        conf[z] = 1;
    }
}
static void remove_v(int v, long long step) {
    inC[v] = 0; dscore[v] = -dscore[v]; stamp[v] = step; conf[v] = 0; cl_del(v);
    for (int i = xadj[v]; i < xadj[v + 1]; i++) {
        int z = adjn[i], e = adje[i];
        if (!inC[z]) { dscore[z] += w[e]; unc_add(e); }
        else dscore[z] -= w[e];
        conf[z] = 1;
    }
}
static void rebuild_scores(void) {
    for (int v = 0; v < n; v++) dscore[v] = 0;
    for (int e = 0; e < m; e++) {
        int a = eu[e], b = ev[e];
        if (inC[a] && !inC[b]) dscore[a] -= w[e];
        else if (inC[b] && !inC[a]) dscore[b] -= w[e];
        else if (!inC[a] && !inC[b]) { dscore[a] += w[e]; dscore[b] += w[e]; }
    }
}

static int read_graph(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) { perror(path); return 0; }
    char line[256]; long long cap = 0; m = 0;
    while (fgets(line, sizeof line, f)) {
        if (line[0] == 'p') { char t[32]; long long mm; sscanf(line, "p %31s %d %lld", t, &n, &mm); cap = mm + 16; eu = malloc(sizeof(int) * cap); ev = malloc(sizeof(int) * cap); }
        else if (line[0] == 'e') {
            int a, b; if (sscanf(line + 1, "%d %d", &a, &b) != 2 || a == b) continue;
            if (m == cap) { cap = cap * 2 + 16; eu = realloc(eu, sizeof(int) * cap); ev = realloc(ev, sizeof(int) * cap); }
            eu[m] = a - 1; ev[m] = b - 1; m++;
        }
    }
    fclose(f);
    /* drop repeated edges: sort pairs (min,max) */
    long long *key = malloc(sizeof(long long) * (m + 1));
    for (int e = 0; e < m; e++) { int a = eu[e] < ev[e] ? eu[e] : ev[e], b = eu[e] < ev[e] ? ev[e] : eu[e]; key[e] = (long long)a * n + b; }
    int cmp(const void *x, const void *y);
    qsort(key, m, sizeof(long long), cmp);
    int k = 0;
    for (int e = 0; e < m; e++) if (e == 0 || key[e] != key[e - 1]) { eu[k] = (int)(key[e] / n); ev[k] = (int)(key[e] % n); k++; }
    m = k; free(key);
    int *dg = calloc(n + 1, sizeof(int));
    for (int e = 0; e < m; e++) { dg[eu[e]]++; dg[ev[e]]++; }
    xadj = malloc(sizeof(int) * (n + 1)); xadj[0] = 0;
    for (int v = 0; v < n; v++) xadj[v + 1] = xadj[v] + dg[v];
    adjn = malloc(sizeof(int) * (size_t)(2 * m + 1)); adje = malloc(sizeof(int) * (size_t)(2 * m + 1));
    int *fill = calloc(n, sizeof(int));
    for (int e = 0; e < m; e++) {
        int a = eu[e], b = ev[e];
        adjn[xadj[a] + fill[a]] = b; adje[xadj[a] + fill[a]++] = e;
        adjn[xadj[b] + fill[b]] = a; adje[xadj[b] + fill[b]++] = e;
    }
    free(dg); free(fill);
    return 1;
}
int cmp(const void *x, const void *y) { long long a = *(const long long *)x, b = *(const long long *)y; return a < b ? -1 : a > b; }

#define TSMAX 100000
static double ts_t[TSMAX]; static int ts_v[TSMAX]; static int ts_n;
static void ts_add(double t, int v) { if (ts_n < TSMAX) { ts_t[ts_n] = t; ts_v[ts_n] = v; ts_n++; } }
static void ts_write(const char *path) {
    FILE *f = fopen(path, "w");
    if (!f) { perror(path); return; }
    fprintf(f, "[[");
    for (int i = 0; i < ts_n; i++) fprintf(f, "%s{\"Time\": %.3f, \"Incumbent\": %d}", i ? ", " : "", ts_t[i], ts_v[i]);
    fprintf(f, "]]\n");
    fclose(f);
}

static int is_cover(const char *c) { for (int e = 0; e < m; e++) if (!c[eu[e]] && !c[ev[e]]) return 0; return 1; }

static void write_sol(const char *path, const char *cover, int size, const char *gname, unsigned long long seed, double t) {
    char tmp[1100]; snprintf(tmp, sizeof tmp, "%s.tmp", path);
    FILE *f = fopen(tmp, "w"); if (!f) { perror(tmp); return; }
    fprintf(f, "# Objective value = %d\n# Stable set for %s found by qob_numvc (vertex cover local search), seed %llu, %.1f s\n", size, gname, seed, t);
    for (int v = 0; v < n; v++) if (!cover[v]) fprintf(f, "%d\n", v + 1);
    fclose(f); rename(tmp, path);
}

int main(int argc, char **argv) {
    if (argc < 6) { fprintf(stderr, "usage: %s GRAPH seconds seed target OUT [gamma] [rho]\n", argv[0]); return 2; }
    const char *gpath = argv[1]; double secs = atof(argv[2]); unsigned long long seed = strtoull(argv[3], 0, 10);
    int target = atoi(argv[4]); const char *out = argv[5];
    rs = seed * 0x9E3779B97F4A7C15ULL + 0x2545F4914F6CDD1DULL; if (!rs) rs = 1;
    double t0 = now();
    if (!read_graph(gpath)) return 2;
    double gamma = argc > 6 ? atof(argv[6]) : 0, rho = argc > 7 ? atof(argv[7]) : 0;
    if (gamma <= 0) gamma = 0.5 * n;
    if (rho <= 0) rho = 0.3;
    int stopsize = argc > 8 ? atoi(argv[8]) : 0;
    if (stopsize <= 0) stopsize = n + 1;
    const char *tspath = argc > 9 ? argv[9] : 0;
    w = malloc(sizeof(int) * m); for (int e = 0; e < m; e++) w[e] = 1;
    inC = calloc(n, 1); conf = malloc(n); memset(conf, 1, n);
    dscore = calloc(n, sizeof(long long)); stamp = calloc(n, sizeof(long long));
    unc = malloc(sizeof(int) * (m + 1)); upos = malloc(sizeof(int) * (m + 1)); nunc = 0;
    cl = malloc(sizeof(int) * n); cpos = malloc(sizeof(int) * n); ncl = 0;
    char *best = malloc(n);
    fprintf(stderr, "qob_numvc %s n=%d m=%d seed=%llu target=%d gamma=%.0f rho=%.2f (read %.1f s)\n", gpath, n, m, seed, target, gamma, rho, now() - t0);
    /* start: all edges uncovered, greedy cover by highest dscore with random ties */
    for (int e = 0; e < m; e++) { unc_add(e); dscore[eu[e]]++; dscore[ev[e]]++; }
    while (nunc > 0) {
        int e = unc[rndint(nunc)];
        int v = dscore[eu[e]] >= dscore[ev[e]] ? eu[e] : ev[e];
        add_v(v, 0);
    }
    /* drop redundant vertices */
    for (int i = ncl - 1; i >= 0; i--) { int v = cl[i]; if (dscore[v] == 0) remove_v(v, 0); }
    memcpy(best, inC, n);
    int bestsize = n - ncl;
    printf("BEST %d t=%.2f step=0\n", bestsize, now() - t0); fflush(stdout);
    ts_add(now() - t0, bestsize);
    if (bestsize > target) write_sol(out, best, bestsize, gpath, seed, now() - t0);
    long long totw = m;  /* sum of weights */
    int lastadded = -1;
    long long step;
    for (step = 1;; step++) {
        if ((step & 1023) == 0 && now() - t0 > secs) break;
        if (nunc == 0) {
            if (n - ncl > bestsize) {
                if (!is_cover(inC)) { fprintf(stderr, "internal: not a cover\n"); return 3; }
                memcpy(best, inC, n); bestsize = n - ncl;
                printf("BEST %d t=%.2f cpu=%.2f step=%lld\n", bestsize, now() - t0, cpu(), step); fflush(stdout);
                ts_add(now() - t0, bestsize);
                if (bestsize > target) write_sol(out, best, bestsize, gpath, seed, now() - t0);
                if (bestsize >= stopsize) break;
            }
            /* remove the vertex of C with the highest dscore */
            int u = -1;
            for (int i = 0; i < ncl; i++) { int v = cl[i]; if (u < 0 || dscore[v] > dscore[u] || (dscore[v] == dscore[u] && stamp[v] < stamp[u])) u = v; }
            remove_v(u, step);
            lastadded = -1;
            continue;
        }
        int u = -1;
        for (int i = 0; i < ncl; i++) {
            int v = cl[i]; if (v == lastadded) continue;
            if (u < 0 || dscore[v] > dscore[u] || (dscore[v] == dscore[u] && stamp[v] < stamp[u])) u = v;
        }
        if (u >= 0) remove_v(u, step);
        int e = unc[rndint(nunc)];
        int a = eu[e], b = ev[e], v;
        if (!conf[a]) v = b; else if (!conf[b]) v = a;
        else v = (dscore[a] > dscore[b] || (dscore[a] == dscore[b] && stamp[a] < stamp[b])) ? a : b;
        add_v(v, step); lastadded = v;
        for (int i = 0; i < nunc; i++) { int f = unc[i]; w[f]++; dscore[eu[f]]++; dscore[ev[f]]++; }
        totw += nunc;
        if ((double)totw >= gamma * m) {
            totw = 0;
            for (int f = 0; f < m; f++) { w[f] = (int)(rho * w[f]); if (w[f] < 1) w[f] = 1; totw += w[f]; }
            rebuild_scores();
        }
    }
    ts_add(now() - t0, bestsize);
    if (tspath) ts_write(tspath);
    printf("END best %d steps %lld time %.1f cpu %.1f\n", bestsize, step, now() - t0, cpu());
    return 0;
}
