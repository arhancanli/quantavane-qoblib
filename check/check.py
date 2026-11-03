#!/usr/bin/env python3
"""Independent checker for QOBLIB solutions found by the qob_* searches (qob_topo.c, qob_numvc.c and others).

A separate code path from the searches: plain Python, integers only, every objective recomputed from the solution
file by its definition, nothing taken from the search's own bookkeeping.

  labs  N SOLUTION [CLAIM]        LABS energy E = sum_{k=1}^{N-1} C_k^2, C_k = sum_i s_i s_{i+k}, s_i = +-1.
                                  SOLUTION: QOBLIB's 0/1 vector format (one bit per line, '#' lines ignored).
  topology N DEG SOLUTION [CLAIM] diameter of the graph ("p edge n m" + "e u v", 1-based; .gz/.xz read too);
                                  checks node count, no loops or repeated edges, max degree <= DEG, connected.
  mis GRAPH SOLUTION [CLAIM]      size of the stable set (1-based node indices, '#' lines ignored); checks every
                                  index in range, no repeats and no edge of GRAPH inside the set.
  network N DEMAND SOLUTION [CLAIM]  largest arc load of a network design ("z", "x#i#j", "f#k#i#j" lines; demand matrix
                                  = leading N x N block of DEMAND, 1000 flow units per demand unit); checks out- and
                                  in-degree 2 everywhere, flow conservation of every commodity k at every node i != k
                                  (inflow minus outflow, not counting arcs back into k, equals 1000 t_ki), no flow on
                                  absent arcs, and that the stated z equals the recomputed largest load.

Prints one line "OK <kind> <value> ..." and exits 0 when the file is valid (and, with CLAIM, when the value equals
CLAIM); prints "FAIL ..." and exits 1 otherwise.
"""
import gzip
import lzma
import sys
from collections import deque


def read_text(path):
    if path.endswith('.gz'):
        with gzip.open(path, 'rt') as f:
            return f.read()
    if path.endswith('.xz'):
        with lzma.open(path, 'rt') as f:
            return f.read()
    with open(path) as f:
        return f.read()


def fail(msg):
    print('FAIL ' + msg)
    sys.exit(1)


def labs(n, path):
    bits = []
    for line in read_text(path).splitlines():
        if line.startswith('#'):
            continue
        for ch in line:
            if ch in '01':
                bits.append(ch)
            elif not (ch.isspace() or ch in ',;'):
                fail('labs: unexpected character %r' % ch)
    if len(bits) != n:
        fail('labs: %d bits, expected %d' % (len(bits), n))
    s = [1 if b == '1' else -1 for b in bits]
    energy = 0
    for k in range(1, n):
        c = 0
        for i in range(n - k):
            c += s[i] * s[i + k]
        energy += c * c
    return energy, ''


def read_graph(path):
    n = m = None
    edges = []
    for line in read_text(path).splitlines():
        t = line.split()
        if not t or t[0] == 'c':
            continue
        if t[0] == 'p':
            n, m = int(t[2]), int(t[3])
        elif t[0] == 'e':
            edges.append((int(t[1]), int(t[2])))
        else:
            fail('graph: unexpected line %r' % line[:60])
    if n is None:
        fail('graph: no p line')
    return n, m, edges


def topology(n, deg, path):
    gn, gm, edges = read_graph(path)
    if gn != n:
        fail('topology: %d nodes, expected %d' % (gn, n))
    if gm != len(edges):
        fail('topology: p line says %d edges, file has %d' % (gm, len(edges)))
    adj = [set() for _ in range(n + 1)]
    for u, v in edges:
        if not (1 <= u <= n and 1 <= v <= n):
            fail('topology: node out of range in edge %d %d' % (u, v))
        if u == v:
            fail('topology: loop at %d' % u)
        if v in adj[u]:
            fail('topology: repeated edge %d %d' % (u, v))
        adj[u].add(v)
        adj[v].add(u)
    maxdeg = max(len(a) for a in adj[1:])
    if maxdeg > deg:
        fail('topology: max degree %d > %d' % (maxdeg, deg))
    diameter = 0
    for src in range(1, n + 1):
        dist = [-1] * (n + 1)
        dist[src] = 0
        q = deque([src])
        while q:
            x = q.popleft()
            for y in adj[x]:
                if dist[y] < 0:
                    dist[y] = dist[x] + 1
                    q.append(y)
        if min(dist[1:]) < 0:
            fail('topology: disconnected (from node %d)' % src)
        diameter = max(diameter, max(dist[1:]))
    return diameter, 'edges=%d maxdeg=%d' % (len(edges), maxdeg)


def mis(graph_path, path):
    n, _, edges = read_graph(graph_path)
    chosen = []
    for line in read_text(path).splitlines():
        if line.startswith('#') or not line.strip():
            continue
        for tok in line.split():
            v = int(tok)
            if not 1 <= v <= n:
                fail('mis: index %d out of 1..%d' % (v, n))
            chosen.append(v)
    if not chosen:
        fail('mis: empty solution')
    inset = set(chosen)
    if len(inset) != len(chosen):
        fail('mis: repeated index')
    for u, v in edges:
        if u in inset and v in inset:
            fail('mis: edge %d %d inside the set' % (u, v))
    return len(inset), 'n=%d edges=%d' % (n, len(edges))


def network(n, demand_path, path):
    t = {}
    for line in read_text(demand_path).splitlines():
        parts = line.strip().split('|')
        if len(parts) < 3 or not parts[1].strip().isdigit():
            continue
        r = int(parts[1].strip())
        if 1 <= r <= n:
            vals = [int(v) for v in parts[2].split(',')[:n]]
            for c in range(n):
                t[(r, c + 1)] = vals[c]
    if len(t) != n * n:
        fail('network: demand block incomplete')
    x, f, z = {}, {}, None
    for line in read_text(path).splitlines():
        line = line.strip()
        if not line or line.startswith('#'):
            continue
        name, val = line.split()
        val = int(val)
        if val < 0:
            fail('network: negative value %s' % line)
        if name == 'z':
            z = val
        elif name.startswith('x#'):
            i, j = map(int, name[2:].split('#'))
            x[(i, j)] = val
        elif name.startswith('f#'):
            k, i, j = map(int, name[2:].split('#'))
            f[(k, i, j)] = val
        else:
            fail('network: unknown line %r' % line)
    nodes = range(1, n + 1)
    for i in nodes:
        if sum(x.get((i, j), 0) for j in nodes if j != i) != 2:
            fail('network: out-degree of %d is not 2' % i)
        if sum(x.get((j, i), 0) for j in nodes if j != i) != 2:
            fail('network: in-degree of %d is not 2' % i)
    for (k, i, j), v in f.items():
        if v and not x.get((i, j), 0):
            fail('network: flow on absent arc %d->%d' % (i, j))
    for k in nodes:
        for i in nodes:
            if i == k:
                continue
            inflow = sum(f.get((k, j, i), 0) for j in nodes if j != i)
            outflow = sum(f.get((k, i, j), 0) for j in nodes if j != i and j != k)
            if inflow - outflow != 1000 * t[(k, i)]:
                fail('network: conservation fails for commodity %d at node %d' % (k, i))
    load = 0
    for i in nodes:
        for j in nodes:
            if i != j:
                load = max(load, sum(f.get((k, i, j), 0) for k in nodes if k != j))
    if z != load:
        fail('network: stated z %s differs from the largest load %d' % (z, load))
    return load, 'arcs=%d' % sum(1 for v in x.values() if v)


def main(argv):
    if len(argv) < 2:
        print(__doc__)
        sys.exit(2)
    kind = argv[1]
    if kind == 'labs':
        value, info = labs(int(argv[2]), argv[3])
        claim = argv[4] if len(argv) > 4 else None
    elif kind == 'topology':
        value, info = topology(int(argv[2]), int(argv[3]), argv[4])
        claim = argv[5] if len(argv) > 5 else None
    elif kind == 'mis':
        value, info = mis(argv[2], argv[3])
        claim = argv[4] if len(argv) > 4 else None
    elif kind == 'network':
        value, info = network(int(argv[2]), argv[3], argv[4])
        claim = argv[5] if len(argv) > 5 else None
    else:
        print(__doc__)
        sys.exit(2)
    if claim is not None and int(claim) != value:
        fail('%s value %d differs from claim %s' % (kind, value, claim))
    print('OK %s %d %s' % (kind, value, info))


if __name__ == '__main__':
    main(sys.argv)
