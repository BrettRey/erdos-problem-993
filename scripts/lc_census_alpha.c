/* lc_census_alpha.c — alpha-filtered LC census for trees, n <= 38
 * (2026-08-19; extends lc_census.c for the depth-3 dist-3 window).
 *
 * Reads gentreeg -p -q output (parent arrays, 1-indexed, parent[1]=0,
 * parent[i] < i) on stdin. For each tree, first computes the matching
 * number mu by the exact greedy leaf-matching (process v = n..2, match v
 * to par[v] iff both unmatched; children have larger indices, so this is
 * the classical exact tree algorithm), hence alpha = n - mu by Konig.
 * Only trees with alpha <= ALPHA_MAX get the full independence-polynomial
 * DP + LC/unimodality checks. Every DP'd tree cross-checks the greedy
 * alpha against deg I(T); a mismatch aborts (exit 4), so a window run
 * validates its own filter on every passed tree.
 *
 * Overflow: coefficients are independence counts <= C(38,19) < 2^35, and
 * every partial product fi*merged[j] counts a subset of independent sets,
 * so it is <= a final coefficient and fits uint64. Only the Turan
 * products c[k-1]*c[k+1] (up to ~2^69) need 128-bit.
 *
 * Build:  cc -O3 -march=native -o lc_census_alpha lc_census_alpha.c
 * Usage:  gentreeg -p -q 35 0/64 | ./lc_census_alpha 35 19
 *         (second arg = ALPHA_MAX; omit for no filter)
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 39
#define MAXC (MAXN + 1)

static int n;                      /* vertices per tree */
static int alpha_max;              /* filter threshold */
static int par[MAXN];              /* parent array, 1-indexed */
static uint64_t F[MAXN][MAXC];     /* poly, vertex excluded */
static uint64_t G[MAXN][MAXC];     /* poly, vertex included */
static int lf[MAXN], lg[MAXN];     /* lengths (degree+1) */

static unsigned long long trees = 0, passed = 0, failures = 0, alarms = 0;
static unsigned long long pass_by_alpha[MAXN];
static int verbose = 0;

static void report(const char *tag, const uint64_t *c, int len)
{
    printf("%s n=%d alpha=%d par=", tag, n, len - 1);
    for (int i = 1; i <= n; i++) printf(i == 1 ? "%d" : ",%d", par[i]);
    printf(" poly=");
    for (int i = 0; i < len; i++)
        printf(i == 0 ? "%llu" : ",%llu", (unsigned long long)c[i]);
    printf(" breaks=");
    int first = 1;
    for (int k = 1; k < len - 1; k++)
        if ((__uint128_t)c[k - 1] * c[k + 1] > (__uint128_t)c[k] * c[k]) {
            int thr = (2 * (len - 1) + 1) / 3;   /* ceil((2a-1)/3) */
            printf(first ? "%d(dist%d)" : ";%d(dist%d)", k, k - thr);
            first = 0;
        }
    printf("\n");
    fflush(stdout);
}

static void process_tree(void)
{
    trees++;

    /* exact greedy matching: children (larger indices) before parents */
    static unsigned char matched[MAXN];
    memset(matched, 0, (size_t)(n + 1));
    int mu = 0;
    for (int v = n; v >= 2; v--) {
        int p = par[v];
        if (!matched[v] && !matched[p]) { matched[v] = matched[p] = 1; mu++; }
    }
    int alpha = n - mu;              /* Konig: trees are bipartite */
    if (alpha > alpha_max) return;
    passed++;
    pass_by_alpha[alpha]++;

    for (int v = 1; v <= n; v++) {
        F[v][0] = 1; lf[v] = 1;
        G[v][0] = 0; G[v][1] = 1; lg[v] = 2;
    }
    /* children have larger indices: fold v into par[v], v = n..2 */
    for (int v = n; v >= 2; v--) {
        int p = par[v];
        uint64_t merged[MAXC], tmp[MAXC];
        int lm = lf[v] > lg[v] ? lf[v] : lg[v];
        for (int i = 0; i < lm; i++)
            merged[i] = (i < lf[v] ? F[v][i] : 0) + (i < lg[v] ? G[v][i] : 0);
        /* F[p] *= merged */
        int lt = lf[p] + lm - 1;
        memset(tmp, 0, (size_t)lt * sizeof(uint64_t));
        for (int i = 0; i < lf[p]; i++) {
            uint64_t fi = F[p][i];
            if (!fi) continue;
            for (int j = 0; j < lm; j++) tmp[i + j] += fi * merged[j];
        }
        memcpy(F[p], tmp, (size_t)lt * sizeof(uint64_t));
        lf[p] = lt;
        /* G[p] *= F[v] */
        lt = lg[p] + lf[v] - 1;
        memset(tmp, 0, (size_t)lt * sizeof(uint64_t));
        for (int i = 0; i < lg[p]; i++) {
            uint64_t gi = G[p][i];
            if (!gi) continue;
            for (int j = 0; j < lf[v]; j++) tmp[i + j] += gi * F[v][j];
        }
        memcpy(G[p], tmp, (size_t)lt * sizeof(uint64_t));
        lg[p] = lt;
    }
    uint64_t c[MAXC];
    int len = lf[1] > lg[1] ? lf[1] : lg[1];
    for (int i = 0; i < len; i++)
        c[i] = (i < lf[1] ? F[1][i] : 0) + (i < lg[1] ? G[1][i] : 0);

    if (len - 1 != alpha) {          /* filter self-check, every tree */
        fprintf(stderr, "ALPHA MISMATCH greedy=%d dp=%d tree #%llu\n",
                alpha, len - 1, trees);
        report("ERROR", c, len);
        exit(4);
    }

    if (verbose) report("TREE", c, len);

    /* log-concavity (128-bit products: c_k up to C(38,19) ~ 2^34) */
    int broken = 0;
    for (int k = 1; k < len - 1; k++)
        if ((__uint128_t)c[k - 1] * c[k + 1] > (__uint128_t)c[k] * c[k])
            { broken = 1; break; }
    if (broken) { failures++; report("FAIL", c, len); }

    /* unimodality */
    int rising = 1, uni = 1;
    for (int i = 1; i < len && uni; i++) {
        if (rising) { if (c[i] < c[i - 1]) rising = 0; }
        else if (c[i] > c[i - 1]) uni = 0;
    }
    if (!uni) { alarms++; report("ALARM_NONUNIMODAL", c, len); }
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: lc_census_alpha n [alpha_max]\n"); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 38) { fprintf(stderr, "n out of range [2,38]\n"); return 2; }
    alpha_max = (argc > 2) ? atoi(argv[2]) : n;
    if (argc > 3 && strcmp(argv[3], "all") == 0) verbose = 1;

    static char buf[1 << 20];
    size_t got;
    int field = 0, val = -1;
    /* stream ints; every n-th completes one tree */
    while ((got = fread(buf, 1, sizeof buf, stdin)) > 0) {
        for (size_t i = 0; i < got; i++) {
            char ch = buf[i];
            if (ch >= '0' && ch <= '9') {
                val = (val < 0 ? 0 : val) * 10 + (ch - '0');
            } else {
                if (val >= 0) {
                    par[++field] = val;
                    val = -1;
                    if (field == n) { process_tree(); field = 0; }
                }
            }
        }
    }
    if (val >= 0) {
        par[++field] = val;
        if (field == n) { process_tree(); field = 0; }
    }
    if (field != 0) {
        fprintf(stderr, "TRUNCATED INPUT: %d leftover fields\n", field);
        return 3;
    }
    printf("STATS n=%d alpha_max=%d trees=%llu passed=%llu failures=%llu alarms=%llu\n",
           n, alpha_max, trees, passed, failures, alarms);
    for (int a = 0; a <= n; a++)
        if (pass_by_alpha[a])
            printf("PASS_ALPHA n=%d alpha=%d count=%llu\n",
                   n, a, pass_by_alpha[a]);
    return alarms ? 42 : 0;
}
