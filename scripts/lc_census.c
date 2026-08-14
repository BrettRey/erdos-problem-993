/* lc_census.c — exact LC-failure census for trees (2026-08-14).
 *
 * Reads gentreeg -p -q output (parent arrays, 1-indexed, parent[1]=0,
 * parent[i] < i) on stdin, computes the FULL independence polynomial of
 * every tree in exact uint64 arithmetic, and reports:
 *   - every tree whose coefficient sequence has a log-concavity break
 *     (i_{k-1} i_{k+1} > i_k^2), with all break positions;
 *   - any non-unimodal sequence (a #993 counterexample) as ALARM;
 *   - a STATS trailer with tree/failure/alarm counts (completeness check
 *     against A000055 is the driver's job).
 *
 * Safety: n <= 32 enforced. i_k <= C(32,16) = 601,080,390 < 2^30, so
 * coefficients and the Turan products i_{k-1} i_{k+1} < 2^60 fit in uint64
 * with room to spare.
 *
 * Build:  cc -O3 -march=native -o lc_census lc_census.c
 * Usage:  gentreeg -p -q 29 0/64 | ./lc_census 29
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 33
#define MAXC (MAXN + 1)

static int n;                      /* vertices per tree */
static int par[MAXN];              /* parent array, 1-indexed */
static uint64_t F[MAXN][MAXC];     /* poly, vertex excluded */
static uint64_t G[MAXN][MAXC];     /* poly, vertex included */
static int lf[MAXN], lg[MAXN];     /* lengths (degree+1) */

static unsigned long long trees = 0, failures = 0, alarms = 0;
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
        if (c[k - 1] * c[k + 1] > c[k] * c[k]) {
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

    if (verbose) report("TREE", c, len);

    /* log-concavity */
    int broken = 0;
    for (int k = 1; k < len - 1; k++)
        if (c[k - 1] * c[k + 1] > c[k] * c[k]) { broken = 1; break; }
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
    if (argc < 2) { fprintf(stderr, "usage: lc_census n\n"); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 32) { fprintf(stderr, "n out of range [2,32]\n"); return 2; }
    if (argc > 2 && strcmp(argv[2], "all") == 0) verbose = 1;

    static char buf[1 << 20];
    size_t got, pos = 0;
    int field = 0, val = -1, neg = 0;
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
        (void)pos; (void)neg;
    }
    if (val >= 0) {
        par[++field] = val;
        if (field == n) { process_tree(); field = 0; }
    }
    if (field != 0) {
        fprintf(stderr, "TRUNCATED INPUT: %d leftover fields\n", field);
        return 3;
    }
    printf("STATS n=%d trees=%llu failures=%llu alarms=%llu\n",
           n, trees, failures, alarms);
    return alarms ? 42 : 0;
}
