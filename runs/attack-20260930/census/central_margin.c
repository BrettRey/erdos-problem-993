/* central_margin.c -- exact census of central-window log-concavity margins
 * for trees (2026-09-30).  Adapted from scripts/lc_census.c (DP unchanged).
 *
 * Reads gentreeg -p -q output (parent arrays, 1-indexed, parent[1]=0,
 * parent[i] < i) on stdin, computes the FULL independence polynomial
 * i_0..i_alpha of every tree in exact uint64, and for
 *
 *     q   = ceil((2 alpha - 1)/3) = (2 alpha + 1) / 3   (integer division)
 *     W5  = { k : max(1,ceil(n/5)) <= k <= q }
 *     W4  = { k : max(1,ceil(n/4)) <= k <= q }
 *     delta_k = (i_k^2 - i_{k-1} i_{k+1}) / i_k^2        (exact rational)
 *     m5(T) = min_{k in W5} delta_k,   m4(T) = min_{k in W4} delta_k
 *
 * records (exactly, by __int128 cross-multiplication):
 *   - the minimum of m5 and m4 over all trees, with a top-TOPK list of
 *     minimising trees (parent array, alpha, k, poly);
 *   - the number of trees with some delta_k < 0 (resp. == 0) inside W5, W4;
 *   - (diagnostic, double precision, clearly labelled FLOAT) for every k in
 *     W5 the fugacity lambda_k with hard-core mean k, V_k = Var at lambda_k,
 *     and the minimum over trees and k of V_k * delta_k, with its minimiser.
 *   - optional Bernoulli sample (hash of the parent array, rate 1/P) of trees,
 *     printed as SAMPLE lines with their polynomial, for downstream Python.
 *
 * Arithmetic safety: n <= 32 enforced, i_k <= C(32,16) < 2^30, so
 * i_k^2 and i_{k-1} i_{k+1} < 2^60; numerators |num| < 2^60 and
 * denominators < 2^60, so cross products < 2^120 fit in signed __int128.
 *
 * Build:  cc -O3 -march=native -o central_margin central_margin.c -lm
 * Usage:  gentreeg -p -q 20 0/5 | ./central_margin 20 [P]
 *         (P > 0: emit SAMPLE lines for trees with hash % P == 0)
 */
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 33
#define MAXC (MAXN + 1)
#define TOPK 12

typedef __int128 i128;

static int n;
static int par[MAXN];
static uint64_t F[MAXN][MAXC];
static uint64_t G[MAXN][MAXC];
static int lf[MAXN], lg[MAXN];
static uint64_t sample_P = 0;

static unsigned long long trees = 0;
static unsigned long long viol5 = 0, viol4 = 0, zero5 = 0, zero4 = 0;
static unsigned long long nonuni = 0, sampled = 0;

typedef struct {
    int valid;
    i128 num, den;            /* margin = num/den, den > 0 */
    int k, alpha;
    int par[MAXN];
    uint64_t poly[MAXC];
    int len;
} rec_t;

static rec_t top5[TOPK], top4[TOPK];

/* single exact-minimum records (margin num/den, minimiser) */
typedef struct {
    int valid;
    i128 num, den;
    int k, alpha;
    int par[MAXN];
} one_t;
static one_t perk[MAXC];      /* min delta_k over trees with k in W5, by absolute k */
static one_t pera[MAXC];      /* min m5(T) over trees with independence number alpha */
static one_t topq;            /* min delta_q over trees (k = q, the window top) */
static one_t newton;          /* min Newton ratio r_k over trees, k in W5 (exact) */
static unsigned long long newton_fail_trees = 0;   /* trees with some r_k < 1 in W5 */
static unsigned long long newton_fail_by_dist[MAXC]; /* first failing k: count by q-k */

/* FLOAT diagnostic: min V_k*delta_k binned by V_k */
#define NVB 5
static const double vbin_lo[NVB] = {0.0, 2.0, 4.0, 8.0, 16.0};
typedef struct {
    int valid; double vd, V, lambda, delta; int k, alpha; int par[MAXN];
} vdrec_t;
static vdrec_t vbin[NVB];

/* float diagnostic: min over trees of min_{k in W5} V_k*delta_k */
static double vd_min = 1e300;
static int vd_k = -1, vd_alpha = -1;
static int vd_par[MAXN];
static uint64_t vd_poly[MAXC];
static int vd_len = 0;
static double vd_lambda = 0, vd_V = 0, vd_delta = 0;

/* a < b ? for fractions with positive denominators */
static int frac_lt(i128 an, i128 ad, i128 bn, i128 bd) { return an * bd < bn * ad; }

static void top_insert(rec_t *top, i128 num, i128 den, int k, int alpha,
                       const uint64_t *c, int len)
{
    /* find the slot: list sorted ascending by margin; ties kept in arrival
       order (the first TOPK trees with the smallest margins). */
    int pos = TOPK;
    for (int i = 0; i < TOPK; i++) {
        if (!top[i].valid || frac_lt(num, den, top[i].num, top[i].den)) { pos = i; break; }
    }
    if (pos == TOPK) return;
    for (int i = TOPK - 1; i > pos; i--) top[i] = top[i - 1];
    rec_t *r = &top[pos];
    r->valid = 1; r->num = num; r->den = den; r->k = k; r->alpha = alpha;
    memcpy(r->par, par, sizeof par);
    memcpy(r->poly, c, (size_t)len * sizeof(uint64_t));
    r->len = len;
}

static void one_update(one_t *r, i128 num, i128 den, int k, int alpha)
{
    if (!r->valid || frac_lt(num, den, r->num, r->den)) {
        r->valid = 1; r->num = num; r->den = den; r->k = k; r->alpha = alpha;
        memcpy(r->par, par, sizeof par);
    }
}

static uint64_t tree_hash(void)
{
    uint64_t h = 0x9E3779B97F4A7C15ULL ^ (uint64_t)n;
    for (int i = 1; i <= n; i++) {
        h ^= (uint64_t)(par[i] + 1) * 0xBF58476D1CE4E5B9ULL;
        h = (h ^ (h >> 31)) * 0x94D049BB133111EBULL;
        h ^= h >> 29;
    }
    return h;
}

/* hard-core mean and variance of the size at t = log(lambda) */
static void moments(const double *lc, int len, double t, double *mu, double *var)
{
    double m = -1e300;
    double lw[MAXC];
    for (int j = 0; j < len; j++) { lw[j] = lc[j] + j * t; if (lw[j] > m) m = lw[j]; }
    double s0 = 0, s1 = 0, s2 = 0;
    for (int j = 0; j < len; j++) {
        double w = exp(lw[j] - m);
        s0 += w; s1 += j * w; s2 += (double)j * j * w;
    }
    *mu = s1 / s0;
    *var = s2 / s0 - (*mu) * (*mu);
}

/* solve mean(t) = k by safeguarded Newton (mean is increasing in t) */
static double solve_t(const double *lc, int len, double k, double t0, double *Vout)
{
    double lo = -60, hi = 60, t = t0, mu, V;
    for (int it = 0; it < 200; it++) {
        moments(lc, len, t, &mu, &V);
        if (fabs(mu - k) < 1e-12 * (1 + k)) break;
        if (mu < k) lo = t; else hi = t;
        double tn = t - (mu - k) / (V > 1e-300 ? V : 1e-300);
        if (!(tn > lo && tn < hi)) tn = 0.5 * (lo + hi);
        t = tn;
    }
    moments(lc, len, t, &mu, &V);
    *Vout = V;
    return t;
}

static void process_tree(void)
{
    trees++;
    for (int v = 1; v <= n; v++) {
        F[v][0] = 1; lf[v] = 1;
        G[v][0] = 0; G[v][1] = 1; lg[v] = 2;
    }
    for (int v = n; v >= 2; v--) {
        int p = par[v];
        uint64_t merged[MAXC], tmp[MAXC];
        int lm = lf[v] > lg[v] ? lf[v] : lg[v];
        for (int i = 0; i < lm; i++)
            merged[i] = (i < lf[v] ? F[v][i] : 0) + (i < lg[v] ? G[v][i] : 0);
        int lt = lf[p] + lm - 1;
        memset(tmp, 0, (size_t)lt * sizeof(uint64_t));
        for (int i = 0; i < lf[p]; i++) {
            uint64_t fi = F[p][i];
            if (!fi) continue;
            for (int j = 0; j < lm; j++) tmp[i + j] += fi * merged[j];
        }
        memcpy(F[p], tmp, (size_t)lt * sizeof(uint64_t));
        lf[p] = lt;
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
    while (len > 1 && c[len - 1] == 0) len--;
    int alpha = len - 1;

    /* sanity: i_0 = 1, i_1 = n */
    if (c[0] != 1 || c[1] != (uint64_t)n) { fprintf(stderr, "DP SANITY FAIL\n"); exit(5); }

    int q = (2 * alpha + 1) / 3;
    int L5 = (n + 4) / 5, L4 = (n + 3) / 4;
    if (L5 < 1) L5 = 1;
    if (L4 < 1) L4 = 1;
    if (q > alpha - 1) { fprintf(stderr, "q > alpha-1 at alpha=%d\n", alpha); exit(6); }

    int has5 = 0, has4 = 0, bad5 = 0, bad4 = 0, z5 = 0, z4 = 0;
    i128 m5n = 0, m5d = 1, m4n = 0, m4d = 1;
    int k5 = -1, k4 = -1, nfail = 0;
    for (int k = L5; k <= q; k++) {
        i128 den = (i128)c[k] * c[k];
        i128 num = den - (i128)c[k - 1] * c[k + 1];
        if (num < 0) bad5 = 1;
        if (num == 0) z5 = 1;
        if (!has5 || frac_lt(num, den, m5n, m5d)) { m5n = num; m5d = den; k5 = k; has5 = 1; }
        one_update(&perk[k], num, den, k, alpha);
        if (k == q) one_update(&topq, num, den, k, alpha);
        {   /* Newton ratio r_k = delta_k (k+1)(alpha-k+1)/(alpha+1); exact for n <= 28 */
            i128 rn = num * (i128)((k + 1) * (alpha - k + 1));
            i128 rd = den * (i128)(alpha + 1);
            one_update(&newton, rn, rd, k, alpha);
            if (rn < rd && !nfail) { nfail = 1; newton_fail_by_dist[q - k]++; }
        }
        if (k >= L4) {
            if (num < 0) bad4 = 1;
            if (num == 0) z4 = 1;
            if (!has4 || frac_lt(num, den, m4n, m4d)) { m4n = num; m4d = den; k4 = k; has4 = 1; }
        }
    }
    if (nfail) newton_fail_trees++;
    if (has5) one_update(&pera[alpha], m5n, m5d, k5, alpha);
    if (bad5) viol5++;
    if (bad4) viol4++;
    if (z5) zero5++;
    if (z4) zero4++;
    if (has5) top_insert(top5, m5n, m5d, k5, alpha, c, len);
    if (has4) top_insert(top4, m4n, m4d, k4, alpha, c, len);
    if (bad5) {
        printf("CENTRAL_VIOLATION n=%d alpha=%d par=", n, alpha);
        for (int i = 1; i <= n; i++) printf(i == 1 ? "%d" : ",%d", par[i]);
        printf("\n"); fflush(stdout);
    }

    /* unimodality alarm */
    int rising = 1, uni = 1;
    for (int i = 1; i < len && uni; i++) {
        if (rising) { if (c[i] < c[i - 1]) rising = 0; }
        else if (c[i] > c[i - 1]) uni = 0;
    }
    if (!uni) {
        nonuni++;
        printf("ALARM_NONUNIMODAL n=%d par=", n);
        for (int i = 1; i <= n; i++) printf(i == 1 ? "%d" : ",%d", par[i]);
        printf("\n"); fflush(stdout);
    }

    /* FLOAT diagnostic: V_k * delta_k on W5 */
    double lc[MAXC];
    for (int j = 0; j < len; j++) lc[j] = log((double)c[j]);
    double t = 0;
    for (int k = L5; k <= q; k++) {
        double V;
        t = solve_t(lc, len, (double)k, t, &V);
        i128 den = (i128)c[k] * c[k];
        i128 num = den - (i128)c[k - 1] * c[k + 1];
        double d = (double)num / (double)den;
        double vd = V * d;
        {
            int b = NVB - 1;
            while (b > 0 && V < vbin_lo[b]) b--;
            vdrec_t *r = &vbin[b];
            if (!r->valid || vd < r->vd) {
                r->valid = 1; r->vd = vd; r->V = V; r->lambda = exp(t); r->delta = d;
                r->k = k; r->alpha = alpha; memcpy(r->par, par, sizeof par);
            }
        }
        if (vd < vd_min) {
            vd_min = vd; vd_k = k; vd_alpha = alpha;
            memcpy(vd_par, par, sizeof par);
            memcpy(vd_poly, c, (size_t)len * sizeof(uint64_t));
            vd_len = len; vd_lambda = exp(t); vd_V = V; vd_delta = d;
        }
    }

    if (sample_P && tree_hash() % sample_P == 0) {
        sampled++;
        printf("SAMPLE n=%d alpha=%d par=", n, alpha);
        for (int i = 1; i <= n; i++) printf(i == 1 ? "%d" : ",%d", par[i]);
        printf(" poly=");
        for (int i = 0; i < len; i++) printf(i == 0 ? "%llu" : ",%llu", (unsigned long long)c[i]);
        printf("\n");
    }
}

static void print_i128(i128 x)
{
    char buf[64]; int p = 63; buf[p] = 0;
    int neg = x < 0;
    unsigned __int128 u = neg ? (unsigned __int128)(-x) : (unsigned __int128)x;
    if (u == 0) buf[--p] = '0';
    while (u) { buf[--p] = (char)('0' + (int)(u % 10)); u /= 10; }
    if (neg) buf[--p] = '-';
    fputs(buf + p, stdout);
}

static void dump_top(const char *tag, rec_t *top)
{
    for (int i = 0; i < TOPK; i++) {
        if (!top[i].valid) break;
        rec_t *r = &top[i];
        printf("%s rank=%d n=%d alpha=%d k=%d q=%d num=", tag, i, n, r->alpha, r->k, (2 * r->alpha + 1) / 3);
        print_i128(r->num);
        printf(" den=");
        print_i128(r->den);
        printf(" par=");
        for (int j = 1; j <= n; j++) printf(j == 1 ? "%d" : ",%d", r->par[j]);
        printf(" poly=");
        for (int j = 0; j < r->len; j++) printf(j == 0 ? "%llu" : ",%llu", (unsigned long long)r->poly[j]);
        printf("\n");
    }
}

static void dump_one(const char *tag, one_t *r)
{
    printf("%s n=%d alpha=%d k=%d q=%d num=", tag, n, r->alpha, r->k, (2 * r->alpha + 1) / 3);
    print_i128(r->num);
    printf(" den=");
    print_i128(r->den);
    printf(" par=");
    for (int j = 1; j <= n; j++) printf(j == 1 ? "%d" : ",%d", r->par[j]);
    printf("\n");
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: central_margin n [P]\n"); return 2; }
    n = atoi(argv[1]);
    /* n <= 28: i_k <= C(27,13) < 2^25, so i_k^2 < 2^50 and the Newton-ratio
       cross products (< 2^50 * 2^50 * 2^13) fit in signed __int128. */
    if (n < 4 || n > 28) { fprintf(stderr, "n out of range [4,28]\n"); return 2; }
    if (argc > 2) sample_P = strtoull(argv[2], NULL, 10);

    static char buf[1 << 20];
    size_t got;
    int field = 0, val = -1;
    while ((got = fread(buf, 1, sizeof buf, stdin)) > 0) {
        for (size_t i = 0; i < got; i++) {
            char ch = buf[i];
            if (ch >= '0' && ch <= '9') {
                val = (val < 0 ? 0 : val) * 10 + (ch - '0');
            } else if (val >= 0) {
                par[++field] = val;
                val = -1;
                if (field == n) { process_tree(); field = 0; }
            }
        }
    }
    if (val >= 0) {
        par[++field] = val;
        if (field == n) { process_tree(); field = 0; }
    }
    if (field != 0) { fprintf(stderr, "TRUNCATED INPUT: %d leftover fields\n", field); return 3; }

    dump_top("TOP5", top5);
    dump_top("TOP4", top4);
    printf("VDMIN_FLOAT n=%d value=%.12g k=%d alpha=%d q=%d lambda=%.12g V=%.12g delta=%.12g par=",
           n, vd_min, vd_k, vd_alpha, (2 * vd_alpha + 1) / 3, vd_lambda, vd_V, vd_delta);
    for (int j = 1; j <= n; j++) printf(j == 1 ? "%d" : ",%d", vd_par[j]);
    printf(" poly=");
    for (int j = 0; j < vd_len; j++) printf(j == 0 ? "%llu" : ",%llu", (unsigned long long)vd_poly[j]);
    printf("\n");
    for (int k = 0; k < MAXC; k++) if (perk[k].valid) dump_one("PERK", &perk[k]);
    for (int a = 0; a < MAXC; a++) if (pera[a].valid) dump_one("PERALPHA", &pera[a]);
    dump_one("TOPQ", &topq);
    dump_one("NEWTON", &newton);
    for (int d = 0; d < MAXC; d++) if (newton_fail_by_dist[d])
        printf("NEWTONFAIL_DIST n=%d q_minus_k=%d trees=%llu\n", n, d, newton_fail_by_dist[d]);
    for (int b = 0; b < NVB; b++) if (vbin[b].valid) {
        vdrec_t *r = &vbin[b];
        printf("VDBIN_FLOAT n=%d bin_lo=%g value=%.12g k=%d alpha=%d q=%d lambda=%.12g V=%.12g delta=%.12g par=",
               n, vbin_lo[b], r->vd, r->k, r->alpha, (2 * r->alpha + 1) / 3, r->lambda, r->V, r->delta);
        for (int j = 1; j <= n; j++) printf(j == 1 ? "%d" : ",%d", r->par[j]);
        printf("\n");
    }
    printf("STATS n=%d trees=%llu viol5=%llu viol4=%llu zero5=%llu zero4=%llu nonunimodal=%llu sampled=%llu newton_fail_trees=%llu\n",
           n, trees, viol5, viol4, zero5, zero4, nonuni, sampled, newton_fail_trees);
    return nonuni ? 42 : 0;
}
