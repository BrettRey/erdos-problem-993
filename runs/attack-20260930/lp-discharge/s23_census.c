/* r1_census.c -- exhaustive exact census of lemma R1 (closed-neighbourhood
 * averaging of the free-count defects) and of lemma PV, on the window
 * W(T) = [max(1,ceil(n/4)), q],  q = ceil((2 alpha - 1)/3) = (2 alpha + 1)/3.
 * (2026-09-30, Erdos #993 attack, lane r1-exhaustive.)
 *
 * Notation: j^v = i(T - N[v]),
 *   D_v(k) = k i_{k-1}(T) j^v_k - (k+1) i_k(T) j^v_{k-1}
 *   L_u(k) = sum_{v in N[u]} D_v(k) / (deg v + 1)
 * R1: L_u(k) <= 0 for all u, k in W.   PV: D_v(k) <= 0 for all v, k in W.
 *
 * Input: gentreeg -p -q n [res/mod] (parent arrays, 1-indexed, par[1] = 0,
 * par[i] < i), parsed exactly as in census/central_margin.c.
 *
 * Polynomials (all exact uint64, every coefficient counts independent sets of
 * an induced subforest of T, hence <= C(n, n/2) <= C(28,14) < 2^26):
 *   down DP (rooted at vertex 1, as in central_margin.c):
 *     F[v] = v excluded, subtree of v;  G[v] = v included;  A = F + G.
 *   up DP (v = 2..n, p = par v): the component of T - subtree(v) containing
 *   p, rooted at p:  UF[v] = p excluded, UG[v] = p included, UA = UF + UG,
 *   with UF[1] = UA[1] = 1 (empty outside):
 *     UF[v] = UA[p] * F[p] / A[v],   UG[v] = UF[p] * G[p] / F[v].
 *   Divisions are exact power-series divisions by polynomials with constant
 *   term 1: q[i] = P[i] - sum_{j>=1} D[j] q[i-j]; every partial value lies
 *   in [q[i], P[i]] since all terms are >= 0, so signed int64 is safe; we
 *   still assert q[i] >= 0 and that the division is exact (zero remainder).
 *   J_v = UF[v] * (G[v] / x).
 * Per-tree assertions (abort on failure):
 *   (a) k i_k(T) = sum_v j^v_{k-1}  for k = 1..alpha+1 (double count);
 *   (b) sum_v D_v(k) = k(k+1)(i_{k-1} i_{k+1} - i_k^2) for every window k;
 *   (c) sum_u LAM L_u(k) = LAM sum_v D_v(k).
 *
 * Exactness of L_u: LAM = lcm(1..n) is a common multiple of every deg v + 1
 * (<= n), so LL_u = LAM * L_u = sum_{v in N[u]} D_v * (LAM/(deg v+1)) is an
 * integer (this is equivalent to clearing by the per-u lcm; sign identical).
 * Headroom (signed __int128, max 2^127 - 1), n <= 28:
 *   i_k, j_k < 2^26;  k+1 <= 15 < 2^4;  |D_v| < 2^4 * 2^26 * 2^26 = 2^56;
 *   LAM <= lcm(1..28) = 80,313,433,200 < 2^37;  |D_v| * LAM/(d+1) < 2^93;
 *   |N[u]| <= 28 < 2^5 terms:  |LL_u| < 2^98;  sum_u |LL_u| < 2^103.
 *   PV ratio comparisons: lhs, rhs < 2^56, cross products < 2^112.
 *   (at n = 24: i_k <= C(24,12) = 2,704,156 < 2^22, LAM = 5,354,228,880
 *    < 2^33, |D_v| < 2^48, |LL_u| < 2^86.)
 * FLOAT (double) is used ONLY to rank tightest cases (labelled FLOAT in the
 * output); the exact values of every reported case are recomputed in Python
 * (Fraction) by merge.py.  Every sign/inequality claim (R1, PV) is exact.
 *
 * Build: cc -O3 -o r1_census r1_census.c
 * Usage: gentreeg -p -q n r/m | ./r1_census n [dump]
 */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXN 29
#define MAXC (MAXN + 2)
#define TOPK 20
#define MAXPRINT 3000

typedef __int128 i128;
typedef int64_t i64;

static int n, dump = 0;
static int par[MAXN];
static i64 F[MAXN][MAXC], G[MAXN][MAXC], A[MAXN][MAXC];
static int lf[MAXN], lg[MAXN], la[MAXN];
static i64 UF[MAXN][MAXC], UG[MAXN][MAXC], UA[MAXN][MAXC];
static int luf[MAXN], lug[MAXN], lua[MAXN];
static i64 J[MAXN][MAXC];
static int lj[MAXN];
static int deg[MAXN];
static i128 LAM;
static i128 W[MAXN];

static unsigned long long trees = 0, empty_window = 0, window_levels = 0, window_vk = 0;
static unsigned long long r1_viol = 0, r1_viol_trees = 0, r1_zero = 0, r1_zero_trees = 0;
static unsigned long long pv_viol = 0, pv_viol_trees = 0, printed = 0, abs_zero_levels = 0;
static unsigned long long r1_viol_by_deg[MAXN], pv_viol_by_deg[MAXN];

typedef struct {
    int valid;
    double val;     /* FLOAT ranking value */
    int k, u, degu, alpha, lo, q;
    i128 LL;        /* exact LAM * L_u */
    int par[MAXN];
} rec_t;
static rec_t top1[TOPK], top2[TOPK];

/* exact max PV ratio lhs/rhs in window (rhs > 0) */
static int pvmax_valid = 0;
static i128 pvmax_l = 0, pvmax_r = 1;
static int pvmax_k, pvmax_v, pvmax_deg, pvmax_alpha, pvmax_par[MAXN];

static void die(const char *msg)
{
    fprintf(stderr, "FATAL %s n=%d tree#%llu par=", msg, n, trees);
    for (int i = 1; i <= n; i++) fprintf(stderr, i == 1 ? "%d" : ",%d", par[i]);
    fprintf(stderr, "\n");
    exit(7);
}

static int mul(const i64 *a, int la_, const i64 *b, int lb, i64 *out)
{
    int lo = la_ + lb - 1;
    i64 tmp[2 * MAXC];
    memset(tmp, 0, sizeof(i64) * (size_t)lo);
    for (int i = 0; i < la_; i++) {
        i64 ai = a[i];
        if (!ai) continue;
        for (int j = 0; j < lb; j++) tmp[i + j] += ai * b[j];
    }
    memcpy(out, tmp, sizeof(i64) * (size_t)lo);
    return lo;
}

static int add(const i64 *a, int la_, const i64 *b, int lb, i64 *out)
{
    int lo = la_ > lb ? la_ : lb;
    for (int i = 0; i < lo; i++) out[i] = (i < la_ ? a[i] : 0) + (i < lb ? b[i] : 0);
    return lo;
}

/* exact division P / D, D[0] = 1; returns quotient length; aborts if inexact */
static int divx(const i64 *P, int lp, const i64 *D, int ld, i64 *q)
{
    if (D[0] != 1) die("divisor constant term != 1");
    int lq = lp - ld + 1;
    if (lq < 1) die("division length");
    i64 r[2 * MAXC];
    memcpy(r, P, sizeof(i64) * (size_t)lp);
    for (int i = 0; i < lq; i++) {
        i64 qi = r[i];
        if (qi < 0) die("negative quotient coefficient");
        q[i] = qi;
        if (qi) for (int j = 1; j < ld; j++) r[i + j] -= qi * D[j];
    }
    for (int i = lq; i < lp; i++) if (r[i] != 0) die("inexact division");
    return lq;
}

static void print_par(FILE *f, const int *p)
{
    for (int i = 1; i <= n; i++) fprintf(f, i == 1 ? "%d" : ",%d", p[i]);
}

static void print_i128(FILE *f, i128 x)
{
    char buf[64]; int p = 63; buf[p] = 0;
    int neg = x < 0;
    unsigned __int128 u = neg ? (unsigned __int128)(-x) : (unsigned __int128)x;
    if (u == 0) buf[--p] = '0';
    while (u) { buf[--p] = (char)('0' + (int)(u % 10)); u /= 10; }
    if (neg) buf[--p] = '-';
    fputs(buf + p, f);
}

static void top_insert(rec_t *top, double val, int k, int u, int alpha, int lo, int q, i128 LL)
{
    int pos = TOPK;
    for (int i = 0; i < TOPK; i++)
        if (!top[i].valid || val > top[i].val) { pos = i; break; }
    if (pos == TOPK) return;
    for (int i = TOPK - 1; i > pos; i--) top[i] = top[i - 1];
    rec_t *r = &top[pos];
    r->valid = 1; r->val = val; r->k = k; r->u = u; r->degu = deg[u];
    r->alpha = alpha; r->lo = lo; r->q = q; r->LL = LL;
    memcpy(r->par, par, sizeof par);
}

static void process_tree(void)
{
    trees++;
    /* degrees */
    for (int v = 1; v <= n; v++) deg[v] = (v == 1) ? 0 : 1;
    for (int v = 2; v <= n; v++) deg[par[v]]++;
    /* down DP */
    for (int v = 1; v <= n; v++) {
        F[v][0] = 1; lf[v] = 1;
        G[v][0] = 0; G[v][1] = 1; lg[v] = 2;
    }
    for (int v = n; v >= 2; v--) {
        int p = par[v];
        la[v] = add(F[v], lf[v], G[v], lg[v], A[v]);
        lf[p] = mul(F[p], lf[p], A[v], la[v], F[p]);
        lg[p] = mul(G[p], lg[p], F[v], lf[v], G[p]);
    }
    la[1] = add(F[1], lf[1], G[1], lg[1], A[1]);
    /* up DP */
    UF[1][0] = 1; luf[1] = 1; UA[1][0] = 1; lua[1] = 1; lug[1] = 0;
    for (int v = 2; v <= n; v++) {
        int p = par[v];
        i64 t[2 * MAXC];
        int lt = mul(UA[p], lua[p], F[p], lf[p], t);
        luf[v] = divx(t, lt, A[v], la[v], UF[v]);
        lt = mul(UF[p], luf[p], G[p], lg[p], t);
        lug[v] = divx(t, lt, F[v], lf[v], UG[v]);
        lua[v] = add(UF[v], luf[v], UG[v], lug[v], UA[v]);
    }
    /* J_v = UF[v] * G[v]/x */
    for (int v = 1; v <= n; v++) {
        if (G[v][0] != 0) die("G[v][0] != 0");
        lj[v] = mul(UF[v], luf[v], G[v] + 1, lg[v] - 1, J[v]);
        while (lj[v] > 1 && J[v][lj[v] - 1] == 0) lj[v]--;
    }
    const i64 *I = A[1];
    int len = la[1];
    while (len > 1 && I[len - 1] == 0) len--;
    int alpha = len - 1;
    if (I[0] != 1 || I[1] != n) die("DP sanity i0/i1");
#define CO(arr, l, k) (((k) >= 0 && (k) < (l)) ? (arr)[k] : 0)
    /* (a) double count k i_k = sum_v j^v_{k-1}, k = 1..alpha+1 */
    for (int k = 1; k <= alpha + 1; k++) {
        i64 s = 0;
        for (int v = 1; v <= n; v++) s += CO(J[v], lj[v], k - 1);
        if (s != (i64)k * CO(I, len, k)) die("double count (a)");
    }
    int q = (2 * alpha + 1) / 3;
    int lo = (n + 3) / 4;
    if (lo < 1) lo = 1;
    if (lo > q) { empty_window++; return; }

    int tree_r1 = 0, tree_zero = 0, tree_pv = 0;
    double best1 = -1e300, best2 = -1e300;
    int b1k = -1, b1u = -1, b2k = -1, b2u = -1;
    i128 b1LL = 0, b2LL = 0;
    for (int k = lo; k <= q; k++) {
        window_levels++;
        i128 D[MAXN], DW[MAXN], LL[MAXN];
        i128 ikm1 = CO(I, len, k - 1), ik = CO(I, len, k), ikp1 = CO(I, len, k + 1);
        i128 S = 0, absS = 0;
        for (int v = 1; v <= n; v++) {
            i128 lhs = (i128)k * ikm1 * CO(J[v], lj[v], k);
            i128 rhs = (i128)(k + 1) * ik * CO(J[v], lj[v], k - 1);
            D[v] = lhs - rhs;
            S += D[v];
            absS += D[v] < 0 ? -D[v] : D[v];
            window_vk++;
            if (D[v] > 0) {
                pv_viol++; tree_pv = 1;
                if (deg[v] < MAXN) pv_viol_by_deg[deg[v]]++;
                if (printed < MAXPRINT) {
                    printed++;
                    printf("PVVIOL n=%d alpha=%d lo=%d q=%d k=%d v=%d deg=%d lhs=", n, alpha, lo, q, k, v, deg[v]);
                    print_i128(stdout, lhs); printf(" rhs="); print_i128(stdout, rhs);
                    printf(" par="); print_par(stdout, par); printf("\n");
                }
            }
            if (rhs > 0 && (!pvmax_valid || lhs * pvmax_r > pvmax_l * rhs)) {
                pvmax_valid = 1; pvmax_l = lhs; pvmax_r = rhs; pvmax_k = k; pvmax_v = v;
                pvmax_deg = deg[v]; pvmax_alpha = alpha; memcpy(pvmax_par, par, sizeof par);
            }
            DW[v] = D[v] * W[deg[v]];   /* S23: neighbour share (1/3)/deg v, scaled by 3*LAM */
        }
        /* (b) */
        if (S != (i128)k * (k + 1) * (ikm1 * ikp1 - ik * ik)) die("sum identity (b)");
        if (absS == 0) abs_zero_levels++;
        for (int v = 1; v <= n; v++) LL[v] = 2 * LAM * D[v];   /* S23: own share 2/3, scaled by 3*LAM */
        for (int v = 2; v <= n; v++) { LL[v] += DW[par[v]]; LL[par[v]] += DW[v]; }
        /* (c) */
        i128 SL = 0;
        for (int u = 1; u <= n; u++) SL += LL[u];
        if (SL != 3 * LAM * S) die("sum identity (c) S23");
        double norm1 = 3.0 * (double)LAM * (double)k * (double)ikm1 * (double)ik;  /* FLOAT */
        double norm2 = absS > 0 ? (double)LAM * (double)absS / (double)n : 0.0;  /* FLOAT */
        for (int u = 1; u <= n; u++) {
            if (LL[u] > 0) {
                r1_viol++; tree_r1 = 1;
                if (deg[u] < MAXN) r1_viol_by_deg[deg[u]]++;
                if (printed < MAXPRINT) {
                    printed++;
                    printf("S23VIOL n=%d alpha=%d lo=%d q=%d k=%d u=%d deg=%d LAM=", n, alpha, lo, q, k, u, deg[u]);
                    print_i128(stdout, LAM); printf(" LL="); print_i128(stdout, LL[u]);
                    printf(" par="); print_par(stdout, par); printf("\n");
                }
            } else if (LL[u] == 0) {
                r1_zero++; tree_zero = 1;
            }
            double v1 = (double)LL[u] / norm1;  /* FLOAT ranking */
            if (v1 > best1) { best1 = v1; b1k = k; b1u = u; b1LL = LL[u]; }
            if (absS > 0) {
                double v2 = (double)LL[u] / norm2;
                if (v2 > best2) { best2 = v2; b2k = k; b2u = u; b2LL = LL[u]; }
            }
        }
        if (dump) {
            printf("DUMP k=%d par=", k); print_par(stdout, par);
            printf(" D=");
            for (int v = 1; v <= n; v++) { if (v > 1) printf(","); print_i128(stdout, D[v]); }
            printf(" LL=");
            for (int v = 1; v <= n; v++) { if (v > 1) printf(","); print_i128(stdout, LL[v]); }
            printf(" LAM="); print_i128(stdout, LAM);
            printf(" J=");
            for (int v = 1; v <= n; v++) {
                if (v > 1) printf(";");
                for (int i = 0; i < lj[v]; i++) printf(i ? ",%lld" : "%lld", (long long)J[v][i]);
            }
            printf(" I=");
            for (int i = 0; i < len; i++) printf(i ? ",%lld" : "%lld", (long long)I[i]);
            printf("\n");
        }
    }
    if (tree_r1) r1_viol_trees++;
    if (tree_zero) r1_zero_trees++;
    if (tree_pv) pv_viol_trees++;
    if (b1k >= 0) top_insert(top1, best1, b1k, b1u, alpha, lo, q, b1LL);
    if (b2k >= 0) top_insert(top2, best2, b2k, b2u, alpha, lo, q, b2LL);
}

static void dump_top(const char *tag, rec_t *top)
{
    for (int i = 0; i < TOPK; i++) {
        if (!top[i].valid) break;
        rec_t *r = &top[i];
        printf("%s rank=%d n=%d alpha=%d lo=%d q=%d k=%d u=%d degu=%d float=%.12g LL=", tag, i, n,
               r->alpha, r->lo, r->q, r->k, r->u, r->degu, r->val);
        print_i128(stdout, r->LL);
        printf(" par="); print_par(stdout, r->par); printf("\n");
    }
}

int main(int argc, char **argv)
{
    if (argc < 2) { fprintf(stderr, "usage: r1_census n [dump]\n"); return 2; }
    n = atoi(argv[1]);
    if (n < 2 || n > 28) { fprintf(stderr, "n out of range [2,28]\n"); return 2; }
    if (argc > 2 && strcmp(argv[2], "dump") == 0) dump = 1;
    /* LAM = lcm(1..n) */
    LAM = 1;
    for (int i = 2; i <= n; i++) {
        i128 a = LAM, b = i;
        while (b) { i128 t = a % b; a = b; b = t; }
        LAM = LAM / a * i;
    }
    for (int d = 1; d <= n; d++) W[d] = LAM / d;

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

    dump_top("TOP_NORM1", top1);
    dump_top("TOP_NORM2", top2);
    if (pvmax_valid) {
        printf("PVMAX n=%d alpha=%d k=%d v=%d deg=%d lhs=", n, pvmax_alpha, pvmax_k, pvmax_v, pvmax_deg);
        print_i128(stdout, pvmax_l); printf(" rhs="); print_i128(stdout, pvmax_r);
        printf(" par="); print_par(stdout, pvmax_par); printf("\n");
    }
    for (int d = 0; d < MAXN; d++) {
        if (r1_viol_by_deg[d]) printf("R1VIOL_BYDEG n=%d deg=%d count=%llu\n", n, d, r1_viol_by_deg[d]);
        if (pv_viol_by_deg[d]) printf("PVVIOL_BYDEG n=%d deg=%d count=%llu\n", n, d, pv_viol_by_deg[d]);
    }
    printf("STATS n=%d trees=%llu empty_window=%llu window_levels=%llu window_vk=%llu "
           "r1_viol=%llu r1_viol_trees=%llu r1_zero=%llu r1_zero_trees=%llu "
           "pv_viol=%llu pv_viol_trees=%llu abs_zero_levels=%llu LAM=",
           n, trees, empty_window, window_levels, window_vk, r1_viol, r1_viol_trees,
           r1_zero, r1_zero_trees, pv_viol, pv_viol_trees, abs_zero_levels);
    print_i128(stdout, LAM);
    printf("\n");
    return 0;
}
