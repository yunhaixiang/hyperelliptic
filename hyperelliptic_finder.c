#include <flint/flint.h>
#include <flint/nmod_poly.h>

#include <errno.h>
#include <limits.h>
#include <signal.h>
#include <stdbool.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_DEGREE 128
#define MAX_GENUS 63

typedef unsigned long ulong;

typedef struct {
    int len;
    ulong c[MAX_DEGREE + 1];
} Poly;

typedef struct {
    bool used;
    ulong *key;
    int value;
} HashEntry;

typedef struct {
    HashEntry *entries;
    size_t capacity;
    size_t size;
    int key_len;
} HashMap;

typedef struct {
    int index;
    ulong *f;
    long long middle_coefficient;
    int canonical_index;
    ulong *canonical_f;
    bool reused_l_polynomial;
} CurveResult;

typedef struct {
    int p;
    int g;
    int key_len;
    int max_curves;
    bool allow_nonmonic;
    bool use_hasse_witt;
    bool quiet;
    bool verbose;
    char reduction[16];
    char output_txt[PATH_MAX];
    char output_json[PATH_MAX];
} Config;

typedef struct {
    long long considered;
    long long skipped_by_reduction;
    long long checked;
    long long rejected_by_hasse_witt;
    long long rejected_by_early_l_coefficient;
} Stats;

typedef struct {
    int p;
    int degree;
    long long order;
    Poly modulus;
} Field;

static Config config;
static CurveResult *results = NULL;
static size_t result_count = 0;
static size_t result_capacity = 0;
static volatile sig_atomic_t interrupted = 0;
static char search_status[32] = "incomplete";
static bool complete_list = false;

static ulong mod_add(ulong a, ulong b, ulong p) { return (a + b) % p; }
static ulong mod_sub(ulong a, ulong b, ulong p) { return (a + p - (b % p)) % p; }
static ulong mod_mul(ulong a, ulong b, ulong p) { return (ulong)(((unsigned long long)a * b) % p); }

static ulong mod_pow(ulong a, unsigned long long e, ulong p) {
    ulong r = 1 % p;
    a %= p;
    while (e) {
        if (e & 1) r = mod_mul(r, a, p);
        a = mod_mul(a, a, p);
        e >>= 1;
    }
    return r;
}

static ulong mod_inv(ulong a, ulong p) {
    long long t = 0, new_t = 1;
    long long r = (long long)p, new_r = (long long)(a % p);
    while (new_r != 0) {
        long long q = r / new_r;
        long long tmp = t - q * new_t;
        t = new_t;
        new_t = tmp;
        tmp = r - q * new_r;
        r = new_r;
        new_r = tmp;
    }
    if (t < 0) t += (long long)p;
    return (ulong)t;
}

static bool is_prime_int(int n) {
    if (n < 2) return false;
    if (n == 2) return true;
    if (n % 2 == 0) return false;
    for (int d = 3; (long long)d * d <= n; d += 2) {
        if (n % d == 0) return false;
    }
    return true;
}

static long long ipow_ll(int base, int exp) {
    long long r = 1;
    for (int i = 0; i < exp; i++) {
        if (r > LLONG_MAX / base) {
            fprintf(stderr, "integer overflow computing %d^%d\n", base, exp);
            exit(1);
        }
        r *= base;
    }
    return r;
}

static void poly_normalize(Poly *a) {
    while (a->len > 0 && a->c[a->len - 1] == 0) a->len--;
}

static void poly_zero(Poly *a) {
    a->len = 0;
    memset(a->c, 0, sizeof(a->c));
}

static void poly_one(Poly *a, ulong p) {
    poly_zero(a);
    a->len = 1;
    a->c[0] = 1 % p;
}

static void poly_from_coeffs(Poly *a, const ulong *coeffs, int len, ulong p) {
    poly_zero(a);
    a->len = len;
    for (int i = 0; i < len; i++) a->c[i] = coeffs[i] % p;
    poly_normalize(a);
}

static void poly_add_to(Poly *out, const Poly *a, const Poly *b, ulong p) {
    Poly r;
    poly_zero(&r);
    r.len = a->len > b->len ? a->len : b->len;
    for (int i = 0; i < r.len; i++) {
        ulong av = i < a->len ? a->c[i] : 0;
        ulong bv = i < b->len ? b->c[i] : 0;
        r.c[i] = mod_add(av, bv, p);
    }
    poly_normalize(&r);
    *out = r;
}

static void poly_sub_to(Poly *out, const Poly *a, const Poly *b, ulong p) {
    Poly r;
    poly_zero(&r);
    r.len = a->len > b->len ? a->len : b->len;
    for (int i = 0; i < r.len; i++) {
        ulong av = i < a->len ? a->c[i] : 0;
        ulong bv = i < b->len ? b->c[i] : 0;
        r.c[i] = mod_sub(av, bv, p);
    }
    poly_normalize(&r);
    *out = r;
}

static void poly_mul_to(Poly *out, const Poly *a, const Poly *b, ulong p) {
    Poly r;
    poly_zero(&r);
    if (a->len == 0 || b->len == 0) {
        *out = r;
        return;
    }
    if (a->len + b->len - 2 > MAX_DEGREE) {
        fprintf(stderr, "polynomial degree exceeds C limit %d\n", MAX_DEGREE);
        exit(1);
    }
    r.len = a->len + b->len - 1;
    for (int i = 0; i < a->len; i++) {
        for (int j = 0; j < b->len; j++) {
            r.c[i + j] = (r.c[i + j] + mod_mul(a->c[i], b->c[j], p)) % p;
        }
    }
    poly_normalize(&r);
    *out = r;
}

static void poly_divmod(const Poly *a, const Poly *b, ulong p, Poly *q, Poly *rem) {
    poly_zero(q);
    *rem = *a;
    if (b->len == 0) {
        fprintf(stderr, "polynomial division by zero\n");
        exit(1);
    }
    if (rem->len < b->len) return;
    q->len = rem->len - b->len + 1;
    ulong inv_lead = mod_inv(b->c[b->len - 1], p);
    while (rem->len >= b->len && rem->len > 0) {
        int shift = rem->len - b->len;
        ulong coeff = mod_mul(rem->c[rem->len - 1], inv_lead, p);
        q->c[shift] = coeff;
        for (int i = 0; i < b->len; i++) {
            rem->c[shift + i] = mod_sub(rem->c[shift + i], mod_mul(coeff, b->c[i], p), p);
        }
        poly_normalize(rem);
    }
    poly_normalize(q);
}

static void poly_mod_to(Poly *out, const Poly *a, const Poly *modulus, ulong p) {
    Poly q, r;
    poly_divmod(a, modulus, p, &q, &r);
    *out = r;
}

static void poly_derivative_to(Poly *out, const Poly *a, ulong p) {
    poly_zero(out);
    if (a->len <= 1) return;
    out->len = a->len - 1;
    for (int i = 1; i < a->len; i++) out->c[i - 1] = mod_mul((ulong)i, a->c[i], p);
    poly_normalize(out);
}

static void poly_gcd_to(Poly *out, const Poly *a, const Poly *b, ulong p) {
    Poly x = *a, y = *b, r;
    while (y.len != 0) {
        poly_mod_to(&r, &x, &y, p);
        x = y;
        y = r;
    }
    if (x.len != 0) {
        ulong inv = mod_inv(x.c[x.len - 1], p);
        for (int i = 0; i < x.len; i++) x.c[i] = mod_mul(x.c[i], inv, p);
    }
    *out = x;
}

static void poly_pow_mod_to(Poly *out, const Poly *base, unsigned long long exponent, const Poly *modulus, ulong p) {
    Poly result, power, tmp;
    poly_one(&result, p);
    poly_mod_to(&power, base, modulus, p);
    while (exponent) {
        if (exponent & 1) {
            poly_mul_to(&tmp, &result, &power, p);
            poly_mod_to(&result, &tmp, modulus, p);
        }
        exponent >>= 1;
        if (exponent) {
            poly_mul_to(&tmp, &power, &power, p);
            poly_mod_to(&power, &tmp, modulus, p);
        }
    }
    *out = result;
}

static void poly_exact_div_to(Poly *out, const Poly *a, const Poly *b, ulong p) {
    Poly q, r;
    poly_divmod(a, b, p, &q, &r);
    if (r.len != 0) {
        fprintf(stderr, "expected exact polynomial division\n");
        exit(1);
    }
    *out = q;
}

static void poly_pow_plain_to(Poly *out, const Poly *base, int exponent, ulong p) {
    Poly result, power, tmp;
    poly_one(&result, p);
    power = *base;
    while (exponent) {
        if (exponent & 1) {
            poly_mul_to(&tmp, &result, &power, p);
            result = tmp;
        }
        exponent >>= 1;
        if (exponent) {
            poly_mul_to(&tmp, &power, &power, p);
            power = tmp;
        }
    }
    *out = result;
}

static bool has_square_root(ulong value, ulong p) {
    value %= p;
    return value == 0 || mod_pow(value, (p - 1) / 2, p) == 1;
}

static ulong nonsquare_representative(ulong p) {
    for (ulong v = 2; v < p; v++) {
        if (!has_square_root(v, p)) return v;
    }
    fprintf(stderr, "could not find nonsquare in F_%lu\n", p);
    exit(1);
}

static int key_degree(const ulong *key, int len) {
    for (int i = len - 1; i >= 0; i--) {
        if (key[i] != 0) return i;
    }
    return -1;
}

static bool normalize_leading_square_class(ulong *key, int key_len, ulong p, bool allow_nonmonic) {
    int degree = key_degree(key, key_len);
    if (degree < 0) return false;
    ulong leading = key[degree] % p;
    if (leading == 0) return false;
    ulong target;
    if (has_square_root(leading, p)) {
        target = 1;
    } else if (allow_nonmonic) {
        target = nonsquare_representative(p);
    } else {
        return false;
    }
    ulong factor = mod_mul(target, mod_inv(leading, p), p);
    for (int i = 0; i < key_len; i++) key[i] = mod_mul(key[i], factor, p);
    return true;
}

static bool poly_squarefree_key(const ulong *key, int key_len, ulong p) {
    Poly f, derivative, gcd;
    poly_from_coeffs(&f, key, key_len, p);
    poly_derivative_to(&derivative, &f, p);
    if (derivative.len == 0) return false;
    poly_gcd_to(&gcd, &f, &derivative, p);
    return gcd.len == 1 && gcd.c[0] == 1;
}

static uint64_t hash_key(const ulong *key, int len) {
    uint64_t h = 1469598103934665603ULL;
    for (int i = 0; i < len; i++) {
        h ^= (uint64_t)key[i] + 0x9e3779b97f4a7c15ULL + ((uint64_t)i << 32);
        h *= 1099511628211ULL;
    }
    return h;
}

static void map_init(HashMap *map, int key_len) {
    map->capacity = 1024;
    map->size = 0;
    map->key_len = key_len;
    map->entries = calloc(map->capacity, sizeof(HashEntry));
    if (!map->entries) {
        perror("calloc");
        exit(1);
    }
}

static bool key_equal(const ulong *a, const ulong *b, int len) {
    return memcmp(a, b, sizeof(ulong) * (size_t)len) == 0;
}

static void map_rehash(HashMap *map) {
    HashEntry *old = map->entries;
    size_t old_cap = map->capacity;
    map->capacity *= 2;
    map->entries = calloc(map->capacity, sizeof(HashEntry));
    if (!map->entries) {
        perror("calloc");
        exit(1);
    }
    map->size = 0;
    for (size_t i = 0; i < old_cap; i++) {
        if (!old[i].used) continue;
        uint64_t h = hash_key(old[i].key, map->key_len);
        for (size_t j = 0;; j++) {
            size_t pos = (h + j) & (map->capacity - 1);
            if (!map->entries[pos].used) {
                map->entries[pos] = old[i];
                map->size++;
                break;
            }
        }
    }
    free(old);
}

static int map_get(const HashMap *map, const ulong *key) {
    uint64_t h = hash_key(key, map->key_len);
    for (size_t j = 0; j < map->capacity; j++) {
        size_t pos = (h + j) & (map->capacity - 1);
        if (!map->entries[pos].used) return 0;
        if (key_equal(map->entries[pos].key, key, map->key_len)) return map->entries[pos].value;
    }
    return 0;
}

static void map_setdefault(HashMap *map, const ulong *key, int value) {
    if ((map->size + 1) * 10 > map->capacity * 7) map_rehash(map);
    uint64_t h = hash_key(key, map->key_len);
    for (size_t j = 0;; j++) {
        size_t pos = (h + j) & (map->capacity - 1);
        if (!map->entries[pos].used) {
            map->entries[pos].used = true;
            map->entries[pos].key = malloc(sizeof(ulong) * (size_t)map->key_len);
            if (!map->entries[pos].key) {
                perror("malloc");
                exit(1);
            }
            memcpy(map->entries[pos].key, key, sizeof(ulong) * (size_t)map->key_len);
            map->entries[pos].value = value;
            map->size++;
            return;
        }
        if (key_equal(map->entries[pos].key, key, map->key_len)) return;
    }
}

static void append_result(const ulong *f, long long middle, int canonical_index, const ulong *canonical_f, bool reused) {
    if (result_count == result_capacity) {
        result_capacity = result_capacity ? result_capacity * 2 : 64;
        results = realloc(results, result_capacity * sizeof(CurveResult));
        if (!results) {
            perror("realloc");
            exit(1);
        }
    }
    CurveResult *r = &results[result_count];
    r->index = (int)result_count + 1;
    r->f = malloc(sizeof(ulong) * (size_t)config.key_len);
    r->canonical_f = malloc(sizeof(ulong) * (size_t)config.key_len);
    if (!r->f || !r->canonical_f) {
        perror("malloc");
        exit(1);
    }
    memcpy(r->f, f, sizeof(ulong) * (size_t)config.key_len);
    memcpy(r->canonical_f, canonical_f, sizeof(ulong) * (size_t)config.key_len);
    r->middle_coefficient = middle;
    r->canonical_index = canonical_index;
    r->reused_l_polynomial = reused;
    result_count++;
}

static void format_polynomial(const ulong *key, int key_len, char *buf, size_t buflen) {
    buf[0] = '\0';
    bool first = true;
    for (int i = 0; i < key_len; i++) {
        ulong coeff = key[i];
        if (coeff == 0) continue;
        char term[64];
        if (i == 0) {
            snprintf(term, sizeof(term), "%lu", coeff);
        } else if (i == 1) {
            snprintf(term, sizeof(term), coeff == 1 ? "x" : "%lux", coeff);
        } else {
            if (coeff == 1) snprintf(term, sizeof(term), "x^%d", i);
            else snprintf(term, sizeof(term), "%lux^%d", coeff, i);
        }
        if (!first) strncat(buf, " + ", buflen - strlen(buf) - 1);
        strncat(buf, term, buflen - strlen(buf) - 1);
        first = false;
    }
    if (first) snprintf(buf, buflen, "0");
}

static void write_key_array(FILE *fp, const ulong *key, int key_len) {
    int degree = key_degree(key, key_len);
    if (degree < 0) degree = 0;
    fprintf(fp, "[");
    for (int i = 0; i <= degree; i++) {
        if (i) fprintf(fp, ", ");
        fprintf(fp, "%lu", key[i]);
    }
    fprintf(fp, "]");
}

static int reduction_class_count(void) {
    int count = 0;
    for (size_t i = 0; i < result_count; i++) {
        if (results[i].canonical_index == results[i].index) count++;
    }
    return count;
}

static void output_paths(const char *input, char *txt, size_t txt_len, char *json, size_t json_len) {
    size_t n = strlen(input);
    if (n >= 5 && strcmp(input + n - 5, ".json") == 0) {
        snprintf(json, json_len, "%s", input);
        snprintf(txt, txt_len, "%.*s.txt", (int)(n - 5), input);
    } else if (n >= 4 && strcmp(input + n - 4, ".txt") == 0) {
        snprintf(txt, txt_len, "%s", input);
        snprintf(json, json_len, "%.*s.json", (int)(n - 4), input);
    } else {
        snprintf(txt, txt_len, "%s", input);
        snprintf(json, json_len, "%s.json", input);
    }
}

static void write_results(void) {
    FILE *txt = fopen(config.output_txt, "w");
    if (!txt) {
        perror(config.output_txt);
        return;
    }
    int classes = reduction_class_count();
    fprintf(txt, "Hyperelliptic curves with trinomial L-polynomial\n");
    fprintf(txt, "========================================================\n");
    fprintf(txt, "p = %d\n", config.p);
    fprintf(txt, "g = %d\n", config.g);
    fprintf(txt, "base field: F_%d\n", config.p);
    fprintf(txt, "presentation reduction: %s\n", config.reduction);
    fprintf(txt, "monic enforced: %s\n", config.allow_nonmonic ? "False" : "True");
    fprintf(txt, "nonmonic leading square class included: %s\n", config.allow_nonmonic ? "True" : "False");
    fprintf(txt, "hasse-witt prefilter: %s\n", config.use_hasse_witt ? "True" : "False");
    fprintf(txt, "reduction class count: %d\n", classes);
    fprintf(txt, "isomorphism class count: %d\n", classes);
    fprintf(txt, "search status: %s\n", search_status);
    fprintf(txt, "complete list: %s\n", complete_list ? "True" : "False");
    fprintf(txt, "total presentations found: %zu\n\n", result_count);
    for (size_t i = 0; i < result_count; i++) {
        char fbuf[4096], cbuf[4096];
        format_polynomial(results[i].f, config.key_len, fbuf, sizeof(fbuf));
        fprintf(txt, "[%d]\n", results[i].index);
        fprintf(txt, "f(x) = %s\n", fbuf);
        fprintf(txt, "middle_coefficient_a_%d = %lld\n", config.g, results[i].middle_coefficient);
        fprintf(txt, "canonical_presentation_index = %d\n", results[i].canonical_index);
        fprintf(txt, "\n");
    }
    fclose(txt);

    FILE *json = fopen(config.output_json, "w");
    if (!json) {
        perror(config.output_json);
        return;
    }
    fprintf(json, "{\n");
    fprintf(json, "  \"p\": %d,\n  \"g\": %d,\n  \"field\": \"F_%d\",\n", config.p, config.g, config.p);
    fprintf(json, "  \"presentation_reduction\": \"%s\",\n", config.reduction);
    fprintf(json, "  \"monic_enforced\": %s,\n", config.allow_nonmonic ? "false" : "true");
    fprintf(json, "  \"allow_nonmonic\": %s,\n", config.allow_nonmonic ? "true" : "false");
    fprintf(json, "  \"hasse_witt_prefilter\": %s,\n", config.use_hasse_witt ? "true" : "false");
    fprintf(json, "  \"reduction_class_count\": %d,\n", classes);
    fprintf(json, "  \"isomorphism_class_count\": %d,\n", classes);
    fprintf(json, "  \"search_status\": \"%s\",\n", search_status);
    fprintf(json, "  \"complete_list\": %s,\n", complete_list ? "true" : "false");
    fprintf(json, "  \"total_presentations_found\": %zu,\n", result_count);
    fprintf(json, "  \"curves\": [\n");
    for (size_t i = 0; i < result_count; i++) {
        char fbuf[4096], cbuf[4096];
        format_polynomial(results[i].f, config.key_len, fbuf, sizeof(fbuf));
        format_polynomial(results[i].canonical_f, config.key_len, cbuf, sizeof(cbuf));
        fprintf(json, "    {\n");
        fprintf(json, "      \"index\": %d,\n", results[i].index);
        fprintf(json, "      \"f_coeffs\": ");
        write_key_array(json, results[i].f, config.key_len);
        fprintf(json, ",\n      \"f_polynomial\": \"%s\",\n", fbuf);
        fprintf(json, "      \"middle_coefficient\": %lld,\n", results[i].middle_coefficient);
        fprintf(json, "      \"canonical_presentation_index\": %d,\n", results[i].canonical_index);
        fprintf(json, "      \"canonical_f_coeffs\": ");
        write_key_array(json, results[i].canonical_f, config.key_len);
        fprintf(json, ",\n      \"canonical_f_polynomial\": \"%s\",\n", cbuf);
        fprintf(json, "      \"reused_l_polynomial\": %s\n", results[i].reused_l_polynomial ? "true" : "false");
        fprintf(json, "    }%s\n", i + 1 == result_count ? "" : ",");
    }
    fprintf(json, "  ]\n}\n");
    fclose(json);
}

static void emit_result(const CurveResult *r) {
    if (config.quiet) return;
    char fbuf[4096];
    format_polynomial(r->f, config.key_len, fbuf, sizeof(fbuf));
    printf("Accepted [%d]\n", r->index);
    printf("f(x) = %s\n", fbuf);
    printf("middle_coefficient_a_%d = %lld\n", config.g, r->middle_coefficient);
    printf("canonical_presentation_index = %d\n", r->canonical_index);
    printf("presentations so far = %zu\n", result_count);
    printf("isomorphism classes so far = %d\n", reduction_class_count());
    fflush(stdout);
}

static bool polynomial_matrix_determinant(Poly *out, Poly matrix[MAX_GENUS][MAX_GENUS], int n, ulong p) {
    Poly work[MAX_GENUS][MAX_GENUS];
    for (int i = 0; i < n; i++) for (int j = 0; j < n; j++) work[i][j] = matrix[i][j];
    Poly previous_pivot;
    poly_one(&previous_pivot, p);
    int sign = 1;

    for (int k = 0; k < n - 1; k++) {
        int pr = -1, pc = -1;
        for (int r = k; r < n && pr < 0; r++) {
            for (int c = k; c < n; c++) {
                if (work[r][c].len != 0) {
                    pr = r;
                    pc = c;
                    break;
                }
            }
        }
        if (pr < 0) {
            poly_zero(out);
            return true;
        }
        if (pr != k) {
            for (int c = 0; c < n; c++) {
                Poly tmp = work[k][c];
                work[k][c] = work[pr][c];
                work[pr][c] = tmp;
            }
            sign = -sign;
        }
        if (pc != k) {
            for (int r = 0; r < n; r++) {
                Poly tmp = work[r][k];
                work[r][k] = work[r][pc];
                work[r][pc] = tmp;
            }
            sign = -sign;
        }
        Poly pivot = work[k][k];
        for (int r = k + 1; r < n; r++) {
            for (int c = k + 1; c < n; c++) {
                Poly a, b, numerator;
                poly_mul_to(&a, &work[r][c], &pivot, p);
                poly_mul_to(&b, &work[r][k], &work[k][c], p);
                poly_sub_to(&numerator, &a, &b, p);
                if (k > 0) poly_exact_div_to(&work[r][c], &numerator, &previous_pivot, p);
                else work[r][c] = numerator;
            }
        }
        previous_pivot = pivot;
        for (int r = k + 1; r < n; r++) poly_zero(&work[r][k]);
        for (int c = k + 1; c < n; c++) poly_zero(&work[k][c]);
    }
    *out = work[n - 1][n - 1];
    if (sign < 0) {
        for (int i = 0; i < out->len; i++) out->c[i] = out->c[i] ? p - out->c[i] : 0;
    }
    poly_normalize(out);
    return true;
}

static bool passes_hasse_witt_prefilter(const ulong *f) {
    int p = config.p, g = config.g;
    Poly fp, powered;
    poly_from_coeffs(&fp, f, config.key_len, (ulong)p);
    poly_pow_plain_to(&powered, &fp, (p - 1) / 2, (ulong)p);

    Poly mat[MAX_GENUS][MAX_GENUS];
    for (int r = 0; r < g; r++) {
        for (int c = 0; c < g; c++) {
            poly_zero(&mat[r][c]);
            int coeff_index = p * (r + 1) - (c + 1);
            ulong entry = (coeff_index >= 0 && coeff_index < powered.len) ? powered.c[coeff_index] % (ulong)p : 0;
            if (r == c) {
                mat[r][c].len = 1;
                mat[r][c].c[0] = 1;
                if (entry) {
                    mat[r][c].len = 2;
                    mat[r][c].c[1] = entry ? (ulong)p - entry : 0;
                }
            } else if (entry) {
                mat[r][c].len = 2;
                mat[r][c].c[0] = 0;
                mat[r][c].c[1] = (ulong)p - entry;
            }
            poly_normalize(&mat[r][c]);
        }
    }
    Poly det;
    polynomial_matrix_determinant(&det, mat, g, (ulong)p);
    for (int k = 1; k < g; k++) {
        if (k < det.len && det.c[k] % (ulong)p != 0) return false;
    }
    return true;
}

static Field make_field(int p, int degree) {
    Field field;
    field.p = p;
    field.degree = degree;
    field.order = ipow_ll(p, degree);
    poly_zero(&field.modulus);
    if (degree == 1) {
        field.modulus.len = 2;
        field.modulus.c[0] = 0;
        field.modulus.c[1] = 1;
        return field;
    }
    nmod_poly_t modulus;
    nmod_poly_init(modulus, (ulong)p);
    nmod_poly_minimal_irreducible(modulus, degree);
    poly_zero(&field.modulus);
    field.modulus.len = degree + 1;
    for (int i = 0; i <= degree; i++) field.modulus.c[i] = nmod_poly_get_coeff_ui(modulus, i) % (ulong)p;
    poly_normalize(&field.modulus);
    nmod_poly_clear(modulus);
    return field;
}

static void field_coeffs(const Field *field, long long value, Poly *out) {
    poly_zero(out);
    out->len = field->degree;
    for (int i = 0; i < field->degree; i++) {
        out->c[i] = (ulong)(value % field->p);
        value /= field->p;
    }
    poly_normalize(out);
}

static long long field_element(const Field *field, const Poly *poly) {
    long long value = 0;
    long long multiplier = 1;
    for (int i = 0; i < field->degree; i++) {
        ulong coeff = i < poly->len ? poly->c[i] % (ulong)field->p : 0;
        value += (long long)coeff * multiplier;
        multiplier *= field->p;
    }
    return value;
}

static long long field_add(const Field *field, long long a, long long b) {
    Poly pa, pb, sum;
    field_coeffs(field, a, &pa);
    field_coeffs(field, b, &pb);
    poly_add_to(&sum, &pa, &pb, (ulong)field->p);
    return field_element(field, &sum);
}

static long long field_mul(const Field *field, long long a, long long b) {
    if (field->degree == 1) return (a * b) % field->p;
    Poly pa, pb, prod, reduced;
    field_coeffs(field, a, &pa);
    field_coeffs(field, b, &pb);
    poly_mul_to(&prod, &pa, &pb, (ulong)field->p);
    poly_mod_to(&reduced, &prod, &field->modulus, (ulong)field->p);
    return field_element(field, &reduced);
}

static long long eval_poly_field(const Field *field, const ulong *f) {
    (void)field;
    (void)f;
    return 0;
}

static long long eval_at_field_element(const Field *field, const ulong *f, long long x) {
    long long result = 0;
    for (int i = config.key_len - 1; i >= 0; i--) {
        result = field_add(field, field_mul(field, result, x), (long long)(f[i] % (ulong)field->p));
    }
    return result;
}

static long long count_points(int k, const ulong *f) {
    Field field = make_field(config.p, k);
    int *square_counts = calloc((size_t)field.order, sizeof(int));
    if (!square_counts) {
        perror("calloc");
        exit(1);
    }
    for (long long y = 0; y < field.order; y++) {
        long long square = field_mul(&field, y, y);
        square_counts[square]++;
    }
    long long finite = 0;
    for (long long x = 0; x < field.order; x++) {
        long long value = eval_at_field_element(&field, f, x);
        finite += square_counts[value];
    }
    free(square_counts);
    int degree = key_degree(f, config.key_len);
    ulong leading = f[degree] % (ulong)config.p;
    long long infinity = 1;
    if (degree % 2 == 0) {
        infinity = 0;
        for (long long y = 0; y < field.order; y++) {
            if (field_mul(&field, y, y) == (long long)leading) infinity++;
        }
    }
    return finite + infinity;
}

static long long next_l_coefficient(const long long *coeffs, const long long *power_sums, int k) {
    long long numerator = power_sums[k - 1];
    for (int i = 1; i < k; i++) numerator += coeffs[i] * power_sums[k - i - 1];
    if (numerator % k != 0) {
        fprintf(stderr, "Newton identity produced nonintegral coefficient at k=%d\n", k);
        exit(1);
    }
    return -numerator / k;
}

static int trinomial_middle_coefficient(const ulong *f, long long *middle) {
    long long coeffs[MAX_GENUS + 1] = {0};
    long long power_sums[MAX_GENUS + 1] = {0};
    coeffs[0] = 1;
    for (int k = 1; k <= config.g; k++) {
        power_sums[k - 1] = ipow_ll(config.p, k) + 1 - count_points(k, f);
        coeffs[k] = next_l_coefficient(coeffs, power_sums, k);
        if (k < config.g && coeffs[k] != 0) return 0;
    }
    *middle = coeffs[config.g];
    return 1;
}

static bool affine_transform(const ulong *poly, int scale, int shift, ulong *out) {
    Poly result, linear, tmp, plus;
    poly_zero(&result);
    poly_zero(&linear);
    linear.len = 2;
    linear.c[0] = (ulong)shift % (ulong)config.p;
    linear.c[1] = (ulong)scale % (ulong)config.p;
    for (int i = config.key_len - 1; i >= 0; i--) {
        poly_mul_to(&tmp, &result, &linear, (ulong)config.p);
        plus = tmp;
        if (poly[i]) {
            if (plus.len < 1) plus.len = 1;
            plus.c[0] = (plus.c[0] + poly[i]) % (ulong)config.p;
        }
        poly_normalize(&plus);
        result = plus;
    }
    memset(out, 0, sizeof(ulong) * (size_t)config.key_len);
    for (int i = 0; i < result.len && i < config.key_len; i++) out[i] = result.c[i];
    int degree = key_degree(out, config.key_len);
    int original_degree = key_degree(poly, config.key_len);
    if (degree != original_degree) return false;
    return normalize_leading_square_class(out, config.key_len, (ulong)config.p, config.allow_nonmonic);
}

static void key_poly_mul_accum(ulong *acc, const Poly *a, const Poly *b, ulong coeff) {
    Poly term;
    poly_mul_to(&term, a, b, (ulong)config.p);
    for (int i = 0; i < term.len && i < config.key_len; i++) {
        acc[i] = (acc[i] + mod_mul(coeff, term.c[i], (ulong)config.p)) % (ulong)config.p;
    }
}

static bool pgl2_transform(const ulong *poly, int a, int b, int c, int d, ulong *out) {
    int binary_degree = 2 * config.g + 2;
    Poly *first = calloc((size_t)binary_degree + 1, sizeof(Poly));
    Poly *second = calloc((size_t)binary_degree + 1, sizeof(Poly));
    Poly l1, l2;
    if (!first || !second) {
        perror("calloc");
        exit(1);
    }
    poly_one(&first[0], (ulong)config.p);
    poly_one(&second[0], (ulong)config.p);
    poly_zero(&l1);
    poly_zero(&l2);
    l1.len = 2;
    l1.c[0] = (ulong)b % (ulong)config.p;
    l1.c[1] = (ulong)a % (ulong)config.p;
    l2.len = 2;
    l2.c[0] = (ulong)d % (ulong)config.p;
    l2.c[1] = (ulong)c % (ulong)config.p;
    for (int i = 1; i <= binary_degree; i++) {
        poly_mul_to(&first[i], &first[i - 1], &l1, (ulong)config.p);
        poly_mul_to(&second[i], &second[i - 1], &l2, (ulong)config.p);
    }
    memset(out, 0, sizeof(ulong) * (size_t)config.key_len);
    for (int i = 0; i < config.key_len && i <= binary_degree; i++) {
        if (poly[i]) key_poly_mul_accum(out, &first[i], &second[binary_degree - i], poly[i]);
    }
    int degree = key_degree(out, config.key_len);
    bool ok = degree == 2 * config.g + 1 || degree == 2 * config.g + 2;
    if (ok) ok = normalize_leading_square_class(out, config.key_len, (ulong)config.p, config.allow_nonmonic);
    free(first);
    free(second);
    return ok;
}

static void insert_orbit_members(HashMap *map, const ulong *poly, int owner) {
    ulong transformed[MAX_DEGREE + 1];
    bool pgl2 = strcmp(config.reduction, "pgl2") == 0 || strcmp(config.reduction, "pgl2save") == 0;
    if (!pgl2) {
        for (int scale = 1; scale < config.p; scale++) {
            for (int shift = 0; shift < config.p; shift++) {
                if (affine_transform(poly, scale, shift, transformed)) map_setdefault(map, transformed, owner);
            }
        }
        return;
    }
    for (int a = 0; a < config.p; a++) {
        for (int b = 0; b < config.p; b++) {
            for (int c = 0; c < config.p; c++) {
                for (int d = 0; d < config.p; d++) {
                    int det = (a * d - b * c) % config.p;
                    if (det < 0) det += config.p;
                    if (det == 0) continue;
                    int entries[4] = {a, b, c, d};
                    int first = 0;
                    while (first < 4 && entries[first] == 0) first++;
                    ulong inv = mod_inv((ulong)entries[first], (ulong)config.p);
                    int na = (int)mod_mul((ulong)a, inv, (ulong)config.p);
                    int nb = (int)mod_mul((ulong)b, inv, (ulong)config.p);
                    int nc = (int)mod_mul((ulong)c, inv, (ulong)config.p);
                    int nd = (int)mod_mul((ulong)d, inv, (ulong)config.p);
                    if (pgl2_transform(poly, na, nb, nc, nd, transformed)) {
                        map_setdefault(map, transformed, owner);
                    }
                }
            }
        }
    }
}

static bool next_lower_coeffs(ulong *f, int degree) {
    for (int i = 0; i < degree; i++) {
        f[i]++;
        if (f[i] < (ulong)config.p) return true;
        f[i] = 0;
    }
    return false;
}

static void handle_signal(int sig) {
    (void)sig;
    interrupted = 1;
}

static void search(void) {
    Stats stats = {0};
    HashMap accepted, seen;
    map_init(&accepted, config.key_len);
    map_init(&seen, config.key_len);
    bool save_repeats = strcmp(config.reduction, "pgl2save") == 0 || strcmp(config.reduction, "affinesave") == 0;
    ulong leading_values[2] = {1, 0};
    int leading_count = 1;
    if (config.allow_nonmonic) {
        leading_values[1] = nonsquare_representative((ulong)config.p);
        leading_count = 2;
    }

    for (int degree_index = 0; degree_index < 2 && !interrupted; degree_index++) {
        int degree = 2 * config.g + 1 + degree_index;
        for (int li = 0; li < leading_count && !interrupted; li++) {
            ulong f[MAX_DEGREE + 1] = {0};
            f[degree] = leading_values[li];
            while (!interrupted) {
                if (poly_squarefree_key(f, config.key_len, (ulong)config.p)) {
                    stats.considered++;
                    if (config.verbose) {
                        char fbuf[4096];
                        format_polynomial(f, config.key_len, fbuf, sizeof(fbuf));
                        printf("[%lld]\nConsidering f(x) = %s\nChecking %s repeats.\n", stats.considered, fbuf, config.reduction);
                        fflush(stdout);
                    }
                    int accepted_owner = map_get(&accepted, f);
                    if (accepted_owner) {
                        if (save_repeats) {
                            CurveResult *owner = &results[accepted_owner - 1];
                            append_result(f, owner->middle_coefficient, owner->canonical_index, owner->canonical_f, true);
                            write_results();
                            emit_result(&results[result_count - 1]);
                            if (config.max_curves > 0 && (int)result_count >= config.max_curves) goto done;
                        } else {
                            stats.skipped_by_reduction++;
                        }
                    } else {
                        int seen_owner = map_get(&seen, f);
                        if (seen_owner) {
                            stats.skipped_by_reduction++;
                        } else {
                            insert_orbit_members(&seen, f, (int)stats.considered);
                            stats.checked++;
                            if (config.use_hasse_witt && !passes_hasse_witt_prefilter(f)) {
                                stats.rejected_by_hasse_witt++;
                            } else {
                                long long middle = 0;
                                if (trinomial_middle_coefficient(f, &middle)) {
                                    append_result(f, middle, (int)result_count + 1, f, false);
                                    insert_orbit_members(&accepted, f, (int)result_count);
                                    write_results();
                                    emit_result(&results[result_count - 1]);
                                    if (config.max_curves > 0 && (int)result_count >= config.max_curves) goto done;
                                } else {
                                    stats.rejected_by_early_l_coefficient++;
                                }
                            }
                        }
                    }
                }
                if (!next_lower_coeffs(f, degree)) break;
                f[degree] = leading_values[li];
            }
        }
    }
done:
    if (interrupted) {
        strcpy(search_status, "interrupted");
        complete_list = false;
    } else if (config.max_curves > 0 && (int)result_count >= config.max_curves) {
        strcpy(search_status, "max_reached");
        complete_list = false;
    } else {
        strcpy(search_status, "complete");
        complete_list = true;
    }
    write_results();
    printf("\nSUMMARY\n");
    printf("Considered %lld squarefree %s presentations.\n", stats.considered, config.allow_nonmonic ? "square-class-normalized" : "monic");
    printf("Checked %lld new reduction-orbit representatives after %s reduction.\n", stats.checked, config.reduction);
    printf("Skipped %lld %s-equivalent presentations.\n", stats.skipped_by_reduction, config.reduction);
    if (config.use_hasse_witt) printf("Hasse-Witt rejected %lld presentations.\n", stats.rejected_by_hasse_witt);
    printf("Early-rejected %lld presentations.\n", stats.rejected_by_early_l_coefficient);
}

static void usage(const char *argv0) {
    fprintf(stderr, "usage: %s p g [--max N] [--output FILE] [--reduction pgl2|pgl2save|affine|affinesave] [--monic-only] [--no-hasse-witt-prefilter] [--quiet|--verbose]\n", argv0);
}

int main(int argc, char **argv) {
    memset(&config, 0, sizeof(config));
    config.max_curves = 0;
    config.allow_nonmonic = true;
    config.use_hasse_witt = true;
    strcpy(config.reduction, "pgl2save");
    output_paths("hyperelliptic_results_c.txt", config.output_txt, sizeof(config.output_txt), config.output_json, sizeof(config.output_json));

    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }
    config.p = atoi(argv[1]);
    config.g = atoi(argv[2]);
    if (!is_prime_int(config.p) || config.p == 2) {
        fprintf(stderr, "Error: p must be an odd prime\n");
        return 1;
    }
    if (config.g < 1 || config.g > MAX_GENUS) {
        fprintf(stderr, "Error: g must be between 1 and %d\n", MAX_GENUS);
        return 1;
    }
    config.key_len = 2 * config.g + 3;
    if (config.key_len > MAX_DEGREE + 1) {
        fprintf(stderr, "Error: degree exceeds C limit %d\n", MAX_DEGREE);
        return 1;
    }

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--max") == 0 && i + 1 < argc) {
            config.max_curves = atoi(argv[++i]);
        } else if (strcmp(argv[i], "--output") == 0 && i + 1 < argc) {
            output_paths(argv[++i], config.output_txt, sizeof(config.output_txt), config.output_json, sizeof(config.output_json));
        } else if (strcmp(argv[i], "--reduction") == 0 && i + 1 < argc) {
            snprintf(config.reduction, sizeof(config.reduction), "%s", argv[++i]);
            if (strcmp(config.reduction, "pgl2") && strcmp(config.reduction, "pgl2save") &&
                strcmp(config.reduction, "affine") && strcmp(config.reduction, "affinesave")) {
                fprintf(stderr, "unknown reduction: %s\n", config.reduction);
                return 1;
            }
        } else if (strcmp(argv[i], "--monic-only") == 0) {
            config.allow_nonmonic = false;
        } else if (strcmp(argv[i], "--no-hasse-witt-prefilter") == 0) {
            config.use_hasse_witt = false;
        } else if (strcmp(argv[i], "--quiet") == 0) {
            config.quiet = true;
        } else if (strcmp(argv[i], "--verbose") == 0) {
            config.verbose = true;
        } else {
            usage(argv[0]);
            return 1;
        }
    }
    if (config.max_curves < 0) {
        fprintf(stderr, "Error: --max must be nonnegative\n");
        return 1;
    }

    signal(SIGINT, handle_signal);
    signal(SIGTERM, handle_signal);

    printf("Searching over F_%d.\n", config.p);
    printf("Using %s presentation reduction.\n", config.reduction);
    printf("Saving detailed results to %s and %s.\n", config.output_txt, config.output_json);
    fflush(stdout);
    write_results();
    search();
    printf("Done. Saved %zu presentations to %s and %s.\n", result_count, config.output_txt, config.output_json);
    return interrupted ? 130 : 0;
}
