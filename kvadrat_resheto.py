#!/usr/bin/python3
from math import sqrt, log2, ceil, floor, isqrt
import random
from math import gcd
import sympy

# Параметры
SIQS_TRIAL_DIVISION_EPS = 25
SIQS_MIN_PRIME_POLYNOMIAL = 400
SIQS_MAX_PRIME_POLYNOMIAL = 4000

# Классы факторной базы и полинома
class FactorBasePrime:
    def __init__(self, p, tmem, lp):
        self.p = p; self.tmem = tmem; self.lp = round(log2(p))
        self.soln1 = self.soln2 = self.ainv = None

class Polynomial:
    def __init__(self, coeff, a, b):
        self.coeff = coeff; self.a = a; self.b = b
    def eval(self, x):
        res = 0
        for c in reversed(self.coeff):
            res = res*x + c
        return res

# Используем sympy для всех криптоарифметических операций
sqrt_mod_prime = lambda a, p: sympy.sqrt_mod(a, p, all_roots=False)
modinv = lambda a, m: pow(a, -1, m)
is_quad = lambda a, p: sympy.legendre_symbol(a, p) == 1

# Вместо lowest_set_bit: используем (a & -a).bit_length()-1
lowest_set_bit = lambda a: (a & -a).bit_length()-1

# Вычисление квадратных корней для найденных зависимостей
def siqs_calc_sqrts(indices, rel):
    X = 1; Y = 1
    for i in indices:
        X *= rel[i][0]
        Y *= rel[i][1]
    return [X, isqrt(Y)]
def siqs_factor_from_square(n, indices, rel):
    s1, s2 = siqs_calc_sqrts(indices, rel)
    return gcd(abs(s1-s2), n)
def siqs_find_factors(n, squares, rel):
    fac = []; rem = n
    for inds in squares:
        f = siqs_factor_from_square(n, inds, rel)
        if 1 < f < rem:
            while rem % f == 0:
                fac.append(f); rem //= f
    if rem > 1: fac.append(rem)
    return fac

# Формирование факторной базы с использованием списка маленьких простых
def siqs_factor_base_primes(n, nf, small_primes):
    fb = []
    for p in small_primes:
        if p == 2:  # пропускаем 2
            continue
        if is_quad(n, p):
            fb.append(FactorBasePrime(p, sqrt_mod_prime(n % p, p), log2(p)))
            if len(fb) >= nf:
                break
    return fb


# Генерация первого полинома
def siqs_find_first_poly(n, m, fb):
    p_min = p_max = None
    for i, f in enumerate(fb):
        if p_min is None and f.p >= SIQS_MIN_PRIME_POLYNOMIAL: p_min = i
        if p_max is None and f.p > SIQS_MAX_PRIME_POLYNOMIAL: p_max = i-1; break
    if p_max is None: p_max = len(fb)-1
    if p_min is None or (p_max-p_min < 20): p_min = min(p_min if p_min is not None else 0, 5)
    target = sqrt(2*n)/m
    target1 = target/(((fb[p_min].p+fb[p_max].p)/2)**0.5)
    best_a = best_ratio = None; best_q = None
    for _ in range(30):
        a = 1; q = []
        while a < target1:
            p_i = random.randint(p_min, p_max)
            if p_i in q: continue
            a *= fb[p_i].p; q.append(p_i)
        ratio = a/target
        if best_ratio is None or (ratio>=0.9 and ratio<best_ratio) or (best_ratio<0.9 and ratio>best_ratio):
            best_ratio = ratio; best_a = a; best_q = q
    a, q = best_a, best_q; s = len(q)
    B = []
    for l in range(s):
        f = fb[q[l]]
        gamma = (f.tmem * modinv(a//f.p, f.p)) % f.p
        if gamma > f.p//2: gamma = f.p - gamma
        B.append(a//f.p * gamma)
    b = sum(B) % a; b_orig = b
    if 2*b > a: b = a - b
    assert 0 < b and 2*b <= a and ((b*b - n) % a == 0)
    g = Polynomial([b*b - n, 2*a*b, a*a], a, b_orig)
    h = Polynomial([b, a], a, b_orig)
    for f in fb:
        if a % f.p:
            f.ainv = modinv(a, f.p)
            f.soln1 = (f.ainv * (f.tmem - b)) % f.p
            f.soln2 = (f.ainv * (-f.tmem - b)) % f.p
    return g, h, B

def siqs_find_next_poly(n, fb, i, g, B):
    v = lowest_set_bit(i)+1
    z = -1 if ceil(i/(2**v)) % 2==1 else 1
    b = (g.b + 2*z*B[v-1]) % g.a; a = g.a; b_orig = b
    if 2*b > a: b = a - b
    assert ((b*b - n) % a == 0)
    g = Polynomial([b*b - n, 2*a*b, a*a], a, b_orig)
    h = Polynomial([b, a], a, b_orig)
    for f in fb:
        if a % f.p:
            f.soln1 = (f.ainv * (f.tmem - b)) % f.p
            f.soln2 = (f.ainv * (-f.tmem - b)) % f.p
    return g, h

# Просеивание
def siqs_sieve(fb, m):
    sieve = [0]*(2*m+1)
    for f in fb:
        if f.soln1 is None: continue
        p = f.p
        i1 = -((m+f.soln1)//p)
        start1 = f.soln1 + i1*p + m
        for a in range(start1, 2*m+1, p): sieve[a] += f.lp
        i2 = -((m+f.soln2)//p)
        start2 = f.soln2 + i2*p + m
        for a in range(start2, 2*m+1, p): sieve[a] += f.lp
    return sieve

def siqs_trial_divide(a, fb):
    div_idx = []
    for i, f in enumerate(fb):
        if a % f.p == 0:
            exp = 0
            while a % f.p == 0:
                a //= f.p; exp += 1
            div_idx.append((i, exp))
        if a==1: return div_idx
    return None

def siqs_trial_division(n, sieve, fb, smooth, g, h, m, req):
    limit = log2(m*sqrt(n)) - SIQS_TRIAL_DIVISION_EPS
    for i, sa in enumerate(sieve):
        if sa >= limit:
            x = i - m; gx = g.eval(x)
            d = siqs_trial_divide(gx, fb)
            if d is not None:
                u = h.eval(x); v = gx
                smooth.append((u, v, d))
                if len(smooth) >= req: return True
    return False

# Линейная алгебра (минимальный вариант – без оптимизаций)
def siqs_build_matrix(fb, smooth):
    M = []
    for u, v, div in smooth:
        row = [0]*len(fb)
        for j, exp in div: row[j] = exp % 2
        M.append(row)
    return M

def siqs_build_matrix_opt(M):
    m = len(M[0])
    cols = ["" for _ in range(m)]
    for row in M:
        for j, bit in enumerate(row):
            cols[j] += "1" if bit else "0"
    return [int(c[::-1],2) for c in cols], len(M), m

def add_column_opt(M_opt, tgt, src):
    M_opt[tgt] ^= M_opt[src]
def find_pivot_column_opt(M_opt, j):
    return None if M_opt[j]==0 else lowest_set_bit(M_opt[j])
def siqs_solve_matrix_opt(M_opt, n, m):
    row_mark = [False]*n; pivots = [-1]*m
    for j in range(m):
        i = find_pivot_column_opt(M_opt, j)
        if i is not None:
            pivots[j] = i; row_mark[i] = True
            for k in range(m):
                if k != j and ((M_opt[k]>>i)&1):
                    add_column_opt(M_opt, k, j)
    squares = []
    for i in range(n):
        if not row_mark[i]:
            inds = [i]
            for j in range(m):
                if (M_opt[j]>>i)&1: inds.append(pivots[j])
            squares.append(inds)
    return squares

def siqs_choose_nf_m(d):
    if d<=34: return 200,65536
    if d<=36: return 300,65536
    if d<=38: return 400,65536
    if d<=40: return 500,65536
    if d<=42: return 600,65536
    if d<=44: return 700,65536
    if d<=48: return 1000,65536
    if d<=52: return 1200,65536
    if d<=56: return 2000,65536*3
    if d<=60: return 4000,65536*3
    if d<=66: return 6000,65536*3
    if d<=74: return 10000,65536*3
    if d<=80: return 30000,65536*3
    if d<=88: return 50000,65536*3
    if d<=94: return 60000,65536*9
    return 100000,65536*9

# Инициализация маленьких простых через решето
def trial_div_init_primes(n, upper):
    global small_primes
    is_prime = [True]*(upper+1)
    is_prime[0:2] = [False, False]
    fac = []; small_primes = []
    rem = n; maxi = isqrt(upper)
    for i in range(2, maxi+1):
        if is_prime[i]:
            small_primes.append(i)
            while rem % i == 0:
                rem //= i; fac.append(i)
                if sympy.isprime(rem): 
                    fac.append(rem); return fac, 1
            for j in range(i*i, upper+1, i): is_prime[j] = False
    for i in range(maxi+1, upper+1):
        if is_prime[i]:
            small_primes.append(i)
            while rem % i == 0:
                rem //= i; fac.append(i)
                if sympy.isprime(rem): 
                    fac.append(rem); return fac, 1
    return fac, rem

def siqs_factorise(n):
    d = len(str(n))
    nf, m = siqs_choose_nf_m(d)
    fb = siqs_factor_base_primes(n, nf, small_primes)
    req_ratio = 1.05; smooth = []; i_poly = 0; prev = 0
    while True:
        req = round(len(fb)*req_ratio)
        while True:
            if i_poly == 0:
                g, h, B = siqs_find_first_poly(n, m, fb)
            else:
                g, h = siqs_find_next_poly(n, fb, i_poly, g, B)
            i_poly += 1
            if i_poly >= 2 ** (len(B) - 1):
                i_poly = 0
            sieve = siqs_sieve(fb, m)
            if siqs_trial_division(n, sieve, fb, smooth, g, h, m, req):
                break
        M = siqs_build_matrix(fb, smooth)
        M_opt, M_n, M_m = siqs_build_matrix_opt(M)
        squares = siqs_solve_matrix_opt(M_opt, M_n, M_m)
        fac = siqs_find_factors(n, squares, smooth)
        if len(fac)>1: return fac
        req_ratio += 0.05

if __name__=='__main__':
    n = 904114958460043420393193689861
    print(n)
    _, _ = trial_div_init_primes(n, 1000000)
    print(siqs_factorise(n))
