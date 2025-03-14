import math
from sympy import factorint, mod_inverse

def adleman_discrete_log(a, b, p):
    """Алгоритм Адлемана для дискретного логарифмирования."""
    if math.gcd(a, p) != 1 or math.gcd(b, p) != 1:
        raise ValueError("a и b должны быть взаимно просты с p")

    # Факторизация p-1
    factors = factorint(p - 1)
    prime_powers = [(q, e) for q, e in factors.items()]

    # Шаг 3: Решение для каждого q^e
    x_mods = []
    for q, e in prime_powers:
        m = q**e
        gamma = pow(a, (p-1) // q, p)
        beta = pow(b, (p-1) // q, p)

        # Поиск x_i mod q^e
        x_i = 0
        for k in range(e):
            exponent = (p-1) // q**(k+1)
            lhs = pow(b * mod_inverse(pow(a, x_i, p), 1, p))
            term = pow(lhs, exponent, p)
            
            # Поиск d_k такого, что gamma^(d_k) ≡ term mod p
            d_k = 0
            while pow(gamma, d_k, p) != term:
                d_k += 1
            
            x_i += d_k * q**k

        x_mods.append((x_i, m))

    # Шаг 4: Китайская теорема об остатках
    x = 0
    total_mod = 1
    for residue, mod in x_mods:
        total_mod *= mod
        term = residue * (total_mod // mod) * mod_inverse(total_mod // mod, mod)
        x += term
    return x % total_mod

# Пример использования
a = 2
b = 9
p = 11  # Простое число, p-1 = 10 = 2 * 5 (гладкое)
x = adleman_discrete_log(a, b, p)
print(f"Решение: {x} (проверка: {a}^{x} mod {p} = {pow(a, x, p)} = {b})")