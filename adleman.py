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
        # beta вычисляем, но далее не используем (можно использовать для проверки)
        beta = pow(b, (p-1) // q, p)

        # Поиск x_i по модулю q^e
        x_i = 0
        for k in range(e):
            exponent = (p - 1) // (q**(k+1))
            # Вычисляем lhs = (b * (a^(x_i))^(-1)) mod p
            lhs = (b * mod_inverse(pow(a, x_i, p), p)) % p
            term = pow(lhs, exponent, p)
            
            # Находим такое d_k, что gamma^(d_k) ≡ term (mod p)
            d_k = 0
            while pow(gamma, d_k, p) != term:
                d_k += 1
            x_i += d_k * (q**k)
        x_mods.append((x_i, m))

    # Шаг 4: Применяем Китайскую теорему об остатках
    # Сначала вычисляем полный модуль (N)
    N = 1
    for (_, mod) in x_mods:
        N *= mod
    x = 0
    for residue, mod in x_mods:
        N_i = N // mod
        inv = mod_inverse(N_i, mod)
        x += residue * N_i * inv
    return x % N

# Пример использования
if __name__ == "__main__":
    a = 2
    b = 9
    p = 11  # Простое число, p-1 = 10 = 2 * 5 (гладкое)
    x = adleman_discrete_log(a, b, p)
    print(f"Решение: {x} (проверка: {a}^{x} mod {p} = {pow(a, x, p)} = {b})")
