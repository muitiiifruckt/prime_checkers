from sympy import randprime

def generate_prime(bits):
    lower_bound = 2**(bits - 1)
    upper_bound = 2**bits - 1
    return randprime(lower_bound, upper_bound)


def generate_N(bits):
    p = generate_prime(bits//2)
    q = generate_prime(bits//2)
    N = p * q
    return N  
if __name__ == "__main__":
    print(generate_N(80))