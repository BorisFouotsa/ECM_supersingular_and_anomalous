import random
import time
#random.seed(time.time())

load('anomalous_curves.sage')

DISC_SUPERSINGULAR = [3,4,7,8,11,19,43,67,163]
SUPERSINGULAR_J = [54000, 287496,-3375, 8000, -32768, -884736,  -884736000,-147197952000, -262537412640768000]

DISC_ANOMALOUS = [3,11,19,27,43,67,163]
ANOMALOUS_J = [0, -32768, -884736, -12288000, -884736000,-147197952000, -262537412640768000]

CURVES_ANOMALOUS = [
    CURVES_ANOMALOUS_3,
    CURVES_ANOMALOUS_11,
    CURVES_ANOMALOUS_19,
    CURVES_ANOMALOUS_27,
    CURVES_ANOMALOUS_43,
    CURVES_ANOMALOUS_67,
    CURVES_ANOMALOUS_163
]

def logbasetwo(n):

    m =1;
    while 2^m <= n:
        m = m+1;
    return m;

def extract_first_integer_in_text(txt):
    for s in txt.split():
        if s.isdigit():
            return(int(s))
    return 0

def gen_modulus_anomalous(m, D = 0, secret = False):

    if D not in DISC_ANOMALOUS:
        D = random.choice(DISC_ANOMALOUS)

    m = m+randint(1,m)

    p = D*m*(m+1)+Integer((1+D)/4)

    while is_prime(p) == False:
        m = m+1
        p = D*m*(m+1)+Integer((1+D)/4)

    n = logbasetwo(p)
    q = randint(2^(n-1),2^n)
    q = next_prime(q)

    if secret == True:
        return p,q, D
    return p*q, D



def gen_modulus_supersingular(B, sizep, secret=False):

    p_approx = 2^sizep
    p = 2*random_prime(B)
    while p < p_approx:
        p = p*random_prime(B)
    i = 5
    while is_prime(p*i-1) == False:
        i = next_prime(i+1)
    p = p*i-1


    n = logbasetwo(p)
    q = randint(2^(n-1),2^n)
    q = next_prime(q)

    bound = max(i,B)   # max(i, B) is  the smoothness bound of p+1

    if secret == True:
        return p,q, bound
    return p*q, bound



def gen_curve_j0(n):

    Zn = Integers(n)
    u = Zn(randint(1,n-1))
    v = Zn(randint(1,n-1))
    while v^2 == u^3:
        v = Zn(randint(1,n-1))

    B = Zn(v^2 - u^3)
    E = EllipticCurve([0, B])
    P = E(u,v)

    return P



def scalar_mult(P, r):

    E = P.curve()
    Q = P
    R = E(0)
    while r > 0:
        if Mod(r,2) == 1:
            R = R + Q
        Q = Q + Q
        r = floor(r/2)

    return R



def try_point_ECM_anomalous(P, n):

    try:
        Q = scalar_mult(P, n)
    except ZeroDivisionError as e:
        txt = getattr(e, 'message', str(e))
        m = extract_first_integer_in_text(txt)
        d = gcd(m,n)
        if d != 1:
            return d
    return 0



def ECM_supersingular(n, B, B1=1, Jinv=SUPERSINGULAR_J):
    '''B is the smoothness bound of p+1 and B1<B is a bound for which (B!)^{log B}
    is calculated beforre entering perfoming the scalar multiplications'''

    Zn = Integers(n)
    t = logbasetwo(B)
    k0 = factorial(B1)^t

    for j in Jinv:
        if j not in SUPERSINGULAR_J:
            continue  # 'Your j-invaraints must be in SUPERSINGULAR_J'
        a = Zn(27*j/(4*(1728 - j)))
        E = EllipticCurve([a,-a])
        Q = E(Zn(1),Zn(1))

        for k in range(B1,B+1):
            try:
                if k == B1:
                    Q = scalar_mult(Q, k0)
                Q = scalar_mult(Q, k^t)
            except ZeroDivisionError as e:
                txt = getattr(e, 'message', str(e))
                del(k0)
                m = extract_first_integer_in_text(txt)
                d = gcd(m,n)
                if d != 1:
                    return d
    del(k0)
    return 0



def ECM_anomalous(n, Disc=DISC_ANOMALOUS):

    Zn = Integers(n)

    for D in Disc:

        if D == 3:    # j = 0
            for t in range(60):
                P = gen_curve_j0(n)
                d = try_point_ECM_anomalous(P, n)
                if d != 0:
                    return d

        else:
            if D not in DISC_ANOMALOUS:
                return f'Your discriminants must be in {DISC_ANOMALOUS}'
            index = DISC_ANOMALOUS.index(D)
            curves = CURVES_ANOMALOUS[index]

            for [a_coefs,[u,v]] in curves:
                E = EllipticCurve(a_coefs)
                E = E.change_ring(Zn)
                P = E(Zn(u),Zn(v))
                d = try_point_ECM_anomalous(P, n)
                if d != 0:
                    return d
    return 0


def test_supersingular(bound=2^10, sizep=512):

    print(f'You are running a test on the ECM based p+1  factoring algorithm.\n')
    print(f'A {2*sizep} bits "vulnerable" RSA nodulus is being generated...\n')

    p,q,B = gen_modulus_supersingular(bound, sizep, secret = True)
    leg_symbols = [legendre_symbol(-D,p) for D in DISC_SUPERSINGULAR]
    while -1 not in leg_symbols:
        p,q,B = gen_modulus_supersingular(bound, sizep, secret = True)
        leg_symbols = [legendre_symbol(-D,p) for D in DISC_SUPERSINGULAR]

    n = p*q
    t = leg_symbols.index(-1)
    D = DISC_SUPERSINGULAR[t]
    j = SUPERSINGULAR_J[t]


    print(f'The RSA modulus is \n n = {n}, \n')
    print(f'n = p*q where p+1 is B = {B} smooth and (-D/p) = -1 with D = {D}.\n')
    print(f"Let's run the ECM p+1 on n with the curve having j-invaraint j = {j}.\n")

    d = ECM_supersingular(n, B, B1=1, Jinv = [j])

    if 1<d<n and Mod(n,d) == 0:
        print(f'A factor d of the RSA modulus  n was recovered,  \n d = {d}\n')
    else:
        print(f'the algorithm failed to factor the RSA modulus  n.')
    return 0



def test_anomalous(sizep=1024, D=0):

    print(f'You are running a test on the ECM with anomalous curve factoring algorithm.\n')
    print(f'A {2*sizep} bits "vulnerable" RSA nodulus is being generated...\n')
    print(f'The modulus  generation may last tens of seconds, especially for large size n...\n')

    n, D = gen_modulus_anomalous(2^sizep, secret = False)
    t = DISC_ANOMALOUS.index(D)
    j = ANOMALOUS_J[t]

    print(f'The RSA modulus is \n \n n = {n}, \n')
    print(f'n = p*q where p = Dm(m+1)+ (D+1)/4 with D = {D}.\n')
    print(f"Let's run the ECM on n with anomalous curves  having j-invariant j = {j}.\n")
    t = time.time()
    d = ECM_anomalous(n, Disc = [D])
    t1 = time.time() - t

    if 1<d<n and Mod(n,d) == 0:
        print(f'A factor d of the RSA modulus  n was recovered in {t1} seconds,  \n\n d = {d}\n')
    else:
        print(f'the algorithm failed to factor the RSA modulus  n. This means that none of the curves used was anomalous. Try again ;)')
    return 0


print('\n')
test_supersingular()
print('\n')
print('\n')
test_anomalous()
print('\n')
