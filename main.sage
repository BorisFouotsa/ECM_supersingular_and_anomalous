import random
import time
random.seed(time.time())


#load('supersingular_curves.sage')
load('anomalous_curves.sage')

SUPERSINGULAR_J = [5400, 287496,-3375, 8000, -32768, -884736,  -884736000,-147197952000, -262537412640768000]

def logbasetwo(n):

    m =1;
    while 2^m <= n:
        m = m+1;
    return m;


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
        return p,q
    return p*q



def gen_modulus_supersingular(B, sizep, secret = False):

    p = 2*random_prime(B)
    while logbasetwo(p) < sizep:
        p = p*random_prime(B)
    i = 5
    while is_prime(p*i-1) == False:
        i = next_prime(i+1)
    p = p*i-1


    n = logbasetwo(p)
    q = randint(2^(n-1),2^n)
    q = next_prime(q)

    bound = max(i,B)   # max(i, B) is  the smoothness bound of p+1
    #leg_symbols = [legendre_symbol(-D,p) for D in DISC_SUPERSINGULAR]  # this indicates which curves will be able to factor n

    if secret == True:
        return p,q, bound  #, leg_symbols
    return p*q, bound  #, leg_symbols



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
        s = getattr(e, 'message', str(e))
        return s
    return 0



def ECM_supersingular(n, B, B1=1):
    '''B is the smoothness bound of p+1 and B1<B is a bound for which (B!)^{log B}
    is calculated beforre entering perfoming the scalar multiplications'''

    Zn = Integers(n)
    t = logbasetwo(B)
    k0 = factorial(B1)^t

    for j in SUPERSINGULAR_J:
        a = Zn(27*j/(4*(1728 - j)))
        E = EllipticCurve([a,-a])
        Q = E(Zn(1),Zn(1))

        for k in range(B1,B+1):
            try:
                if k == B1:
                    Q = scalar_mult(Q, k0)
                Q = scalar_mult(Q, k^t)
            except ZeroDivisionError as e:
                ans = getattr(e, 'message', str(e))
                del(k0)
                return ans, index
    del(k0)
    return 0



def ECM_anomalous(n, Disc = DISC_ANOMALOUS):

    Zn = Integers(n)

    for D in Disc:

        if D == 3:    # j = 0
            for t in range(60):
                P = gen_curve_j0(n)
                s = try_point_ECM_anomalous(P, n)
                if s != 0:
                    return s

        else:
            if D not in DISC_ANOMALOUS:
                return 'Your discriminants must be in [3,11,19,27,43,67,163]'
            index = DISC_ANOMALOUS.index(D)
            curves = CURVES_ANOMALOUS[index]

            for [a_coefs,[u,v]] in curves:
                E = EllipticCurve(a_coefs)
                E = E.change_ring(Zn)
                P = E(Zn(u),Zn(v))
                s = try_point_ECM_anomalous(P, n)
                if s != 0:
                    return s
    return 0
