import random
import time
random.seed(time.time())


def logbasetwo(n):

    m =1;
    while 2^m <= n:
        m = m+1;
    return m;


def gen_modulus_anomalous(D,m, secret = False):

    if D not in [3, 11, 19, 27, 43, 67, 163]:
        return 'D must be in [3,11,19,27,43,67,163]'

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


def gen_modulus_supersingular(B, sizep, resMod12 = False, secret = False):

    if resMod12 not in [False, 3,2,1]:
        return 'resMod12 must be in [False, 3,2,1]'

    p = 2*random_prime(B)
    while logbasetwo(p) < sizep:
        p = p*random_prime(B)
    i = 7
    while is_prime(p*i-1) == False:
        i = next_prime(i+1)
    p = p*i-1


    if resMod12 == 3:
        while Mod(p,4) != 3:
            gen_modulus_supersingular(B, sizep, 3, secret)

    if resMod12 == 2:
        while Mod(p,3) != 2:
            gen_modulus_supersingular(B, sizep, 2, secret)

    if resMod12 == 1:
        while Mod(p,12) != 1:
            gen_modulus_supersingular(B, sizep, 1, secret)

    n = logbasetwo(p)
    q = randint(2^(n-1),2^n)
    q = next_prime(q)

    if secret == True:
        return p,q
    return p*q



def gen_modulus_supersingular_no_constraint(B, sizep, secret = False):

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

    if secret == True:
        return p,q, max(i,B)
    return p*q, max(i,B)



def gen_curve_point_from_j(n,j):

    Zn = Integers(n)
    u = randint(1,n-1)
    v = u*randint(1,n-1) - u -1
    if j == 0:
        A = 0
        B = v^2 - u^3
    elif j == 1728:
        A = (v^2 - u^3)/u
        B = 0
    E = EllipticCurve([A, B])
    E = E.change_ring(Zn)
    P = E(Zn(u),Zn(v))

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



def ECM_supersingular(n, B1, B, Disc = DISC_SUPERSINGULAR):

    Zn = Integers(n)
    t = logbasetwo(B)
    k0 = factorial(B1)^t

    for D in Disc:
        index = DISC_SUPERSINGULAR.index(D)
        [a_coefs,[u,v]] = CURVES_SUPERSINGULAR[index]
        E = EllipticCurve(a_coefs)
        E = E.change_ring(Zn)
        Q = E(Zn(u),Zn(v))

        for k in range(B1,B+1):
            try:
                if k == B1+1:
                    Q = scalar_mult(Q, k0)
                Q = scalar_mult(Q, k^t)
            except ZeroDivisionError as e:
                ans = getattr(e, 'message', str(e)) # r =getattr(e, 'message', repr(e))
                del(k0)
                return ans, index
    return 0



def ECM_anomalous(n, Disc = DISC_ANOMALOUS):

    Zn = Integers(n)

    for D in Disc:

        if D == 3:    # j = 0
            for t in range(50):
                P = gen_curve_point_from_j(n,0)
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
