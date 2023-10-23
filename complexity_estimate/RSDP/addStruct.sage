def get_E(q,z):
    assert is_prime(q),"q is not prime!"
    assert (q-1)%z == 0, "z does not divide q-1!"

    F = GF(q)
    alpha = F.primitive_element()
    m = (q-1)/z
    g = alpha**m
    return [g^l for l in range(z)]

def shift_E(E, shift):
    return [e-shift for e in E if (e-shift) != 0]


def get_Eminus(E):
    E_temp = [-e for e in E] 
    E_minus = []
    for a in E_temp:
        if a not in E and a!=0:
            E_minus += [a]

    if len(E_minus) < len(E):
        print("size of -E is smaller than of E!")
    assert len(E_minus) == len(E),"size of -E is smaller than of E!"

    return E_minus
    
def get_Diff(E, onlyMax = False):
    """
    generates difference set {a-b|a,b in E}
    """
    Diff = []
    F = E[0].parent()
    F_list = F.list()
    q = F.cardinality()
    counts = [0 for _ in range(q)]
    for a in E:
        for b in E:
            x = a-b
            Diff += [x]
            counts[F_list.index(x)] += 1

    Diff = list( set(Diff) ) #unique
    Eminus = [-e for e in E]
    for x in F_list:
        if x in E or x in Eminus or x==0:
            counts[F_list.index(x)] = -1
    # only take those that are neither in E or in -E
    D = []
    for d in Diff:
        if d not in E and d not in Eminus and d!=0:
            if onlyMax:
                if counts[F_list.index(d)] == max(counts):
                    D += [d]
            else:
                D += [d]

    return D



def get_nd(E, verb, avg):
    """
    a in E
    b in E
    a+b in E
    """
    z = len(E)
    counts = [0 for _ in range(z)]
    for a in E:
        for b in E:
            x = a+b 
            if x in E:
                counts[E.index(x)] += 1

    irregular = 0
    for l in range(1,z):
        if counts[l] != counts[0]:
            irregular = 1
    
    if irregular:
            if verb:
                print("%Warning, D-repr not regular!")
                print("%",E)
                print("%",counts)
            if avg: return 1.0*sum(counts)/len(counts)
            else:   return max(counts)

    return counts[0]

def get_nb(E, D, verb, avg):
    """
    a in E
    b in D
    a+b in E
    """
    z = len(E)
    counts = [0 for _ in range(z)]
    for a in E:
        for b in D:
            x = a+b 
            if x in E:
                counts[E.index(x)] += 1

    irregular = 0
    for l in range(1,z):
        if counts[l] != counts[0]:
            irregular = 1
    
    if irregular:
            if verb:
                print("%Warning, B-repr not regular!")
                print("%",E)
                print("%",counts)
            if avg: return 1.0*sum(counts)/len(counts)
            else:   return max(counts)

    return counts[0]

def get_nc(E, Eminus, D, verb, avg):
    """
    a in E
    b in E
    a+b in D
    """
    z = len(D)
    counts = [0 for _ in range(z)]
    for a in E:
        for b in Eminus:
            x = a+b 
            if x in D:
                counts[D.index(x)] += 1
    irregular = 0
    for l in range(1,z):
        if counts[l] != counts[0]:
            irregular = 1
    
    if irregular:
            if verb:
                print("%Warning, C-repr not regular!")
                print("%",E)
                print("%",counts)
            if avg: return 1.0*sum(counts)/len(counts)
            else:   return max(counts)

    return counts[0]




