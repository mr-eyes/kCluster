import kProcessor as kp

def shouldbe(kmer):
    return kmer.replace("A","00").replace("C","01").replace("G","11").replace("T","10")

def str_to_int(kmer):
    _map = {"A":0, "C":1, "T":2, "G":3}
    strint = 0
    for n in kmer:
        curr = _map[n]
        strint = strint | curr
        strint = strint << 2
        
    return strint >> 2


with open("idx_min_test.dump" , 'r') as dump:
    next(dump)
    for line in dump:
        line = line.strip().split("\t")
        kmer = line[0]
        hashed = line[-1]
        if hashed != shouldbe(kmer):
            print(kmer)
            break
        else:
            print("OK")