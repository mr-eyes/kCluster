import kProcessor as kp
import random

# Mimic the function in kProcessor of coversion from str_to_int
def str_to_int(kmer):
    _map = {"A": 0, "C": 1, "T": 2, "G": 3}
    strint = 0
    for n in kmer:
        curr = _map[n]
        strint = strint | curr
        strint = strint << 2

    return strint >> 2


def int_to_str(kmer, kSize):
    _map = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
    kmer_str = ""
    for i in range(kSize, 0, -1):
        base = (kmer >> (i*2-2)) & 3
        ch = _map[base]
        kmer_str += ch

    return kmer_str


# Function to match 2 kmers with common first Q*2 bits or First Q Chars
# Return bool value, True if match, False if not
def matchingKmers(kmer1_int, kmer2_int, kSize, Q):
    msk = (~(-1 << Q*2)) << (kSize*2 - Q*2)
    return ((kmer1_int ^ kmer2_int) & msk == 0)


kSize = 10
kList_str = {}
kList = {}
kf = kp.kDataFrameMQF(kSize)
for i in range(100, 1000, 10):
    r = random.randint(1, 50)
    kmer = int_to_str(i, kSize)
    kList_str[kmer] = r
    kList[i] = r
    kf.setCount(kmer, r)


minQ = 4
maxQ = kSize
stepQ = 1

it = kf.begin()
prev_kmer = it.getHashedKmer()

count = 1000
while it != kf.end():
    if count == 0:
        break
    else:
        count -= 1

    it.next()
    curr_kmer = it.getHashedKmer()
    for Q in range(maxQ, minQ, -stepQ):
        if(matchingKmers(curr_kmer, prev_kmer, kf.getkSize(), Q)):
            print("Matching Q%d %s & %s | TRUE" % (Q, int_to_str(
                prev_kmer, kSize), int_to_str(curr_kmer, kSize)))
        else:
            print("Matching Q%d %s & %s | FALSE" % (Q, int_to_str(
                prev_kmer, kSize), int_to_str(curr_kmer, kSize)))
    print("---------------------------------------")

    prev_kmer = curr_kmer
