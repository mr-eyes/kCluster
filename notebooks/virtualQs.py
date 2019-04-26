import kProcessor as kp
import random
import hashlib
import json

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


def mask(kSize, Q):
    return (~(-1 << Q*2)) << (kSize*2 - Q*2)

# Take list of colors, sort it, and create a return a hash value


def create_super_color(colors):
    val = hashlib.md5(str(sorted(list(set(colors)))).encode()).hexdigest()[:9]
    print("Hashing: %s to %s" % (str(sorted(colors)), val))
    return val


kf = kp.kDataFrame.load("idx_min_test")
kSize = kf.getkSize()

minQ = 8
stepQ = 1
maxQ = 8
superColors = {}
temp_superColors = {}
superColorsCount = {}
masks = {}

# Determine minQ and maxQ and get list of masks & superColorsDIct initialization

#for Q in range(maxQ, minQ-1, -stepQ):
for Q in range(minQ, maxQ + 1, stepQ):
    masks[Q] = mask(kSize, Q)
    superColors[Q] = {}
    superColorsCount[Q] = {}
    temp_superColors[Q] = []

#print(temp_superColors)
#print(superColors)
#print(superColorsCount)

#list(map(bin, masks.values()))

it = kf.begin()
prev_kmer = it.getHashedKmer()
prev_kmer_color = it.getKmerCount()

count = 15

# Create list of kmers
while it != kf.end():
    
    # if count == 0:
    #     break
    # count -= 1
    
    it.next()
    curr_kmer = it.getHashedKmer()
    curr_kmer_color = it.getKmerCount()

    # Apply XOR to kmer1 and kmer2 (single time per iteration)
    xor = prev_kmer ^ curr_kmer

    # Apply all masks with all Qs
    for Q, MASK in masks.items():

        matched = not bool(xor & MASK)  # True if there's match, False if not

        if matched:
            print("Matching Q%d %s & %s | TRUE | [prevC:%d, currC=%d]" % (Q, int_to_str(
                prev_kmer, kSize), int_to_str(curr_kmer, kSize), prev_kmer_color, curr_kmer_color))
            temp_superColors[Q] += [prev_kmer_color, curr_kmer_color]

        else:
            print("Matching Q%d %s & %s | FALSE | [prevC:%d, currC=%d]" % (Q, int_to_str(
                prev_kmer, kSize), int_to_str(curr_kmer, kSize),  prev_kmer_color, curr_kmer_color))
            temp_superColors[Q].append(prev_kmer_color)
            super_color_id = create_super_color(temp_superColors[Q])

            # Check if the superColor already exist
            # If yes: increment the count to one
            # If No:  Insert the new superColor and set the count to 1
            if super_color_id not in superColors[Q]:
                if len(temp_superColors[Q]) == 0:
                    print("Inserting an empty temp_superColors[%d]" % Q)
                superColors[Q][super_color_id] = list(set(temp_superColors[Q]))
                superColorsCount[Q][super_color_id] = 1

            else:
                # IF the supercolor already exist, just increment it
                superColorsCount[Q][super_color_id] += 1

            temp_superColors[Q] = [curr_kmer_color]

    #print(superColors[1])
    print("+++++++++++++++++++++++++++++++++++++++++++++++")

    prev_kmer = curr_kmer
    prev_kmer_color = curr_kmer_color

# If the last iteration got a match, push it to the superColors
if 1:
    for Q, colors in temp_superColors.items():
        colors.remove(curr_kmer_color)
        if not len(colors):
            continue
        super_color_id = create_super_color(colors)

        if super_color_id not in superColors[Q]:

            superColors[Q][super_color_id] = list(set(colors))
            superColorsCount[Q][super_color_id] = 1

        else:
            # IF the supercolor already exist, just increment it
            superColorsCount[Q][super_color_id] += 1


print("SUPER COLORS")
print(json.dumps(superColors, sort_keys=True, indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")
print("SUPER COLORS COUNT")
print(json.dumps(superColorsCount, sort_keys=True, indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")
print("TEMP SUPER COLORS")
print(json.dumps(temp_superColors, sort_keys=True, indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")

print(superColors)

print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")

print(superColorsCount)
