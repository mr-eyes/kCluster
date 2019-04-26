import kProcessor as kp
import random
import hashlib
import json
import argparse
import sys
import os

VERBOSE = 0

class virtualQs:

    __kSize = 0

    def __init__(self, index_file_path):
        """VirtualQs class constructor.

        Args:
            index_file_name (string): coloredKDataFrame index file.

        """

        # if ".mqf" in index_file_path:
        #     exit("pass the index file prefix without extension")
        
        if os.path.exists(index_file_path + ".mqf"):
            print("Checking the file: ", index_file_path + ".mqf")
            #sys.stderr.write()
            #exit("index file does not exist")

        self.kf = kp.kDataFrame.load(index_file_path)
        self.__kSize = self.kf.getkSize()

        if self.__kSize == 0:
            exit("error loading the index")

    #@staticmethod
    def __mask(self,Q):
        return (~(-1 << Q*2)) << (self.__kSize*2 - Q*2)


    def set_params(self, minQ, maxQ, stepQ = 1):
        """virtualQs parameters setting.

        Args:
            minQ (int): minimum virtual Q (>= 1).
            maxQ (int): minimum virtual Q (<= kmer size).
            stepQ (int): virtual Q step (< maxQ)

        """

        if maxQ > self.__kSize:
            print("maxQ should't exceed the kmer Size, auto resetting Q=kSize", file=sys.stderr)
            self.__maxQ = self.__kSize
        else:
            self.__maxQ = maxQ
        
        if(minQ < 1):
            print("auto resetting minQ to 1", file=sys.stderr)
            self.__minQ = 1
        else:
            self.__minQ = minQ
        
        if(stepQ < 1):
            print("auto resetting Q step to 1", file=sys.stderr)
            self.__minQ = 1
        else:
            self.__stepQ = stepQ

        self.superColors = {}
        self.temp_superColors = {}
        self.superColorsCount = {}
        self.__masks = {}
        
        # Determine minQ and maxQ and get list of masks & superColorsDIct initialization
        # for Q in range(maxQ, minQ-1, -stepQ):
        for Q in range(self.__minQ, self.__maxQ + 1, self.__stepQ):
            print("َQ:  ", Q)
            self.__masks[Q] = self.__mask(Q)
            self.superColors[Q] = {}
            self.superColorsCount[Q] = {}
            self.temp_superColors[Q] = []


    # Take list of colors, sort it, and create a return a hash value
    @staticmethod
    def create_super_color(colors):
        return hashlib.md5(str(sorted(list(set(colors)))).encode()).hexdigest()[:9]

    @property
    def get_params(self):
        return {"minQ": self.__minQ, "maxQ": self.__maxQ, "stepQ": self.__stepQ}

    @property
    def kSize(self):
        return(self.__kSize)

    @property
    def masks(self):
        return self.__masks

    @staticmethod
    def int_to_str(kmer, kSize):
        _map = {0: 'A', 1: 'C', 2: 'T', 3: 'G'}
        kmer_str = ""
        for i in range(kSize, 0, -1):
            base = (kmer >> (i*2-2)) & 3
            ch = _map[base]
            kmer_str += ch
        
        return kmer_str


VQ = virtualQs(index_file_path="idx_min_test")
VQ.set_params(minQ = 8, maxQ = VQ.kSize, stepQ = 1)
print(VQ.get_params)

it = VQ.kf.begin()
prev_kmer = it.getHashedKmer()
prev_kmer_color = it.getKmerCount()

# Create list of kmers
while it != VQ.kf.end():    
    it.next()
    curr_kmer = it.getHashedKmer()
    curr_kmer_color = it.getKmerCount()

    # Apply XOR to kmer1 and kmer2 (single time per iteration)
    xor = prev_kmer ^ curr_kmer

    # Apply all masks with all Qs
    for Q, MASK in VQ.masks.items():
        matched = not bool(xor & MASK)  # True if there's match, False if not

        if matched:
            print("Matching Q%d %s & %s | TRUE | [prevC:%d, currC=%d]" % (Q, VQ.int_to_str(
                prev_kmer, VQ.kSize), VQ.int_to_str(curr_kmer, VQ.kSize), prev_kmer_color, curr_kmer_color))

            VQ.temp_superColors[Q] += [prev_kmer_color, curr_kmer_color]

        else:
            print("Matching Q%d %s & %s | FALSE | [prevC:%d, currC=%d]" % (Q, VQ.int_to_str(
                prev_kmer, VQ.kSize), VQ.int_to_str(curr_kmer, VQ.kSize),  prev_kmer_color, curr_kmer_color))
            VQ.temp_superColors[Q].append(prev_kmer_color)
            super_color_id = VQ.create_super_color(VQ.temp_superColors[Q])

            # Check if the superColor already exist
            # If yes: increment the count to one
            # If No:  Insert the new superColor and set the count to 1
            if super_color_id not in VQ.superColors[Q]:
                if len(VQ.temp_superColors[Q]) == 0:
                    print("Inserting an empty temp_superColors[%d]" % Q)
                VQ.superColors[Q][super_color_id] = list(set(VQ.temp_superColors[Q]))
                VQ.superColorsCount[Q][super_color_id] = 1

            else:
                # IF the supercolor already exist, just increment it
                VQ.superColorsCount[Q][super_color_id] += 1

            VQ.temp_superColors[Q] = [curr_kmer_color]

    #print(superColors[1])
    print("+++++++++++++++++++++++++++++++++++++++++++++++")

    prev_kmer = curr_kmer
    prev_kmer_color = curr_kmer_color

# If the last iteration got a match, push it to the superColors

for Q, colors in VQ.temp_superColors.items():
    colors.remove(curr_kmer_color)
    if not len(colors):
        continue
    super_color_id = VQ.create_super_color(colors)

    if super_color_id not in VQ.superColors[Q]:

        VQ.superColors[Q][super_color_id] = list(set(colors))
        VQ.superColorsCount[Q][super_color_id] = 1

    else:
        # IF the supercolor already exist, just increment it
        VQ.superColorsCount[Q][super_color_id] += 1


print("SUPER COLORS")
print(json.dumps(VQ.superColors, sort_keys=True, indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")
print("SUPER COLORS COUNT")
print(json.dumps(VQ.superColorsCount, sort_keys=True, indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")
print("TEMP SUPER COLORS")
print(json.dumps(VQ.temp_superColors, sort_keys=True,indent=4, separators=(',', ': ')))
print("\n--------------------------------\n\n")

print(VQ.superColors)

print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")
print("\n--------------------------------\n\n")

print(VQ.superColorsCount)
