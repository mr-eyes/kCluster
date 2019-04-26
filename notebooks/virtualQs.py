import kProcessor as kp
import random
import hashlib
import json
import argparse
import sys
import os
import pickle

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
    def __mask(self, Q):
        return (~(-1 << Q*2)) << (self.__kSize*2 - Q*2)

    def set_params(self, minQ, maxQ, stepQ=1):
        """virtualQs parameters setting.

        Args:
            minQ (int): minimum virtual Q (>= 1).
            maxQ (int): minimum virtual Q (<= kmer size).
            stepQ (int): virtual Q step (< maxQ)

        """

        if maxQ > self.__kSize:
            print(
                "maxQ should't exceed the kmer Size, auto resetting Q=kSize", file=sys.stderr)
            self.__maxQ = self.__kSize

        elif maxQ == 0:
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

    def export_superColors(self, prefix, Q, method="json"):
        if Q not in self.superColors and Q not in self.superColorsCount:
            exit("virtualQ: %d does not exist" % Q)

        if method == "pickle":
            suffix = ".pickle"
        elif method == "json":
            suffix = ".json"
        else:
            exit("export only in [pickle,json]")

        virtualQs_file_name = prefix + "_" + str(Q) + suffix
        virtualQs_count_file_name = prefix + "_" + str(Q) + "_counts" + suffix

        if method == "pickle":
            print("writing virtual Q %d pickles ..." % Q)
            with open(virtualQs_file_name, "wb") as f:
                pickle.dump(self.superColors[Q], f, pickle.HIGHEST_PROTOCOL)

            with open(virtualQs_count_file_name, "wb") as f:
                pickle.dump(
                    self.superColorsCount[Q], f, pickle.HIGHEST_PROTOCOL)

        elif method == "json":
            with open(virtualQs_file_name, "w") as f:
                f.write(json.dumps(
                    self.superColors[Q], sort_keys=True, indent=4, separators=(',', ': ')))

            with open(virtualQs_file_name, "w") as f:
                f.write(json.dumps(
                    self.superColorsCount[Q], sort_keys=True, indent=4, separators=(',', ': ')))

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


def construct_virtualQs(min_q, max_q, step_q, index_prefix, output_prefix, output_type, force_write = True):
    VQ = virtualQs(index_file_path=index_prefix)
    VQ.set_params(minQ=min_q, maxQ=max_q, stepQ=step_q)
    
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
            # True if there's match, False if not
            matched = not bool(xor & MASK)

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
                    VQ.superColors[Q][super_color_id] = list(
                        set(VQ.temp_superColors[Q]))
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
    print(json.dumps(VQ.superColorsCount, sort_keys=True,
                    indent=4, separators=(',', ': ')))
    print("\n--------------------------------\n\n")
    print("TEMP SUPER COLORS")
    print(json.dumps(VQ.temp_superColors, sort_keys=True,
                    indent=4, separators=(',', ': ')))
    print("\n--------------------------------\n\n")

    print(VQ.superColors)

    print("\n--------------------------------\n\n")
    print("\n--------------------------------\n\n")
    print("\n--------------------------------\n\n")
    print("\n--------------------------------\n\n")

    print(VQ.superColorsCount)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', action='store', dest='index_prefix',
                        help='index file prefix')

    parser.add_argument('-m', action='store', dest='minQ',
                        help='minimum Q')

    parser.add_argument('-M', action='store', dest='maxQ',
                        help='maximum Q, 0 for Q=kSize')

    parser.add_argument('-s', action='store', dest='stepQ',
                        help='Q step')

    parser.add_argument('-f', action='store', dest='force',
                        help='force rewrite the written virtualQs files')

    parser.add_argument('-e', action='store', dest='output_type',
                        help='output type [json|pickle] default:json --optional')

    parser.add_argument('-o', action='store', dest='output_prefix',
                        help='virtualQs output files prefix --optional')

    if len(sys.argv)==1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.index_prefix:
        if os.path.isfile(args.index_prefix + ".mqf"):
            index_file_path = args.index_prefix
        else:
            sys.exit('Index prefix ' + args.index_prefix + 'Does not exist!')
    else:
        exit("must provide index file prefix, use -i <index_file_prefix>")

    if args.minQ:
        minQ = int(args.minQ)
    else:
        exit("must provide minimum Q, use -m <int>")

    if args.maxQ:
        maxQ = int(args.maxQ)
    else:
        exit("must provide maximum Q, use -M <int>")

    if args.stepQ:
        stepQ = int(args.stepQ)
    else:
        exit("must provide Q step, use -s <int>")

    if args.output_prefix:
        output_prefix = args.output_prefix
    else:
        from random import randint as rint
        output_prefix = str(rint(1, 100))

    if args.output_type:
        output_type = args.output_type
    else:
        output_type = "json"

    construct_virtualQs(minQ, maxQ, stepQ, index_file_path, output_prefix, output_type)

