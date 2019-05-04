import kProcessor as kp
import random
import hashlib
import json
import argparse
import sys
import os
import pickle

# TODO use logging instead of normal prints.
# TODO use Collections:Defaultdict for optimized performance.
# TODO use NUMBA for optimizing maths and loops

class virtualQs:
    """Holds the superColors and superColorsCount tables."""

    __kSize = None


    def __init__(self, index_file_path: str):
        """VirtualQs class constructor.

        Args:
            index_file_name (string): coloredKDataFrame index file.

        """

        self.kf = kp.kDataFrame.load(index_file_path)
        self.__kSize = self.kf.getkSize()

        if self.__kSize is None:
            print("error loading the index", file = sys.stderr)
            sys.exit(1)

    def __mask(self, Q):
        """create a bit mask given kmer size and Q value."""
        return (~(-1 << Q*2)) << (self.__kSize*2 - Q*2)

    def set_params(self, minQ: int, maxQ: int, stepQ: int = 2):
        """virtualQs parameters setting.

        Args:
            minQ (int): minimum virtual Q (>= 1).
            maxQ (int): minimum virtual Q (<= kmer size).
            stepQ (int): virtual Q step (< maxQ)

        """

        if maxQ > self.__kSize:
            print(
                "*WARNING* maxQ should't exceed the kmer Size, auto reinitializing Q with kSize %d" % (self.__kSize), file=sys.stderr)
            self.__maxQ = self.__kSize

        elif maxQ is 0:
            self.__maxQ = self.__kSize

        else:
            self.__maxQ = maxQ

        if(minQ < 1):
            print("*WARNING* minQ shouldn't be less than 1, auto reinitializing minQ to 5", file=sys.stderr)
            self.__minQ = 5
        elif minQ > self.__maxQ:
            print("*WARNING* minQ shouldn't exceed the maxQ, auto reinitializing minQ to maxQ", file=sys.stderr)
            self.__minQ = self.__maxQ
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
            self.__masks[Q] = self.__mask(Q)
            self.superColors[Q] = {}
            self.superColorsCount[Q] = {}
            self.temp_superColors[Q] = []

    def export_superColors(self, prefix, Q, method="json"):
        """superColors table exporting

        Retrieves Q value needs to be exported and output file format 

        Args:
            prefix: exported file name prefix.
            Q: Q value to be extracted from the superColors tables.
            method: specify the output file format pickle or json.  
        """

        if Q not in self.superColors and Q not in self.superColorsCount:
            print("virtualQ: {} does not exist".format(Q), file=sys.stderr)
            sys.exit(1)
        else:
            print("Writing")
        
        

        if method == "pickle":
            suffix = ".pickle"
        elif method == "json":
            suffix = ".json"
        else:
            print("export only in [pickle,json]", file=sys.stderr)
            sys.exit(1)

        virtualQs_file_name = prefix + "_Q" + str(Q) + suffix
        virtualQs_count_file_name = prefix + "_Q" + str(Q) + "_counts" + suffix

        if method == "pickle":
            print("writing virtual Q %d pickles ..." % Q)
            with open(virtualQs_file_name, "wb") as f:
                pickle.dump(self.superColors[Q], f, pickle.HIGHEST_PROTOCOL)

            with open(virtualQs_count_file_name, "wb") as f:
                pickle.dump(
                    self.superColorsCount[Q], f, pickle.HIGHEST_PROTOCOL)

        elif method == "json":
            with open(virtualQs_file_name, "w") as f:
                print("Writing")
                f.write(json.dumps(
                    self.superColors[Q], sort_keys=True, indent=4, separators=(',', ': ')))

            with open(virtualQs_count_file_name, "w") as f:
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


def construct_virtualQs(min_q, max_q, step_q, index_prefix, output_prefix, output_type, force_write=True):
    VQ = virtualQs(index_file_path=index_prefix)
    VQ.set_params(minQ=min_q, maxQ=max_q, stepQ=step_q)

    #print("Constructing virtualQs with params: ", VQ.get_params, file = sys.stderr)

    it = VQ.kf.begin()
    prev_kmer = it.getHashedKmer()
    prev_kmer_color = it.getKmerCount()

    # Iterate over all kmers.
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
                VQ.temp_superColors[Q] += [prev_kmer_color, curr_kmer_color]

                # print("Matching Q%d %s & %s | TRUE | [prevC:%d, currC=%d]" % (Q, VQ.int_to_str(
                #     prev_kmer, VQ.kSize), VQ.int_to_str(curr_kmer, VQ.kSize), prev_kmer_color, curr_kmer_color))


            else:
                VQ.temp_superColors[Q].append(prev_kmer_color)
                super_color_id = VQ.create_super_color(VQ.temp_superColors[Q])

                # print("Matching Q%d %s & %s | FALSE | [prevC:%d, currC=%d]" % (Q, VQ.int_to_str(
                #     prev_kmer, VQ.kSize), VQ.int_to_str(curr_kmer, VQ.kSize),  prev_kmer_color, curr_kmer_color))


                # Check if the superColor already exist
                # If yes: increment the count to one
                # If No:  Insert the new superColor and set the count to 1
                if super_color_id not in VQ.superColors[Q]:

                    VQ.superColors[Q][super_color_id] = list(set(VQ.temp_superColors[Q]))

                    VQ.superColorsCount[Q][super_color_id] = 1

                else:
                    # IF the supercolor already exist, just increment it
                    VQ.superColorsCount[Q][super_color_id] += 1

                VQ.temp_superColors[Q] = [curr_kmer_color]

        #print(superColors[1])
        #print("+++++++++++++++++++++++++++++++++++++++++++++++")

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

    # print("SUPER COLORS")

    # print(json.dumps(VQ.superColors, sort_keys=True, indent=4, separators=(',', ': ')))
    # print("\n--------------------------------\n\n")
    # print("SUPER COLORS COUNT")
    # print(json.dumps(VQ.superColorsCount, sort_keys=True,
    #                 indent=4, separators=(',', ': ')))
    # print("\n--------------------------------\n\n")
    # print("TEMP SUPER COLORS")
    # print(json.dumps(VQ.temp_superColors, sort_keys=True,
    #                 indent=4, separators=(',', ': ')))
    # print("\n--------------------------------\n\n")

    # print(VQ.superColors)

    # print("\n--------------------------------\n\n")
    # print("\n--------------------------------\n\n")
    # print("\n--------------------------------\n\n")
    # print("\n--------------------------------\n\n")

    # print(VQ.superColorsCount)

    # Save all Qs to files.

    _params = VQ.get_params
    print(_params)
    for Q in range(_params["minQ"], _params["maxQ"] + 1, _params["stepQ"]):
        VQ.export_superColors(output_prefix, Q, output_type)


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', action='store', dest='index_prefix',
                        help='index file prefix <str>')

    parser.add_argument('-m', action='store', dest='minQ',
                        help='minimum Q <int>')

    parser.add_argument('-M', action='store', dest='maxQ',
                        help='maximum Q <int>, 0 for Q=kSize')

    parser.add_argument('-s', action='store', dest='stepQ',
                        help='Q step <int>')

    parser.add_argument('-f', action='store', dest='force',
                        help='force rewrite the written virtualQs files --optional')

    parser.add_argument('-e', action='store', dest='output_type',
                        help='output type <json|pickle> default:json --optional')

    parser.add_argument('-o', action='store', dest='output_prefix',
                        help='virtualQs output files prefix <str> --optional')

    if len(sys.argv) is 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    args = parser.parse_args()

    if args.index_prefix:
        if os.path.isfile(args.index_prefix + ".mqf"):
            index_file_path = args.index_prefix
        else:
            print('Index prefix ' + args.index_prefix + ' Does not exist!', file = sys.stderr)
            sys.exit(1)

    else:
        print("must provide index file prefix, use - i < index_file_prefix >", file=sys.stderr)
        sys.exit(1)

    if args.minQ:
        minQ = int(args.minQ)
    else:
        print("must provide minimum Q, use -m <int>", file=sys.stderr)
        sys.exit(1)

    if args.maxQ:
        maxQ = int(args.maxQ)
    else:
        print("must provide maximum Q, use -M <int>", file = sys.stderr)
        sys.exit(1)

    if args.stepQ:
        stepQ = int(args.stepQ)
    else:
        print("must provide Q step, use -s <int>", file = sys.stderr)
        sys.exit(1)

    if args.output_prefix:
        output_prefix = args.output_prefix
    else:
        from random import randint as rint
        output_prefix = str(rint(1, 100))

    if args.output_type:
        output_type = args.output_type
    else:
        output_type = "json"

    construct_virtualQs(minQ, maxQ, stepQ, index_file_path,
                        output_prefix, output_type)


if __name__ == '__main__':

    print("""
        ██╗  ██╗ ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ 
        ██║ ██╔╝██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗
        █████╔╝ ██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝
        ██╔═██╗ ██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗
        ██║  ██╗╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║
        ╚═╝  ╚═╝ ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝                                            
    \n\n""")

    main()
