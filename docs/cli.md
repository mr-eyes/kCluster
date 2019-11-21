# Command Line Interface

## Options

```
Usage: kSpider [OPTIONS] COMMAND [ARGS]...

Options:
  --version    Show the version and exit.
  -q, --quiet
  --help       Show this message and exit.

Commands:
  index_kmers     FASTA file indexing by Kmers
  pairwise        Generating pairwise matrices for single/multiple...
  cluster         Sequences clustering regarding user-selected virtualQs.
  dump            Dump sqlite database table to the stdout in TSV format.
  index_skipmers  FASTA file indexing by Skipmers
  preprocess_cds  Preprocess protein coding transcript to extract CDS
```

## Sub-Commands

### 1- index_kmers

**index_kmers** subcommand is responsible for indexing FASTA/Q file using default substrings (kmers).

```
Usage: kSpider index_kmers [OPTIONS]

  FASTA file indexing by Kmers

Options:
  -f, --fasta PATH               FASTA file  [required]
  -n, --names PATH               Names file  [required]
  -k, --kmer-size INTEGER RANGE  kmer size  [required]
  -o, --output TEXT              index output file prefix
  --help                         Show this message and exit.
```

### 2- index_skipmers

**index_skipmers** subcommand is responsible for indexing FASTA/Q file using kmers skipping method called [*skipmers*](https://www.biorxiv.org/content/biorxiv/early/2017/09/19/179960.full.pdf)
You can select specific Open Reading Frame or processing all the ORFs by default. 

```
Usage: kSpider index_skipmers [OPTIONS]

  FASTA file indexing by Skipmers

Options:
  -f, --fasta PATH            FASTA file  [required]
  -n, --names PATH            Names file  [required]
  -k, --kmer-size INTEGER     kmer size  [required]
  -m, --cycle-bases INTEGER   used bases per cycle  [required]
  -n, --cycle-length INTEGER  kmer size(cycle length  [required]
  --orf INTEGER               select ORF <1,2,3>
  -o, --output TEXT           index output file prefix
  --help                      Show this message and exit.

```

### 3-Pairwise

The **pairwise** subcommand generates a pairwise distance similarity between each to sequences in a sqlite file.
The pairwise matrix contains number of shared virtualQs and number of kmers in the smalled sequence alongside another meta information.

```
Usage: kSpider pairwise [OPTIONS]

  Generating pairwise  matrices for single/multiple virtualQs.

Options:
  -m, --min-q INTEGER            minimum virtualQ  [default: 5]
  -M, --max-q INTEGER            maximum virtualQ
  -s, --step-q INTEGER           virtualQs range step  [default: 2]
  -i, --index-prefix TEXT        kProcessor index file prefix  [required]
  -o, --output-prefix TEXT       virtualQs output file(s) prefix
  --force                        Overwrite the already proessed virtualQs
  --backup                       Back up old virtualQs
  --export-colors [json|pickle]  export supercolors data [debugging purposes]
  --help                         Show this message and exit.
```

### 4- Dump

The **dump** subcommands exports the sqlite database table into a TSV file.

```
Usage: kSpider dump [OPTIONS]

  Dump sqlite database table to the stdout in TSV format.

Options:
  -d, --db PATH                   sqlite database file  [required]
  -t, --table [virtualQs|meta_info|namesmap] database table to be exported  [default: virtualQs]
  --help                          Show this message and exit.
```

### 5- Clustering

**cluster** subcommand is used for clustering the pairwise matrix generated from the `pairwise` subcommand.

```
Usage: kSpider cluster [OPTIONS]

  Sequences clustering regarding user-selected virtualQs.

Options:
  -m, --min-q INTEGER       minimum virtualQ
  -M, --max-q INTEGER       maximum virtualQ
  -s, --step-q INTEGER      virtualQs range step
  -c, --cutoff FLOAT RANGE  cluster sequences with (similarity > cutoff)
                            [default: 0.0]
  -d, --db PATH             sqlite database file  [required]
  --help                    Show this message and exit.
```


### 6- Preprocess CDS

**preprocess_cds** subcommand process protein coding transcript to extract CDS sequences in FASTA format alongside names_file to be ready for indexing

This command is recommended to be used before **index_skipmers** with `--orf 1`.

```
Usage: kSpider preprocess_cds [OPTIONS]

  Preprocess protein coding transcript to extract CDS

Options:
  -f, --fasta PATH              FASTA file  [required]
  -t, --diff-threshold INTEGER  minimum length difference %
  -l, --cds-len INTEGER         Minimum CDS length
  -s, --strand                  Examine only the top strand
  --help                        Show this message and exit.
```