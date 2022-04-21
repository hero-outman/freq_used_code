##################################################################
                        rgt-hint footprinting
##################################################################

usage: rgt-hint footprinting [-h] [--organism STRING] [--hmm-file FILE]
                             [--bias-table FILE_F,FILE_R] [--paired-end]
                             [--bias-correction] [--bias-type STRING]
                             [--output-location PATH] [--output-prefix STRING]
                             [--atac-seq] [--dnase-seq] [--histone]
                             [reads.bam regions.bed [reads.bam regions.bed ...]]

positional arguments:
  reads.bam regions.bed
                        BAM file of reads and BED files of interesting regions

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19,
                        hg38. mm9, and mm10. DEFAULT: hg19
  --hmm-file FILE       If the argument is not given, then a default HMM will
                        be used.
  --bias-table FILE_F,FILE_R
                        List of files with all possible k-mers (for any k) and
                        their bias estimates. Each line should contain a kmer
                        and the bias estimate separated by tab.
  --paired-end          Set it if your ATAC-seq data is paired-end sequenced.
                        Note that this option is only applied to ATAC-seq
                        data. DEFAULT: False
  --bias-correction     If set, footprint calling will based on bias corrected
                        DNase-seq signal. This option is only applied to
                        DNase-seq. DEFAULT: False
  --bias-type STRING    Type of protocol used to generate the DNase-seq.
                        Available options are: 'SH' (DNase-seq single-hit
                        protocol), 'DH' (DNase-seq double-hit protocol).
                        DEFAULT: SH
  --output-location PATH
                        Path where the output bias table files will be
                        written. DEFAULT: current directory
  --output-prefix STRING
                        The prefix for results files. DEFAULT: footprints
  --atac-seq            If set, footprint calling will based on ATAC-seq
                        model. DEFAULT: False
  --dnase-seq           If set, footprint calling will based on DNase-seq
                        model. DEFAULT: False
  --histone             If set, footprint calling will based on histone
                        modification model. DEFAULT: False


##################################################################
                        rgt-hint tracks
##################################################################

usage: rgt-hint tracks [-h] [--organism STRING] [--bias-table FILE1_F,FILE1_R]
                       [--raw] [--bc] [--norm] [--bigWig] [--strand-specific]
                       [--output-location PATH] [--output-prefix STRING]
                       [reads.bam regions.bed [reads.bam regions.bed ...]]

positional arguments:
  reads.bam regions.bed
                        BAM file of reads and BED files of interesting regions

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19,
                        hg38. mm9, and mm10. DEFAULT: hg19
  --bias-table FILE1_F,FILE1_R
                        Bias table files used to generate bias corrected
                        tracks. DEFAULT: None
  --raw                 If set, the raw signals from DNase-seq or ATAC-seq
                        data will be generated. DEFAULT: False
  --bc                  If set, the bias corrected signals from DNase-seq or
                        ATAC-seq data will be generated. DEFAULT: False
  --norm                If set, the normalised signals from DNase-seq or ATAC-
                        seq data will be generated. DEFAULT: False
  --bigWig              If set, all .wig files will be converted to .bw files.
                        DEFAULT: False
  --strand-specific     If set, the tracks will be splitted into two files,
                        one for forward and another for reverse strand.
                        DEFAULT: False
  --output-location PATH
                        Path where the output bias table files will be
                        written. DEFAULT: current directory
  --output-prefix STRING
                        The prefix for results files. DEFAULT: tracks



##################################################################
                        rgt-motifanalysis matching
##################################################################
 
usage: rgt-motifanalysis matching [-h] [--organism STRING] [--fpr FLOAT]
                                  [--pseudocounts FLOAT]
                                  [--rand-proportion FLOAT] [--norm-threshold]
                                  [--motif-dbs PATH [PATH ...]]
                                  [--remove-strand-duplicates] [--rmdup]
                                  [--filter KEY_VALUE_PATTERN]
                                  [--filter-type {inexact,exact,regex}]
                                  [--target-genes PATH] [--make-background]
                                  [--promoter-length INT]
                                  [--output-location PATH] [--bigbed]
                                  [--normalize-bitscore]
                                  (--input-matrix matrix.txt | --promoters-only | --input-files regions.bed [regions.bed ...])

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19
                        or hg38.
  --fpr FLOAT           False positive rate cutoff.
  --pseudocounts FLOAT  Pseudocounts to be added to raw counts of each PFM.
  --rand-proportion FLOAT
                        If set, a random regions file will be created (eg, for
                        later enrichment analysis). The number of coordinates
                        will be equal to this value times the size of the
                        inputregions. We advise you use a value of at least
                        10.
  --norm-threshold      If this option is used, the thresholds for all PWMs
                        will be normalized by their length. In this scheme,
                        the threshold cutoff is evaluated in the regular way
                        by the given fpr. Then, all thresholds are divided by
                        the length of the motif. The final threshold consists
                        of the average between all normalized motif
                        thresholds. This single threshold will be applied to
                        all motifs.
  --motif-dbs PATH [PATH ...]
                        New 'motif DB' folders to use instead of the ones
                        within the RGTDATA folder. Each folder must contain
                        PWM files.
  --remove-strand-duplicates
                        Certain motifs are 'palindromic', or more specifically
                        they have a palindromic consensus sequence. When this
                        happens, the output MPBS file will have duplicates:
                        same chromosome and initial and final position, but
                        opposing strand. Select this option to only retain the
                        'strand duplicate' with the highest score. Duplicates
                        due to overlapping input regions are NOT affected by
                        this.
  --rmdup               Remove any duplicate region from the input BED files.
  --filter KEY_VALUE_PATTERN
                        List of key-value patterns to select a subset of TFs
                        using the metadata (MTF files), e.g. for Mouse and
                        Human on Selex data use:
                        "species:sapiens,mus;data_source:selex". NB: the
                        DATABASE values must be written in full - exact
                        matching is always performed.Valid key types are
                        "name", "gene_names", "family", "uniprot_ids",
                        "data_source", "tax_group", "species", "database",
                        "name_file" and "gene_names_file"
  --filter-type {inexact,exact,regex}
                        Only useful together with the --filter argument.Exact
                        will only match perfect matching of the value for each
                        key. Inexact will match in case the value pattern is
                        contained within the motif. Regex allows for a more
                        complex pattern use.
  --input-matrix matrix.txt
                        The experimental matrix allows the specification of
                        gene-association rules among input files (see online
                        documentation for details).
  --promoters-only      If you ONLY want to perform promoter matching without
                        providing any input file/matrix. If --target-genes is
                        not provided, then all available promoters will be
                        matched against. Note how this makes '--make-
                        background' redundant.
  --input-files regions.bed [regions.bed ...]
                        BED files to perform motif matching on.

Promoter-regions matching:
  These arguments are only used with the --promoters-only option (for the
  purpose of matching only on the promoters of all or a subset of genes)

  --target-genes PATH   List of genes (one per line) to get the promoter
                        regions from.
  --make-background     If set, it will perform motif matching on the
                        'background regions', composed of the promoters of all
                        available genes (minus the target genes, if
                        specified). It doesn't require --target-genes.
  --promoter-length INT
                        Length of the promoter region (in bp) to be extracted
                        from each gene.

Output:
  Where to put the output files and how to post-process them.

  --output-location PATH
                        Path where the output MPBS files will be written.
                        Defaults to 'match' in the current directory.
  --bigbed              If this option is used, all bed files will be written
                        as bigbed.
  --normalize-bitscore  In order to print bigbed files the scores need to be
                        normalized between 0 and 1000. Don't use this option
                        if real bitscores should be printed in the resulting
                        bed file. Without this option, bigbed files will never
                        be created.                        



##################################################################
                        rgt-motifanalysis enrichment
##################################################################
usage: rgt-motifanalysis enrichment [-h] [--organism STRING]
                                    [--matching-location PATH]
                                    [--use-only-motifs PATH]
                                    [--input-matrix PATH]
                                    [--multiple-test-alpha FLOAT]
                                    [--motif-dbs PATH [PATH ...]]
                                    [--filter KEY_VALUE_PATTERN]
                                    [--filter-type {inexact,exact,regex}]
                                    [--promoter-length INT]
                                    [--maximum-association-length INT]
                                    [--exclude-target-genes]
                                    [--output-location PATH]
                                    [--print-thresh FLOAT] [--bigbed]
                                    [--logo-copy | --logo-embed]
                                    background.bed [input.bed [input.bed ...]]

positional arguments:
  background.bed        BED file containing background regions.
  input.bed             BED files to be enriched against the background.

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19
                        or hg38.
  --matching-location PATH
                        Directory where the matching output containing the
                        MPBS files resides. Defaults to 'match' in the current
                        directory.
  --use-only-motifs PATH
                        Only use the motifs contained within this file (one
                        for each line).
  --input-matrix PATH   If an experimental matrix is provided, the input
                        arguments will be ignored.
  --multiple-test-alpha FLOAT
                        Alpha value for multiple test.
  --motif-dbs PATH [PATH ...]
                        New 'motif DB' folders to use instead of the ones
                        within the RGTDATA folder. Each folder must contain
                        PWM files.
  --filter KEY_VALUE_PATTERN
                        List of key-value patterns to select a subset of TFs
                        using the metadata (MTF files), e.g. for Mouse and
                        Human on Selex data use:
                        "species:sapiens,mus;data_source:selex". NB: the
                        DATABASE values must be written in full - exact
                        matching is always performed.Valid key types are
                        "name", "gene_names", "family", "uniprot_ids",
                        "data_source", "tax_group", "species", "database",
                        "name_file" and "gene_names_file"
  --filter-type {inexact,exact,regex}
                        Only useful together with the --filter argument.Exact
                        will only match perfect matching of the value for each
                        key. Inexact will match in case the value pattern is
                        contained within the motif. Regex allows for a more
                        complex pattern use.
  --logo-copy           The logos are copied to a local directory. The HTML
                        report will contain relative paths to this directory.
  --logo-embed          The logos are embedded directly into the HTML report.

Promoter-regions enrichment:
  Used both for gene set via experimental matrix (see documentation), and
  for reporting the gene names associated to each motif.

  --promoter-length INT
                        Length of the promoter region (in bp) to be extracted
                        from each gene.
  --maximum-association-length INT
                        Maximum distance between a coordinate and a gene (in
                        bp) in order for the former to be considered
                        associated with the latter.
  --exclude-target-genes
                        If set the specified target genes areexcluded from
                        background file

Output:
  Where to put the output files and how to post-process them.

  --output-location PATH
                        Path where the output MPBS files will be written.
                        Defaults to 'enrichment' in the current directory.
  --print-thresh FLOAT  Only MPBSs whose factor's enrichment corrected p-value
                        are less than equal this option are printed. Use 1.0
                        to print all MPBSs.
  --bigbed              If this option is used, all bed files will be written
                        as bigbed.                        


rgt-hint differential -h
usage: rgt-hint differential [-h] [--organism STRING]
                             [--mpbs-files FILE1,FILE2...]
                             [--reads-files FILE1,FILE2...]
                             [--conditions STRING] [--colors STRING]
                             [--window-size INT] [--fdr FLOAT] [--bc]
                             [--nc INT] [--output-location PATH]
                             [--output-prefix STRING] [--standardize]
                             [--output-profiles]

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19,
                        hg38. mm9, and mm10. DEFAULT: hg19
  --mpbs-files FILE1,FILE2...
                        Predicted motif binding sites for each condition.Files
                        should be separated with comma.
  --reads-files FILE1,FILE2...
                        Reads for each condition. Files should be separated
                        with comma.
  --conditions STRING   Name for each condition. DEFAULT:
                        condition1,condition2, ...
  --colors STRING       Set color in line plot. DEFAULT: None, ...
  --window-size INT     The window size for differential analysis. DEFAULT:
                        200
  --fdr FLOAT           The false discovery rate. DEFAULT: 0.05
  --bc                  If set, all analysis will be based on bias corrected
                        signal. DEFAULT: False
  --nc INT              The number of cores. DEFAULT: 1
  --output-location PATH
                        Path where the output bias table files will be
                        written. DEFAULT: current directory
  --output-prefix STRING
                        The prefix for results files. DEFAULT: differential
  --standardize         If set, the signal will be rescaled to (0, 1) for
                        plotting.
  --output-profiles     If set, the footprint profiles will be writen into a
                        text, in which each row is a specific instance of the
                        given motif. DEFAULT: False


rgt-hint differential -h
usage: rgt-hint differential [-h] [--organism STRING]
                             [--mpbs-files FILE1,FILE2...]
                             [--reads-files FILE1,FILE2...]
                             [--conditions STRING] [--colors STRING]
                             [--window-size INT] [--fdr FLOAT] [--bc]
                             [--nc INT] [--output-location PATH]
                             [--output-prefix STRING] [--standardize]
                             [--output-profiles]

optional arguments:
  -h, --help            show this help message and exit
  --organism STRING     Organism considered on the analysis. Must have been
                        setup in the RGTDATA folder. Common choices are hg19,
                        hg38. mm9, and mm10. DEFAULT: hg19
  --mpbs-files FILE1,FILE2...
                        Predicted motif binding sites for each condition.Files
                        should be separated with comma.
  --reads-files FILE1,FILE2...
                        Reads for each condition. Files should be separated
                        with comma.
  --conditions STRING   Name for each condition. DEFAULT:
                        condition1,condition2, ...
  --colors STRING       Set color in line plot. DEFAULT: None, ...
  --window-size INT     The window size for differential analysis. DEFAULT:
                        200
  --fdr FLOAT           The false discovery rate. DEFAULT: 0.05
  --bc                  If set, all analysis will be based on bias corrected
                        signal. DEFAULT: False
  --nc INT              The number of cores. DEFAULT: 1
  --output-location PATH
                        Path where the output bias table files will be
                        written. DEFAULT: current directory
  --output-prefix STRING
                        The prefix for results files. DEFAULT: differential
  --standardize         If set, the signal will be rescaled to (0, 1) for
                        plotting.
  --output-profiles     If set, the footprint profiles will be writen into a
                        text, in which each row is a specific instance of the
                        given motif. DEFAULT: False                                                