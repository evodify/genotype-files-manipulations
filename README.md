# Genotype calls files manipulations

Set of scripts to manipulate tab-delimited genotype calls files as well as to convert them to other popular formats.

All python scripts contain description of input and output data format in a header of each file.
To see possible options, run python script with --help option:
`python script.py --help`

Most of these scripts require the custom python module `calls`, so make sure that you also download and put the file `calls.py` in the same directory where your scripts are.

Examples of a tab-delimited genotype calls file (hereafter, tab file).

Two-character coded table (e.g. produced with [VariantsToTable](https://software.broadinstitute.org/gatk/gatkdocs/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php) from the GATK) :

```
CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   T/A ./. ./. A/A ./. ./. ./. ./.
chr_1   2   C   T/C T/C ./. C/C C/C ./. C/C ./.
chr_1   3   C   C/GCC   C/C ./. C/C C/C C/C C/C C/C
chr_1   4   T   T/T T/T ./. T/T T/T T/T T/T T/T
chr_2   1   A   A/A A/A ./. A/A A/A A/A A/A A/A
chr_2   2   C   C/C C/C ./. C/C C/C C/C C/C C/C
chr_2   3   C   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT   AT/AT
chr_2   4   C   C/C T/T C/C C/C C/C C/C C/C C/C
chr_2   5   T   T/T C/C T/T C/T T/T C/T T/T T/T
chr_3   1   G   G/G ./. ./. G/G ./. ./. ./. ./.
chr_3   2   C   G/C C/C ./. C/C C/C ./. C/C ./.
chr_3   3   CTT CTT/CTT CTT/C   CTT/C   CTT/CTT CTT/CTT CTT/CTT CTT/CTT CTT/CTT
chr_3   4   TA  T/T T/T ./. T/T T/T T/T T/T T/TA
chr_3   5   G   */* G/* ./. G/G G/G G/G C/C G/G

```

One-character coded tab file (heterozygous genotypes are represented by ambiguous characters R, Y, M, K, S, W):

```
CHROM   POS REF sample1 sample2 sample3 sample4 sample5 sample6 sample7 sample8
chr_1   1   A   W   N   N   A   N   N   N   N
chr_1   2   C   Y   Y   N   C   C   N   C   N
chr_1   3   C   N   C   N   C   C   C   C   C
chr_1   4   T   T   T   N   T   T   T   T   T
chr_2   1   A   A   A   N   A   A   A   A   A
chr_2   2   C   C   C   N   C   C   C   C   C
chr_2   3   C   N   N   N   N   N   N   N   N
chr_2   4   C   C   T   C   C   C   C   C   C
chr_2   5   T   T   C   T   Y   T   Y   T   T
chr_3   1   G   G   N   N   G   N   N   N   N
chr_3   2   C   S   C   N   C   C   N   C   N
chr_3   3   N   N   N   N   N   N   N   N   N
chr_3   4   N   T   T   N   T   T   T   T   N
chr_3   5   G   -   N   N   G   G   G   C   G
```

##

[FastaToPhylip.py](FastaToPhylip.py) converts FASTA to PHYLIP.

[FastaToTab.py](FastaToTab.py)  converts FASTA to tab-delimited file with columns: Chr, Pos, REF.

[assessNs.py](assessNs.py) calculates missing data (Ns) per position/sample and visualize the results.

[calculateNsPerWindow.py](calculateNsPerWindow.py) calculates number of positions with missing data (Ns) using the sliding window approach.

[calls.py](calls.py) is a custom python module. It is a dependecy for the most of the scripts listed here.

[calls_to_treeMix_input.py](calls_to_treeMix_input.py) outputs alleles counts file that is required as input for [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).

[callsToFastaPhy.py](callsToFastaPhy.py) converts genotype calls file to FASTA and PHYLIP.

[combine_overlapping_intervals.py](combine_overlapping_intervals.py) combines overlapping genetic intervals in the BED format.

[extractSIFT4Gannotation.py](extractSIFT4Gannotation.py) extracts the [SIFT4G annotation](http://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html) for a given set of samples according to their genotypes.

[filterByNs.py](filterByNs.py) removes all sites that consists of more than a given amount of missing data (Ns).

[find_popSpesificAlleles.py](find_popSpesificAlleles.py) outputs only unique allele of one population relative to another.

[keep_biallelic.py](keep_biallelic.py) removes sites with more than two alleles.

[make_input_MSMC.py](make_input_MSMC.py) makes input for [MSMC](https://github.com/stschiff/msmc).

[makeSweepFinderInput.py](makeSweepFinderInput.py) makes an input file for [SweepFinder](http://people.binf.ku.dk/rasmus/webpage/sf.html).

[mask_tab.py](mask_tab.py) removes the masked sites from a tab file. The masked sites are provided in a BED file.

[make_input_stairway_plot_v1_BS.py](make_input_stairway_plot_v1_BS.py) makes input files including bootstrap replicates for [Stairway](https://sites.google.com/site/jpopgen/stairway-plot) version 1.

[make_input_stairway_plot_v2.py](make_input_stairway_plot_v2.py) makes an input file for [Stairway](https://sites.google.com/site/jpopgen/stairway-plot) version 2.

[merge_SNP_wholeGenome_TabFiles.py](merge_SNP_wholeGenome_TabFiles.py) merges whole genome and SNPs tab files. This is needed because non-polymorphic sites and SNPs are filtered differently with [GATK](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php).

[mergeChrPos.py](mergeChrPos.py) merges all chromosomes into continious genomic coordinates.

[mergeTabFiles.py](mergeTabFiles.py) merges two tab files by their overlapping positions.

[polarizeGT.py](polarizeGT.py) polarizes the genotype data by keeping only derived alleles relative to an outgroup/ancestral sequence.

[pseudoPhasingHetero.py](pseudoPhasingHetero.py) phases the sequences by random split of heterozygous sites.

[MAFtoTAB.py](MAFtoTAB.py) transforms the [MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) file to tab file

[remove_masked_intervals.py](remove_masked_intervals.py) compares a BED interval file with the BED file of masked regions and removes them.

[removeMonomorphic.py](removeMonomorphic.py) removes monomorphic positions, i.e. keeps only SNPs.

[select_genes_by_intervals.py](select_genes_by_intervals.py) extracts gene names from a bed file by provided coordinates.

[select_intervals.py](select_intervals.py) extracts lines from a calls file according to scaffold name, start and end positions.

[selectSamples.py](selectSamples.py) subsamples a genotype calls file by sample names. It also can be used to rearrange samples in a calls file.

[skip_intervals.py](skip_intervals.py) extracts lines from a calls file according to scaffold name, start and end positions.

[slidingWindowSNPs.py](slidingWindowSNPs.py) cuts genotype calls file with the given window size and outputs FASTA files for every window.

[summarySIFT.awk](summarySIFT.awk) summarizes the extracted SIFT4G annotation (output of [extractSIFT4Gannotation.py](extractSIFT4Gannotation.py))

[summarizeTAB.awk](summarizeTAB.awk) summirized the genotyope file by counting homozygot, heterozygot, missing etc.

[tabToBED.py](tabToBED.py) converts a tab-delimited file to a bed file.

[vcf_to_SIFT4G.py](vcf_to_SIFT4G.py) converts a VCF file to SIFT4G input.

[vcfTab_to_callsTab.py](vcfTab_to_callsTab.py) converts the two-character coded table produced with VariantsToTable (GATK) to the one-character coded genotype table (calls format).
