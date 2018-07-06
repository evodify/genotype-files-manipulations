# Genotype calls files manipulations

Set of scripts to manipulate tab-delimited genotype calls files as well as to convert them to other popular formats.

All python scripts contain description of input and output data format in a header of each file.
To see possible options, run python script with --help option:
`python script.py --help`

Most of these scripts require the custom python module `calls`, so make sure that you also download and put the file `calls.py` in the same directory where your scripts are.

Examples of a tab-delimited genotype calls file (hereafter, tab file).

Two-character coded table (e.g. produced with [VariantsToTable](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_variantutils_VariantsToTable.php) from the GATK) :

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

[assessNs_in_callsTab.py](assessNs_in_callsTab.py) calculates missing data (Ns) per position/sample and visualize the results.

[calculateNsPerWindow.py](calculateNsPerWindow.py) calculates number of positions with missing data (Ns) using the sliding window approach.

[calls.py](calls.py) is a custom python module. It is a dependency for the most of the scripts listed here.

[calls_to_ped_map.py](calls_to_ped_map.py) converts genotype calls file to ped and map files suitable for [PLINK](http://zzz.bwh.harvard.edu/plink/).

[calls_to_treeMix_input.py](calls_to_treeMix_input.py) outputs alleles counts file that is required as input for [TreeMix](https://bitbucket.org/nygcresearch/treemix/wiki/Home).

[callsToBED.py](callsToBED.py) converts a tab-delimited file to a bed file.

[callsToFastaPhy_RAM.py](callsToFastaPhy_RAM.py) converts genotype calls file to FASTA and PHYLIP with little RAM consumption.

[callsToFastaPhy_speed.py](callsToFastaPhy_speed.py) converts genotype calls file to FASTA and PHYLIP fast but consumes a lot of RAM.

[combine_overlapping_BEDintervals.py](combine_overlapping_BEDintervals.py) combines overlapping genetic intervals in the BED format.

[extractSIFT4Gannotation.py](extractSIFT4Gannotation.py) extracts the [SIFT4G annotation](http://sift.bii.a-star.edu.sg/sift4g/AnnotateVariants.html) for a given set of samples according to their genotypes.

[FastaToPhylip.py](FastaToPhylip.py) converts FASTA to PHYLIP.

[FastaToTab.py](FastaToTab.py)  converts FASTA to tab-delimited file with columns: Chr, Pos, REF.

[filterByNs_callsTab.py](filterByNs_callsTab.py) removes all sites that consists of more than a given amount of missing data (Ns).

[find_popSpecificAlleles_in_callsTab.py](find_popSpecificAlleles_in_callsTab.py) outputs only unique allele of one population relative to another.

[findCommonAlleles.py](findCommonAlleles.py) outputs common and rare alleles in a given set of samples.

[keep_biallelic_in_callsTab.py](keep_biallelic_in_callsTab.py) removes sites with more than two alleles.

[make_input_MSMC_from_callsTab.py](make_input_MSMC_from_callsTab.py) makes input for [MSMC](https://github.com/stschiff/msmc).

[makeSweepFinderInput_from_callsTab.py](makeSweepFinderInput_from_callsTab.py) makes an input file for [SweepFinder](http://people.binf.ku.dk/rasmus/webpage/sf.html).

[make_input_stairway_plot_v1_BS.py](make_input_stairway_plot_v1_BS.py) makes input files including bootstrap replicates for [Stairway](https://sites.google.com/site/jpopgen/stairway-plot) version 1.

[make_input_stairway_plot_v2.py](make_input_stairway_plot_v2.py) makes an input file for [Stairway](https://sites.google.com/site/jpopgen/stairway-plot) version 2.

[merge_phased_callsTab.py](merge_phased_callsTab.py) merges phased sites into two-character coded genotype file.

[merge_SNP_wholeGenome_TabFiles.py](merge_SNP_wholeGenome_TabFiles.py) merges whole genome and SNPs tab files. This is needed because non-polymorphic sites and SNPs are filtered differently with [GATK](https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_gatk_tools_walkers_filters_VariantFiltration.php).

[mergeChrPos_in_callsTab.py](mergeChrPos_in_callsTab.py) merges all chromosomes into continuous genomic coordinates.

[mergeTabFiles.py](mergeTabFiles.py) merges two tab files by their overlapping positions.

[polarizeGT_in_callsTab.py](polarizeGT_in_callsTab.py) polarizes the genotype data by keeping only derived alleles relative to an outgroup/ancestral sequence.

[pseudoPhasingHetero_in_callsTab.py](pseudoPhasingHetero_in_callsTab.py) phases the sequences by random split of heterozygous sites.

[MAFtoTAB.py](MAFtoTAB.py) transforms the [MAF](https://genome.ucsc.edu/FAQ/FAQformat.html#format5) file to tab file

[remove_Insertions_from_callsTab.py](remove_Insertions_from_callsTab.py) removes insertions of longer than 1 bp and replaces deletions of 1 bp marked as "*" with "-".

[remove_masked_intervals_from_callsTab.py](remove_masked_intervals_from_callsTab.py) removes the masked sites from a tab file. The masked sites are provided in a BED file.

[remove_masked_intervals_fromBED.py](remove_masked_intervals_fromBED.py) compares a BED interval file with the BED file of masked regions and removes them.

[removeMonomorphic_in_callsTab.py](removeMonomorphic_in_callsTab.py) removes monomorphic positions, i.e. keeps only SNPs.

[select_genes_by_intervals.py](select_genes_by_intervals.py) extracts gene names from a bed file by provided coordinates.

[select_intervals_in_callsTab.py](select_intervals_in_callsTab.py) extracts lines from a calls file according to scaffold name, start and end positions.

[selectSamples_in_callsTab.py](selectSamples_in_callsTab.py) subsamples a genotype calls file by sample names. It also can be used to rearrange samples in a calls file.

[slidingWindowSNPs.py](slidingWindowSNPs.py) cuts genotype calls file with the given window size and outputs FASTA files for every window.

[split_calls_by_chromosomes.py](split_calls_by_chromosomes.py) splits a calls file into several files by chromosomes.

[summarySIFT.awk](summarySIFT.awk) summarizes the extracted SIFT4G annotation (output of [extractSIFT4Gannotation.py](extractSIFT4Gannotation.py))

[summarizeTAB.awk](summarizeTAB.awk) summarizes the genotyope file by counting homozygot, heterozygot, missing etc.

[vcf_to_SIFT4G.py](vcf_to_SIFT4G.py) converts a VCF file to SIFT4G input.

[vcfTab_to_callsTab.py](vcfTab_to_callsTab.py) converts the two-character coded table produced with VariantsToTable (GATK) to the one-character coded genotype table (calls format).
