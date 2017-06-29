#!/usr/bin/awk

# counts number of genotypes and missing data (N) in every sample (column)

# Input
# CHROM   POS     REF     12.4.GT 13.16.GT        16.9.GT 39-12-28.GT     5.16.GT 70.5.GT AL87.GT BEL5.GT DL174.GT        FR50.GT GY37.GT HJC419.GT       HRB132.GT  HY85.GT  IRRU2.GT        JO56.GT JZH152.GT       KMB206.GT       NJ219.GT        PL.GT   RK32.GT SE14.GT SE33.GT STA4.GT STJ2.GT TBS195.GT       TR73.GT TY118.GT   VLA3.GT  WAC5.GT XN444.GT
# scaffold_1      191     A       N       N       A       N       A       N       A       C       A       A       N       A       N       N       A       N       A  NA       N       N       N       A       N       N       A       A       N       N       A       A
# scaffold_1      563     T       T       N       W       N       N       T       N       W       N       W       N       N       T       N       W       T       T  NT       T       T       T       W       T       T       N       N       T       W       N       N
# scaffold_1      647     A       C       C       M       C       C       C       C       M       C       M       C       C       C       C       M       C       C  CC       C       C       C       M       C       C       C       C       C       M       C       C
# scaffold_1      669     T       T       T       T       T       T       T       T       T       T       T       T       T       K       T       T       T       T  TT       T       T       T       T       T       T       T       T       T       T       T       T
# scaffold_1      679     C       M       M       M       M       M       M       N       M       M       M       C       M       M       M       M       M       M  MM       M       M       M       M       M       M       M       M       M       M       M       M
# scaffold_1      704     T       Y       Y       Y       Y       Y       Y       N       Y       Y       Y       Y       T       Y       Y       T       Y       Y  YY       Y       Y       Y       Y       Y       Y       Y       Y       Y       T       Y       Y
# scaffold_1      721     T       C       C       C       C       C       C       N       C       C       C       N       C       C       C       C       C       C  CC       C       C       C       C       C       C       C       C       C       C       C       C
# scaffold_1      722     C       Y       Y       Y       Y       Y       Y       N       Y       Y       Y       N       Y       Y       Y       Y       Y       Y  YY       Y       Y       Y       Y       Y       Y       Y       Y       Y       C       Y       Y
# scaffold_1      733     G       K       K       N       K       N       K       N       G       N       N       N       N       K       N       N       K       N  NN       N       K       K       N       K       K       N       K       N       G       K       N

# # or
# CHROM   POS     REF     12.4.GT 13.16.GT        16.9.GT 39-12-28.GT     5.16.GT 70.5.GT AL87.GT BEL5.GT DL174.GT        FR50.GT GY37.GT HJC419.GT       HRB132.GT  HY85.GT  IRRU2.GT        JO56.GT JZH152.GT       KMB206.GT       NJ219.GT        PL.GT   RK32.GT SE14.GT SE33.GT STA4.GT STJ2.GT TBS195.GT       TR73.GT TY118.GT   VLA3.GT  WAC5.GT XN444.GT
# scaffold_1      191     A       ./.     ./.     A/A     ./.     A/A     ./.     A/A     C/C     A/A     A/A     ./.     A/A     ./.     ./.     A/A     ./.     A/A   ./.      A/A     ./.     ./.     ./.     A/A     ./.     ./.     A/A     A/A     ./.     ./.     A/A     A/A
# scaffold_1      563     T       T/T     ./.     T/A     ./.     ./.     T/T     ./.     T/A     ./.     T/A     T/*     T/*     T/T     T/*     T/A     T/T     T/T   T/*      T/T     T/T     T/T     T/T     T/A     T/T     T/T     ./.     ./.     T/T     T/A     ./.     T/*
# scaffold_1      647     A       C/C     C/C     A/C     C/C     C/C     C/C     C/C     A/C     C/C     A/C     C/C     C/C     C/C     C/C     A/C     C/C     C/C   C/C      C/C     C/C     C/C     C/C     A/C     C/C     C/C     C/C     C/C     C/C     A/C     C/C     C/C
# scaffold_1      669     T       T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/G     T/T     T/T     T/T     T/T   T/T      T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T     T/T
# scaffold_1      679     C       C/A     C/A     C/A     C/A     C/A     C/A     ./.     C/A     C/A     C/A     C/C     C/A     C/A     C/A     C/A     C/A     C/A   C/A      C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A     C/A
# scaffold_1      704     T       T/C     T/C     T/C     T/C     T/C     T/C     ./.     T/C     T/C     T/C     T/C     T/T     T/C     T/C     T/T     T/C     T/C   T/C      T/C     T/C     T/C     T/C     T/C     T/C     T/C     T/C     T/C     T/C     T/T     T/C     T/C
# scaffold_1      721     T       C/C     C/C     C/C     C/C     C/C     C/C     ./.     C/C     C/C     C/C     ./.     C/C     C/C     C/C     C/C     C/C     C/C   C/C      C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C     C/C
# scaffold_1      722     C       C/T     C/T     C/T     C/T     C/T     C/T     ./.     C/T     C/T     C/T     ./.     C/T     C/T     C/T     C/T     C/T     C/T   C/T      C/T     C/T     C/T     C/T     C/T     C/T     C/T     C/T     C/T     C/T     C/C     C/T     C/T
# scaffold_1      733     G       G/T     G/T     G/*     G/T     ./.     G/T     ./.     G/G     G/*     G/*     ./.     G/*     G/T     G/*     G/*     G/T     G/*G/*      ./.     G/*     G/T     G/T     G/*     G/T     G/T     G/*     G/T     G/*     G/G     G/T     ./.
# 
# # Output
# REF 9 0 9 0 0
# 12.4.GT 8 4 4 1 0
# 13.16.GT 7 4 3 2 0
# 16.9.GT 8 5 3 1 0
# 39-12-28.GT 7 4 3 2 0
# 5.16.GT 7 3 4 2 0
# 70.5.GT 8 4 4 1 0
# AL87.GT 3 0 3 6 0
# BEL5.GT 9 5 4 0 0
# DL174.GT 7 3 4 2 0
# FR50.GT 8 5 3 1 0
# GY37.GT 4 1 3 5 0
# HJC419.GT 7 2 5 2 0
# HRB132.GT 8 5 3 1 0
# HY85.GT 6 3 3 3 0
# IRRU2.GT 8 4 4 1 0
# JO56.GT 8 4 4 1 0
# JZH152.GT 8 3 5 1 0
# KMB206.GT 6 3 3 3 0
# NJ219.GT 8 3 5 1 0
# PL.GT 7 3 4 2 0
# RK32.GT 8 4 4 1 0
# SE14.GT 8 4 4 1 0
# SE33.GT 8 5 3 1 0
# STA4.GT 8 4 4 1 0
# STJ2.GT 8 4 4 1 0
# TBS195.GT 7 3 4 2 0
# TR73.GT 8 4 4 1 0
# TY118.GT 7 3 4 2 0
# VLA3.GT 8 3 5 1 0
# WAC5.GT 8 4 4 1 0
# XN444.GT 7 3 4 2 0



# awk -f summarySIFT.awk input.tab

{if (NF > maxNF ) {
    for (i = 3; i <= NF; i++)
        countN[i] = 0; countHomo[i] = 0; countHetero[i] = 0; countNA[i] = 0; maxNF = NF;
    }
    if (NR == 1 ) { for (i = 3; i <= NF; i++) samples[i] = $i;}
    else {
    for (i = 3; i <= NF; i++)
        {if ($i == "N" || $i == "./.") countN[i]++;
        else if  ($i == "A" || $i == "T" || $i == "G" || $i == "C" \
            || $i == "A/A" || $i == "T/T" || $i == "G/G" || $i == "C/C") countHomo[i]++;
        else if  ($i == "R" || $i == "Y" || $i == "M" || $i == "K" || $i == "S" || $i == "W" \
            || $i == "G/A" || $i == "T/C" || $i == "A/C" || $i == "G/T" || $i == "C/G" || $i == "A/T" \
            || $i == "A/G" || $i == "C/T" || $i == "C/A" || $i == "T/G" || $i == "G/C" || $i == "T/A") countHetero[i]++;
        else countNA[i]++; # NA -> may be single position indels
        }
    }
}
    END {
        print "Sample", "Genotypes", "Heterozygots", "Homozygots", "Missing", "Unknown";
        for (i = 3; i <= maxNF; i++)
            print samples[i], countHomo[i]+countHetero[i]+0, countHetero[i]+0, countHomo[i]+0,  countN[i]+0, countNA[i]+0;
        }
