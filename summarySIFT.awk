#!/bin/sh

# The categories are hard-coded. You may need to replace them.

# Input:

# #CHROM  POS     sample1 sample2 sample3
# scaffold_1      2068    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED
# scaffold_1      2421    NONSYNONYMOUS|DELETERIOUS       SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED
# scaffold_1      2439    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED
# scaffold_1      2472    SYNONYMOUS|TOLERATED    NONSYNONYMOUS|DELETERIOUS       SYNONYMOUS|TOLERATED
# scaffold_1      2475    SYNONYMOUS|TOLERATED    NONSYNONYMOUS|DELETERIOUS       SYNONYMOUS|TOLERATED
# scaffold_1      2488    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED
# scaffold_1      2502    NA      SYNONYMOUS|TOLERATED    SYNONYMOUS|DELETERIOUS
# scaffold_1      2511    SYNONYMOUS|TOLERATED    NA      NA
# scaffold_1      2527    NONSYNONYMOUS|TOLERATED SYNONYMOUS|TOLERATED    SYNONYMOUS|TOLERATED

# # Output:

# Sample SYNONYMOUS|TOLERATED SYNONYMOUS|DELETERIOUS NONSYNONYMOUS|TOLERATED NONSYNONYMOUS|DELETERIOUS NA
# 12.4_Co 6 0 1 1 1
# 12.4_Cg 6 0 0 2 1
# 13.16_Co 7 1 0 0 1

# command:

# awk -f summarySIFT.awk input.SIFT


{if (NF > maxNF ) {
    for (i = 3; i <= NF; i++)
        countST[i] = 0; countSD[i] = 0; countNT[i] = 0; countND[i] = 0; countNA[i] = 0; maxNF = NF;
    }
    if (NR == 1 ) { for (i = 3; i <= NF; i++) samples[i] = $i;}
    else {
    for (i = 3; i <= NF; i++)
        {if ($i == "SYNONYMOUS|TOLERATED") countST[i]++;
        else if  ($i == "SYNONYMOUS|DELETERIOUS") countSD[i]++;
        else if  ($i == "NONSYNONYMOUS|TOLERATED") countNT[i]++;
        else if  ($i == "NONSYNONYMOUS|DELETERIOUS") countND[i]++;
        else countNA[i]++;
        }
    }
}
    END {
        print "Sample", "SYNONYMOUS|TOLERATED", "SYNONYMOUS|DELETERIOUS", "NONSYNONYMOUS|TOLERATED", "NONSYNONYMOUS|DELETERIOUS", "NA";
        for (i = 3; i <= maxNF; i++)
            print samples[i], countST[i]+0, countSD[i]+0, countNT[i]+0, countND[i]+0, countNA[i]+0;
        }
