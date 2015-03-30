## largeQvalue: A program for calculating FDR estimates with large datasets

This is an implementation of the qvalue package (Alan Dabney, John D. Storey and with assistance from Gregory R. Warnes (). qvalue: Q-value estimation for false discovery rate control. R package version 1.34.0.) which is designed for use with large datasets where memory or computation time may be an issue with R. It has been used to analyse full cis scans of gene expression data, with hundreds of millions of P values. A description of the algorithms and instructions for usage can be found in the accompanying paper: http://biorxiv.org/content/early/2014/10/06/010074. 

Given a whitespace separated table (specified by --input or as the last option on the command line) and a column number (--col, default = 1), the program will write a new file where the equivalent q values are in the last column. The --param flag specifies a file to write parameter values used to estimate these q values, this file can be copied and pasted into R to produce diagnostic plots.

### Binaries

The latest binary can be found in the bin folder on this site. An older version can be found here:

ftp://ftp.sanger.ac.uk/pub/resources/software/largeqvalue/largeQvalue.tar.gz

### Example commands:
```
largeQvalue < QTLresults.txt
```
The file QTLresults.txt contains the p values for association in the first column, the results are written to stdout with q values controlling for multiple testing, calculated using default settings, in a new final column.
```
largeQvalue --header --col 4 --out output --param parameter_file --lambda 0,0.9,0.05 --robust QTLresults.txt
```

The p values can be found in the 4th column of QTLresults.txt, write the results to output and the estimated parameter values to parameter_file. Use values of lambda from 0 to 0.9 in steps of 0.05 to estimate proportion of null hypotheses (standard settings in qvalue) and produce estimates of q values robust for small p values.

### List of options:

```
    Usage: largeQvalue [options]

    Options:

    --help    Print help and quit.

    --version Print version and quit.

    --input CHAR
            File  containing  p  values  to analyse. This can also be specified by the last argument on the command-line after all others have been parsed. If neither are present, it is taken from the
             stdin [stdin].

    --out CHAR
             File to write results to [stdout].

    --param CHAR
             Print out parameter list to specified file.

    --header  Input has header line [FALSE].

    --col INT Column with p values [1].

    --sep CHAR
             Separator to use to separate the column with q values. Specified as either space or tab (which can be shortened to s or t) [tab].

    --issorted
             File has already been sorted with no missing values [FALSE].

    --pi0 DOUBLE
             Use given value of pi0.

    --lambda DOUBLE(,DOUBLE,DOUBLE)
             Either a fixed number or a sequence given as 0,0.9,0.05 (start,end,step) [0,0.9,0.05].

	--robust  More robust values for small p values [FALSE].

    --df DOUBLE
              Number of degrees of freedom used by the spline when estimating pi0 [3].

    --log     Smoothing spline applied to log pi0 values [FALSE].

    --boot    Apply bootstrap method to find pi0 [FALSE].

    --seed DOUBLE
              Set seed for generating bootstrap samples [0].

    --fast DOUBLE
              Report nominal P value threshold for each gene corresponding to given FDR threshold when input is a fastQTL results file.

```
