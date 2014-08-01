## large_q_value: A program for calculating FDR estimates with large datasets

This is an implementation of the qvalue package (Alan Dabney, John D. Storey and with assistance from Gregory R. Warnes (). qvalue: Q-value estimation for false discovery rate control. R package version 1.34.0.) which is designed for use with large datasets where memory or computation time may be an issue with R.

Given a whitespace separated table (specified by --input or as the last option on the command line) and a column number (--col, default = 1), the program will write a new file where the equivalent q values are in the last column. The --param flag specifies a file to write parameter values used to estimate these q values, this file can be copied and pasted into R to produce diagnostic plots.

This program uses a different spline fitting algorithm to that used by the qvalue package. If exact replication of qvalue results is needed, then this can be produced by running the program twice. First, run the code with the --param flag specified. Paste the parameter output into R and then run the program again, this time specifying pi0 with the --pi0 flag.

### Binaries

Latest binaries for 64bit linux can be found here:

https://www.dropbox.com/sh/57gis7mthuc5l2k/AADMIu6MhTzhGFCOWewr-YtIa/large_q_value (version compiled with reference compiler, DMD)

https://www.dropbox.com/sh/57gis7mthuc5l2k/AAA1QqYIlGEqogbIrhTpN-Pwa/large_q_value_ldc (version compile with ldc, runs much faster)

### Example command:

./large_q_value --header --col 4 --out output --param parameter_file --lambda 0,0.9,0.05 --robust QTLresults.txt

The p values can be found in the 4th column of QTLresults.txt, write the results to output and the estimated parameter values to parameter_file. Use values of lambda from 0 to 0.9 in steps of 0.05 to estimate proportion of null hypotheses (standard settings in qvalue) and produce estimates of q values robust for small p values.

### List of options:

Usage: large_q_value [options]
Options:

```
    --help     : Print help and quit
    --version  : Print version and quit
    --header   : Input has header line (default = FALSE)
    --boot     : Apply bootstrap method to find pi0 (default = FALSE)
    --seed     : Set seed for generating bootstrap samples (default = 0, equivalent to GSL default)
    --log      : Smoothing spline applied to log pi0 values (default = FALSE)
    --robust   : More robust values for small p values (default = FALSE)
    --pi0      : Use value of pi0 given (useful for recreating qvalue package results)
    --lambda   : Either a fixed number or a sequence given 0,0.9,0.05 (default = 0,0.9,0.05)
    --param    : Print out parameter list to given file
    --out      : File to write results to (default = stdout)
    --input    : File to take results from (must be specified, if not explicitly, the last parameter after all options have been parsed is used)
    --issorted : File has already been sorted with no missing values (default = FALSE)
    --col      : Column with p values (default = 1)
```
