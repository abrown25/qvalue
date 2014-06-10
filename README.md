## large_q_value: A program for calculating FDR estimates with large datasets

This is an implementation of the qvalue package (Alan Dabney, John D. Storey and with assistance from Gregory R. Warnes (). qvalue: Q-value estimation for false discovery rate control. R package version 1.34.0.) which is designed for use with large datasets where memory or computation time may be an issue with R.

Given a whitespace separated table (specified by --input or as the last option on the command line) and a column number (--col, default = 1), the program will write a new file where the equivalent q values are in the last column. The --param flag specifies a file to write parameter values calculated to estimate these q values, this file can be copied and pasted into R to produce diagnostic plots.

This program uses a different spline fitting algorithm to that used by the qvalue package. If exact replication of qvalue results is need, then this can be produced by running the program twice. First, run the code with the --param flag specified. Paste the parameter output into R and then run the program again, this time specifying pi0 with the --pi0 flag.

### Example command:

./large_q_value --header --col 4 --out output --param parameter_file --lambda 0,0.9,0.05 --robust QTLresults.txt

The p values can be found in the 4th column of QTLresults.txt, write the results to output and the estimated parameter values to parameter_file. Use values of lambda from 0 to 0.9 in steps of 0.05 to estimate proportion of null hypotheses (standard settings in qvalue) and produce estimates of q values robust for small p values.

### Issues:

the bootstrap option gives the wrong answers as the sampling procedure is not correct.

### List of options:

Usage: large_q_value [options]
Options:

```
    --help : Print help and quit
    --header : Input has header line
    --smoother : Smoothing spline applied to log pi0 values
    --robust : More robust values for small p values
    --boot : Bootstrap pi0 estimate
    --pi0 : Use value of pi0 given (useful for recreating qvalue package results)
    --lambda : Either a fixed number or a sequence given 0,0.9,0.05 (used to estimate pi0)
    --param : Print out parameter list to given file
    --out : file to write results to (default stdout)
    --input : file to take results from (must be specified, cannot be taken from stdin. If not given with flag, then last parameter remaining after all flagged options have been parsed is used)
    --col : column with p values (default 1)
```
