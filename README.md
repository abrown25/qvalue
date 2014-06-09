Usage: large_q_value [options]
Options:
    --help     : Print help and quit
    --header   : Input has header line
    --smoother : Smoothing spline applied to log pi0 values
    --robust   : More robust values for small p values
    --pi0      : Use value of pi0 given (useful for recreating qvalue package results)
    --lambda   : Either a fixed number or a sequence given 0,0.9,0.05
    --param    : Print out parameter list to given file
    --out      : file to write results to (default stdout)
    --input    : file to take results from (must be specified, if not explicitly, the last parameter after all options have been parsed is used)
    --col      : column with p values (default 1)
