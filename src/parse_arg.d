module parse_arg;

import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.getopt;
import std.math : isNaN;
import std.stdio : writeln;

static immutable string helpString = "Usage: largeQvalue [options]
Options:
    --help     : Print help and quit
    --version  : Print version and quit
    --input    : File to take results from. Can also be specified by the last argument on the command-line after all others have been parsed. If not present, is taken to be the stdin.
    --out      : File to write results to (default = stdout)
    --param    : Print out parameter list to given file
    --header   : Input has header line (default = FALSE)
    --col      : Column with p values (default = 1)
    --issorted : File has already been sorted with no missing values (default = FALSE)
    --pi0      : Use value of pi0 given
    --lambda   : Either a fixed number or a sequence given 0,0.9,0.05 (default = 0,0.9,0.05)
    --robust   : More robust values for small p values (default = FALSE)
    --df       : Number of degrees of freedom used by the spline to estimate pi0 (default = 3)
    --log      : Smoothing spline applied to log pi0 values (default = FALSE)
    --boot     : Apply bootstrap method to find pi0 (default = FALSE)
    --seed     : Set seed for generating bootstrap samples (default = 0, equivalent to GSL default)
";

static immutable string versionString = "largeQvalue: version 1.0.1";

class InputException : Exception
{
    pure this(string s)
    {
        super(s);
    }
}

class Opts
{

    bool help = false;
    bool version_ = false;
    bool header = false;
    bool writeParam = false;
    bool boot = false;
    bool logSmooth = false;
    bool robust = false;
    bool issorted = false;
    double pi0;
    string lambda = "0,0.9,0.05";
    double lambdaStart;
    double lambdaEnd;
    double lambdaStep;
    double df = 3;
    size_t col = 1;
    size_t seed;
    string input = "";
    string param = "";
    string outF = "";

    this(string[] args)
    {

        try
        {
            getopt(args, "help", &help, "version", &version_, "header",
                &header, "boot", &boot, "log", &logSmooth, "robust", &robust,
                "issorted", &issorted, "pi0", &pi0, "lambda", &lambda, "df",
                &df, "col", &col, "seed", &seed, "input", &input, "param",
                &param, "out", &outF);
        }
        catch (Exception e)
        {
            writeln("Failed to run. ", e.msg);
            exit(0);
        }

        try
        {
        }
        catch (InputException e)
        {
            writeln(e.msg);
            exit(0);
        }

        try
        {
            double[] lambdaOpts = to!(double[])(split(lambda, ","));
            lambdaStart = lambdaOpts[0];
            if (lambdaOpts.length == 1)
            {
                enforce(lambdaStart >= 0 && lambdaStart < 1,
                    new InputException("lambda must be in interval [0, 1)"));
                lambdaEnd = lambdaStart + 1;
                lambdaStep = 1;
            }
            else
            {
                lambdaStep = lambdaOpts[2];
                lambdaEnd = lambdaOpts[1] + lambdaStep;
                enforce(lambdaEnd > lambdaStart + df * lambdaStep,
                    new InputException("Lambda sequence too short to estimate splines"));
                enforce(lambdaStart >= 0 && lambdaEnd < 1, new InputException(
                    "Lambda values must lie within [0, 1) interval"));
		enforce(df >= 1 && df < (lambdaEnd - lambdaStart) / lambdaStep,
			new InputException("df must be between 1 and length of lambda"));
            }
        }
        catch (ConvException e)
        {
            writeln("Non-numeric parameters handed to lambda");
            exit(0);
        }
        catch (InputException e)
        {
            writeln(e.msg);
            exit(0);
        }

        try
        {
            enforce(pi0.isNaN || (pi0 > 0 && pi0 <= 1),
                new InputException("pi0 must be in (0, 1] interval"));
        }
        catch (InputException e)
        {
            writeln(e.msg);
            exit(0);
        }

        if (args.length > 1 && input == "")
            input = args[$ - 1];
        if (param != "")
            writeParam = true;
        if (help)
        {
            writeln(helpString);
            exit(0);
        }
        if (version_)
        {
            writeln(versionString);
            exit(0);
        }
    }
}
