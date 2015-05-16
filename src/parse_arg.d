module parse_arg;

import std.array : split;
import std.c.stdlib : exit;
import std.conv : to, ConvException;
import std.exception : enforce;
import std.getopt;
import std.math : isNaN;
import std.stdio : stderr, writeln;

static immutable string helpString = "
largeQvalue(1)                                                                               Statistical genetics                                                                              largeQvalue(1)



NAME
       largeQvalue: A program for calculating FDR estimates with large datasets.

SYNOPSIS
       largeQvalue [options]


DESCRIPTION
       This  is  an  implementation  of the qvalue package (Alan Dabney, John D. Storey and with assistance from Gregory R. Warnes (). qvalue: Q-value estimation for false discovery rate control. R package
       version 1.34.0.) which is designed for use with large datasets where memory or computation time may be an issue with R. It has been used to analyse full cis scans of gene expression data, with  hun-
       dreds of millions of P values. A description of the algorithms and instructions for usage can be found in the accompanying paper: http://biorxiv.org/content/early/2014/10/06/010074.


OPTIONS
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



largeQvalue-1.0.1                                                                              27th March 2015                                                                                 largeQvalue(1)
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
  string sep = "t";
  double lambdaStart;
  double lambdaEnd;
  double lambdaStep;
  double df = 3;
  size_t col = 1;
  size_t seed;
  string input = "";
  string param = "";
  string outF = "";
  double fast = 2.0;

  this(string[] args)
  {

    try
    {
      getopt(args, "help", &help, "version", &version_, "header", &header,
        "boot", &boot, "log", &logSmooth, "robust", &robust, "issorted",
        &issorted, "pi0", &pi0, "lambda", &lambda, "sep", &sep, "df", &df,
        "col", &col, "seed", &seed, "input", &input, "param", &param,
        "out", &outF, "fast", &fast);
    }
    catch (Exception e)
    {
      stderr.writeln("Failed to run. ", e.msg);
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
        enforce(lambdaStart >= 0 && lambdaEnd < 1,
          new InputException("Lambda values must lie within [0, 1) interval"));
        enforce(df >= 1 && df < (lambdaEnd - lambdaStart) / lambdaStep,
          new InputException("df must be between 1 and length of lambda"));
      }
    }
    catch (ConvException e)
    {
      stderr.writeln("Non-numeric parameters handed to lambda");
      exit(0);
    }
    catch (InputException e)
    {
      stderr.writeln(e.msg);
      exit(0);
    }

    if (fast != 2)
    {
      if (fast < 0 || fast >= 1)
      {
        stderr.writeln("Requested nominal P value threshold is not in [0, 1) interval.");
        exit(0);
      }
      col = 10;
      sep = "s";
    }

    if (sep == "s" || sep == "space")
      sep = " ";
    else
    {
      if (sep != "t" && sep != "tab")
        stderr.writeln("--sep misspecified, defaulting to tab.");
      sep = "\t";
    }

    try
    {
      enforce(pi0.isNaN || (pi0 > 0 && pi0 <= 1),
        new InputException("pi0 must be in (0, 1] interval"));
    }
    catch (InputException e)
    {
      stderr.writeln(e.msg);
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
