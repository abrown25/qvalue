module parse_arg;

import std.getopt;
import std.conv : to, ConvException;
import std.array : split;
import std.math : isNaN;
import std.stdio;
import std.c.stdlib : exit;
import std.exception : enforce;

static immutable string helpString = "Usage: large_q_value [options]
Options:
    --help     : Print help and quit
    --header   : Input has header line
    --smoother : Smoothing spline applied to log pi0 values
    --robust   : More robust values for small p values
    --boot     : Bootstrap pi0 estimate
    --pi0      : Use value of pi0 given (useful for recreating qvalue package results)
    --lambda   : Either a fixed number or a sequence given 0,0.9,0.05
    --param    : Print out parameter list to given file
    --out      : file to write results to (default stdout)
    --input    : file to take results from (must be specified, if not explicitly, the last parameter after all options have been parsed is used)
    --col      : column with p values (default 1)
";

class InputException : Exception {
  pure this(string s) {super(s);}
}

class Opts{

  bool help = false;
  bool header = false;
  bool writeParam = false;
  bool smoother = false;
  bool robust = false;
  bool boot = false;
  double pi0;
  string lambda = "0,0.9,0.05";
  double lambdaStart;
  double lambdaEnd;
  double lambdaStep;
  int ncoeff = 5;
  size_t col = 1;
  string input = "";
  string param = "";
  string outF = "";

  this(string[] args){
    if (args.length==0){
      writeln(helpString);
      exit(0);
    }

    try{
      getopt(
	     args,
	     "help", &help,
	     "header", &header,
	     "smoother", &smoother,
	     "robust", &robust,
	     "boot", &boot,
	     "pi0", &pi0,
	     "lambda", &lambda,
	     "coef", &ncoeff,
	     "col", &col,
	     "input", &input,
	     "param", &param,
	     "out", &outF);
    } catch (Exception e){
      writeln("Failed to run. ", e.msg);
      exit(0);
    }

    try{
      enforce(ncoeff > 3, new InputException("At least 4 coefficients required for splines"));
    } catch (InputException e){
      writeln(e.msg);
      exit(0);
    }

    try{
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
	  enforce(lambdaEnd > lambdaStart + ncoeff * lambdaStep,
	  	  new InputException("Lambda sequence too short to estimate splines"));
	  enforce(lambdaStart >= 0 && lambdaEnd < 1,
		  new InputException("Lambda values must lie within [0, 1) interval"));
	}
    } catch (ConvException e) {
      writeln("Non-numeric parameters handed to lambda");
      exit(0);
    } catch (InputException e) {
      writeln(e.msg);
      exit(0);
    }

    try {
      enforce(pi0.isNaN || (pi0 > 0 && pi0 <=1), new InputException("pi0 must be in (0, 1] interval"));
    } catch (InputException e) {
      writeln(e.msg);
      exit(0);
    }

    if (args.length != 0 && input=="")
      input = args[$ - 1];
    if (param != "")
      writeParam = true;
    if(help){
      writeln(helpString);
      exit(0);
    }
  }
}
