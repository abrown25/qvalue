import std.c.stdlib : exit;
import std.stdio;
import std.conv;
import std.array;
import std.algorithm : find, makeIndex, reduce;
import std.range;
import std.string : chomp;
import std.math : fabs, isNaN, log, exp, pow;
import parse_arg;

extern (C) {
  void splineFit(double* lambda, double* pi0, double* pi0Est, size_t length, int ncoeff);
}

pure double[] pValtoQ(in double[] pVal, ref size_t[] orderIndex, double pi0, bool robust){
  immutable size_t len = pVal.length;
  auto qVal = new double[len];
  size_t dupcount = 0;

  foreach(i, ref e; orderIndex)
    {
      dupcount++;
      if (i == (len - 1) || fabs(pVal[e]) < fabs(pVal[orderIndex[i + 1]]))
	{
	  foreach(ref j; orderIndex[(i - dupcount + 1)..(i + 1)])
	    if (robust)
	      qVal[j] = pi0 * pVal[j] * len / ((i + 1) * (1 - pow(1 - pVal[j], len)));
	    else
	      qVal[j] = pi0 * pVal[j] * len / (i + 1);
	  dupcount = 0;
	}
    }

  return qVal;
}

size_t lowerSearch(size_t[] order, ref double[] values, double value){
  size_t offset = 0;
  while (!order.empty)
    {
      auto i = order.length / 2;
      auto mid = values[order[i]];
      if (mid > value)
	order = order[0 .. i];
      else if (mid < value && ((i == order.length - 1) || values[order[i + 1]] >= value))
	return i + 1 + offset;
      else
	{
	  order = order[i + 1 .. $];
	  offset += i + 1;
	}
    }
  return 0;
}


void main(in string[] args){

  Opts opts = new Opts(cast(string[])args);

  File inFile;
  File outFile;
  File paramFile;

  try{
    inFile = File(opts.input);
    if (opts.outF == "")
      outFile = stdout;
    else
      outFile = File(opts.outF, "w");
    if (opts.writeParam)
      paramFile = File(opts.param, "w");
  } catch (Exception e){
    writeln(e.msg);
    exit(0);
  }


  double[] pVals;

  if (opts.header)
    inFile.readln;

  foreach(ref line; inFile.byLine)
    {
      auto splitLine = line.split;
      try{
	enforce(splitLine.length > opts.col -1,
		new InputException("Column with p values doesn't exist"));
	pVals ~= to!double(splitLine[opts.col - 1]);
      } catch (ConvException e) {
      } catch (InputException e) {
	writeln(e.msg);
	exit(0);
      }
    }

  try{
    enforce(pVals.length > 0,
	    new InputException("No p values read"));
    foreach(ref e; pVals)
      enforce(e <= 1 && e >= 0,
	      new InputException("Some p values are outside [0, 1] interval"));
  } catch (InputException e){
    writeln(e.msg);
    exit(0);
  }

  auto orderIndex = new size_t[pVals.length];
  makeIndex!("a<b")(pVals, orderIndex);

  double pi0Final;

  if (opts.pi0.isNaN)
    {
      double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
      double[] pi0;
      double[] pi0Est = new double[lambda.length];

      pi0 ~= lowerSearch(orderIndex, pVals, lambda[0]);

      foreach(ref e; 1..lambda.length)
	pi0 ~= pi0[$ - 1] + lowerSearch(orderIndex[cast(size_t)pi0[$ - 1] .. $], pVals, lambda[e]);

      foreach(i, ref e; pi0)
	e = (pVals.length - e) / (1 - lambda[i]) / pVals.length;

      if(lambda.length != 1)
	{
	  if (opts.smoother)
	    foreach(ref e; pi0)
	      e = log(e);

	  splineFit(lambda.ptr, pi0.ptr, pi0Est.ptr, lambda.length, opts.ncoeff);

	  if (opts.smoother)
	    {
	      foreach(ref e; pi0)
		e = exp(e);
	      foreach(ref e; pi0Est)
		e = exp(e);
	    }
	}
      else
	pi0Est[0] = pi0[0];

      pi0Final = min(pi0Est[$ - 1], 1);

      try{
	enforce(pi0Final > 0,
		new InputException("Pi0 estimate is <= 0"));
      } catch (InputException e) {
	writeln(e.msg);
	exit(0);
      }

      if (opts.writeParam)
	{
	  paramFile.writeln("#The estimated value of ",to!dchar(0x03C0),"0 is:         ", pi0Final, "\n");
	  paramFile.writeln(
			    "#", to!dchar(0x03BB), " values to calculate this were:      [", join(to!(string[])(lambda), ", "), "]\n\n",
			    "#with the corresponding ", to!dchar(0x03C0), "0 values:     [", join(to!(string[])(pi0), ", "), "]\n\n",
			    "#and spline-smoothed ", to!dchar(0x03C0), "0 values:        [", join(to!(string[])(pi0Est), ", "), "]\n");
	  paramFile.writeln("###R code to produce diagnostic plots and qvalue package estimate of ", to!dchar(0x03C0),"0\n");
	  paramFile.writeln("data <- data.frame(lambda = c(", join(to!(string[])(lambda), ", "), "),
                   pi0 = c(", join(to!(string[])(pi0), ", "), "),
                   pi0Est = c(", join(to!(string[])(pi0Est), ", "), "))\n");
	  paramFile.writeln("qvalEst = smooth.spline(data$lambda, data$pi0, df = 3)$y ## replace 3 if different degrees of freedom is required
print(paste(c(\"Estimate of pi0 from qvalue package is\", qvalEst[length(qvalEst)]), collapse = ' '))\n");
	  paramFile.writeln("### Code to draw diagnostic plots with ggplot2\n");

	  paramFile.writeln("library(ggplot2)
ggplot(data = data, aes(x = lambda, y = pi0)) + geom_point() +
                                                geom_line(aes(x = lambda, y = pi0Est)) +
                                                geom_abline(slope = 0, intercept = data$pi0Est[nrow(data)], col = 'red')
");
	}
    }
  else
    {
      pi0Final = opts.pi0;
      if (opts.writeParam)
	paramFile.writeln("Using specified value of ", to!dchar(0x03C0), "0 = ", pi0Final);
    }

  double[] qVal = pValtoQ(pVals, orderIndex, pi0Final, opts.robust);

  qVal[orderIndex[$-1]] = qVal[orderIndex[$ - 1]] > 1 ? 1 : qVal[orderIndex[$ - 1]];

  foreach(ref e; iota(qVal.length - 2, 0, -1))
    qVal[orderIndex[e]] = qVal[orderIndex[e]] > qVal[orderIndex[e + 1]] ? qVal[orderIndex[e + 1]]
                                                                        : qVal[orderIndex[e]];

  inFile.seek(0);

  if (opts.header)
    outFile.writeln(chomp(inFile.readln), "\tQvalue");

  size_t i = 0;
  foreach(ref line; inFile.byLine)
    {
      auto splitLine = line.split;
      try{
	to!double(splitLine[opts.col - 1]);
	outFile.writeln(line, "\t", qVal[i]);
	i++;
      } catch (ConvException e){
	outFile.writeln(line, "\tNA");
      }
    }
}
