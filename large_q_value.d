import std.c.stdlib : exit;
import std.stdio;
import std.conv;
import std.array;
import std.algorithm : find, makeIndex, reduce;
import std.range;
import std.string : chomp;
import std.math : fabs, isNaN;
import parse_arg;

enum double EPSILON = 0.00000001;

extern (C) {
  void spline_fit(double* lambda, double* pi0, double* pi0Est, size_t length, int ncoeff);
}

pure size_t[] bestRank(in double[] rankArray, ref size_t[] orderIndex){
  immutable size_t len = rankArray.length;
  auto bestIndex = new size_t[len];
  size_t dupcount = 0;

  foreach(i, ref e; orderIndex)
    {
      dupcount++;
      if (i == (len - 1) || fabs(rankArray[e]) < fabs(rankArray[orderIndex[i + 1]]))
	{
	  foreach(ref j; orderIndex[(i - dupcount + 1)..(i + 1)])
	    bestIndex[j] = i + 1;
	  dupcount = 0;
	}
    }
  return bestIndex;
}

void main(in string[] args){

  Opts opts = new Opts(cast(string[])args[1..$]);

  File inFile = File(opts.input);

  File outFile;
  if (opts.outF == "")
    outFile = stdout;
  else
    outFile = File(opts.outF, "w");

  File paramFile;
  if (opts.writeParam)
    paramFile = File(opts.param, "w");

  double[] pVals;

  if (opts.header)
    inFile.readln;

  foreach(ref line; inFile.byLine)
    {
      auto splitLine = line.split;
	{
	  try{
	    pVals ~= to!double(splitLine[opts.col - 1]);
	  } catch (ConvException e){
	  }
	}
    }

  foreach(ref e; pVals)
    assert(e<=1 && e>=0);

  auto orderIndex = new size_t[pVals.length];
  makeIndex!("a<b")(pVals, orderIndex);

  size_t[] bestIndex =  bestRank(pVals, orderIndex);

  double pi0Final;
  if (opts.pi0.isNaN)
    {
      double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
      double[] pi0;
      pi0 ~= 0.0;

      foreach(ref e; 1..lambda.length)
	pi0 ~= pi0[$ -1 ] + map!(a=>pVals[a])(orderIndex[cast(size_t)pi0[$-1]..$])
	  .assumeSorted
	  .lowerBound(lambda[e])
	  .length;
      foreach(i, ref e; pi0)
	e = (pVals.length - e) / (1 - lambda[i]) / pVals.length;
      double[] pi0Est = new double[lambda.length];

      spline_fit(lambda.ptr, pi0.ptr, pi0Est.ptr, lambda.length, opts.ncoeff);
      pi0Final = pi0Est[$ - 1];
      if (opts.writeParam)
	{
	  paramFile.writeln("The estimated value of ",to!dchar(0x03C0),"0 is:         ", pi0Final, "\n");
	  paramFile.writeln(
			    to!dchar(0x03BB), " values to calculate this were:      [", join(to!(string[])(lambda), ", "), "]\n\n",
			    "with the corresponding ", to!dchar(0x03C0), "0 values:     [", join(to!(string[])(pi0), ", "), "]\n\n",
			    "with spline-smoothed", to!dchar(0x03C0), "0 values:     [", join(to!(string[])(pi0Est), ", "), "]");
	  paramFile.writeln("\n###R code to produce diagnostic plots and qvalue package estimate of ", to!dchar(0x03C0),"0\n");
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
    pi0Final = opts.pi0;

  double[] qVal;
  foreach(i, ref e; pVals)
    qVal ~= pi0Final * e * pVals.length / bestIndex[i];

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
