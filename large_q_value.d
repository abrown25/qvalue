import std.c.stdlib : exit;
import std.stdio;
import std.conv;
import std.array;
import std.algorithm : find, makeIndex, reduce;
import std.range;
import std.string : chomp;
import std.math : fabs, isNaN, log, exp, pow;
import parse_arg;

enum double EPSILON = 0.00000001;

extern (C) {
  void splineFit(double* lambda, double* pi0, double* pi0Est, size_t length, int ncoeff);
}

extern (C) {
  void bootSample(size_t* boots, size_t* counts, size_t countSize);
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

double[] countPi0(ref size_t[] counts, ref double[] lambda, size_t total){
  double[] pi0;
  counts[$ - 1] = total - counts[$ - 1];
  for(int i = cast(int)(counts.length - 2); i > -1; i--)
    counts[i] = counts[i + 1] - counts[i];
  foreach(i, ref e; counts)
    pi0 ~= (total - e) / (1 - lambda[i]) / total;
  return pi0;
}

double[] bootError(double[] pi0Real, double[] pi0Boot){
  double  minP = 1;
  foreach(ref e; pi0Real)
    minP = minP > e ? e : minP;
  writeln(minP);
  double[] mse = new double[pi0Real.length];
  mse[] = 0.0;

  foreach(ref e; chunks(pi0Boot, pi0Real.length))
    foreach(i, ref f; e)
      mse[i] += pow(f - minP, 2);

  return mse;
}

void main(in string[] args){

  Opts opts = new Opts(cast(string[])args[1..$]);

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
      if (opts.boot)
	{
	  double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
	  size_t[] counts = new size_t[lambda.length];
	  foreach(ref e; 1..lambda.length)
	    counts[e - 1] = map!(a => pVals[a])(orderIndex)
	      .assumeSorted
	      .upperBound(lambda[e])
	      .length;
	  foreach(ref e; iota(lambda.length - 1, 0, -1))
	    counts[e] = counts[e - 1] - counts[e];
	  counts[0] = pVals.length - counts[0];
	  size_t[] boots = new size_t[lambda.length * 100];
	  bootSample(boots.ptr, counts.ptr, counts.length);
	  double[] pi0Real = countPi0(counts, lambda, pVals.length);
	  double[] pi0Boots;
	  foreach(ref e; chunks(boots, pi0Real.length))
	    pi0Boots ~= countPi0(e, lambda, pVals.length);
	  double[] mse = bootError(pi0Real, pi0Boots);
	  double minMSE = mse[$ - 1];
	  size_t placeMin = mse.length - 1;
	  for(int i = cast(int)(mse.length - 2); i > -1 ; i--)
	    if (mse[i] < minMSE)
	      {
		minMSE = mse[i];
		placeMin = i;
	      }
	  pi0Final = pi0Real[placeMin];
	}
      else
	{
	  double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
	  double[] pi0;
	  double[] pi0Est = new double[lambda.length];

	  pi0 ~= map!(a => pVals[a])(orderIndex)
	    .assumeSorted
	    .lowerBound(lambda[0])
	    .length;

	  foreach(ref e; 1..lambda.length)
	    pi0 ~= pi0[$ - 1] + map!(a => pVals[a])(orderIndex[cast(size_t)pi0[$ - 1]..$])
	      .assumeSorted
	      .lowerBound(lambda[e])
	      .length;

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
    }
      else
	{
	  pi0Final = opts.pi0;
	  if (opts.writeParam)
	    paramFile.writeln("Using specified value of ", to!dchar(0x03C0), "0 = ", pi0Final);
	}

  double[] qVal;

  size_t[] bestIndex =  bestRank(pVals, orderIndex);

  if (opts.robust)
    qVal = map!(a => pi0Final * pVals[a] * pVals.length / (bestIndex[a] * (1 - pow(1 - pVals[a], pVals.length))))(iota(0, pVals.length)).array;
  else
    qVal = map!(a => pi0Final * pVals[a] * pVals.length / bestIndex[a])(iota(0, pVals.length)).array;

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
