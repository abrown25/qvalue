import std.algorithm : assumeSorted, makeIndex, min, reduce, reverse;
import std.array : array, join;
import std.math : exp, fabs, fmin, log, pow;
import std.range : chunks, indexed, iota, zip;
import std.stdio : File, stdout;
import std.string : chomp;

import parse_arg;

extern (C) {
  void bootSample(size_t* bootCount, double* probs, size_t total, size_t countSize, size_t seed);
}

extern (C) {
  void splineFit(double* lambda, double* pi0, double* pi0Est, size_t length, int ncoeff);
}

double getBootPi0(in Opts opts, in double[] pVals, in size_t[] orderIndex, File paramFile){
  double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
  double[] pi0;
  size_t[] pi0Count;

  pi0Count ~= pVals.indexed(orderIndex)
    .assumeSorted
    .lowerBound(lambda[0])
    .length;

  foreach(ref e; 1 .. lambda.length)
    pi0Count ~= pi0Count[$ - 1] + pVals.indexed(orderIndex[pi0Count[$ - 1] .. $])
      .assumeSorted
      .lowerBound(lambda[e])
      .length;

  foreach(i, ref e; pi0Count)
    pi0 ~= (pVals.length - e) / (1 - lambda[i]) / pVals.length;
  immutable double minP = pi0.reduce!min;

  double[] probs;
  probs ~= cast(double)pi0Count[0] / pVals.length;
  foreach(e; zip(pi0Count[0..($ - 1)], pi0Count[1..$]))  
    probs ~= cast(double)(e[1] - e[0]) / pVals.length;

  size_t[] bootCount = new size_t[lambda.length * 100];
  double[] bootPi0;

  bootSample(bootCount.ptr, probs.ptr, pVals.length, lambda.length, opts.seed);
  foreach(ref e; chunks(bootCount, lambda.length))
    foreach(i, ref f; e)
      bootPi0 ~= (pVals.length - f) / (1 - lambda[i]) / pVals.length;

  double[] mse = new double[lambda.length];
  mse[] = 0.0;

  foreach(ref e; chunks(bootPi0, lambda.length))
    foreach(i, ref f; e)
      mse[i] += pow(f - minP, 2);

  double minMSE = mse[0];
  double pi0Final = pi0[0];

  foreach(e; zip(mse[1..$], pi0[1..$]))
    {
      if (e[0] < minMSE)
	{
	  minMSE = e[0];
	  pi0Final = e[1];
	}
      else if (e[0] == minMSE && e[1] < pi0Final)
	pi0Final = e[1];
    }

  if (opts.writeParam)
    {
      paramFile.writeln("#The estimated value of π₀ is:            ", pi0Final, "\n");
      paramFile.writeln("#λ values to calculate this were:         [", join(to!(string[])(lambda), ", "), "]\n\n",
			"#with the corresponding π₀ values:        [", join(to!(string[])(pi0), ", "), "]\n\n",
			"#and mean squared error estimates:        [", join(to!(string[])(mse), ", "), "]\n");
      paramFile.writeln("###R code to produce diagnostic plots for bootstrap estimates of π₀");
      paramFile.writeln("
plot.pi0.data <- data.frame(x = rep(c(", join(to!(string[])(lambda), ", "), "), 100),
                            y = c(", join(to!(string[])(bootPi0), ", "), "))");
      paramFile.writeln("
plot.data <- data.frame(x = c(", join(to!(string[])(lambda), ", "), "),
                        y = c(", join(to!(string[])(pi0), ", "), "),
                        mse = c(", join(to!(string[])(mse), ", "),"),
                        minpi0 = ", minP, ",
                        final = ", pi0Final, ")

library(ggplot2)
plot1 <- ggplot(plot.data, aes(x = x, y = y)) + geom_boxplot(data = plot.pi0.data, aes(x = x, y = y, group = x)) +
                                                geom_point(colour='blue') +
                                                geom_hline(yintercept = plot.data$minpi0, colour = 'blue') +
                                                geom_line(aes(x = x, y = mse), linetype = 'dashed') +
                                                geom_hline(yintercept = plot.data$final, colour = 'red') +
                                                geom_vline(xintercept = plot.data$x[plot.data$mse==min(plot.data$mse)], linetype = 'dashed') +
                                                ylim(0,1) +
                                                labs(x = expression(lambda), y = expression(pi[1])) +
                                                theme(axis.title = element_text(size = rel(2)))
print(plot1)
");
    }
  return pi0Final;
}

double getSmootherPi0(in Opts opts, in double[] pVals, in size_t[] orderIndex, File paramFile){
  double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
  double[] pi0;
  size_t[] pi0Count;
  double[] pi0Est = new double[lambda.length];

  pi0Count ~= pVals.indexed(orderIndex)
    .assumeSorted
    .lowerBound(lambda[0])
    .length;

  foreach(ref e; 1 .. lambda.length)
    pi0Count ~= pi0Count[$ - 1] + pVals.indexed(orderIndex[pi0Count[$ - 1] .. $])
      .assumeSorted
      .lowerBound(lambda[e])
      .length;

  foreach(i, ref e; pi0Count)
    pi0 ~= (pVals.length - e) / (1 - lambda[i]) / pVals.length;

  if(lambda.length != 1)
    {
      if (opts.logSmooth)
	foreach(ref e; pi0)
	  e = log(e);

      splineFit(lambda.ptr, pi0.ptr, pi0Est.ptr, lambda.length, opts.ncoeff);

      if (opts.logSmooth)
	{
	  foreach(ref e; pi0)
	    e = exp(e);
	  foreach(ref e; pi0Est)
	    e = exp(e);
	}
    }
  else
    pi0Est[0] = pi0[0];

  double pi0Final = fmin(pi0Est[$ - 1], 1);

  try{
    enforce(pi0Final > 0,
	    new InputException("Pi0 estimate is <= 0"));
  } catch (InputException e) {
    writeln(e.msg);
    exit(0);
  }

  if (opts.writeParam)
    {
      paramFile.writeln("#The estimated value of π₀ is:         ", pi0Final, "\n");
      paramFile.writeln("#λ values to calculate this were:      [", join(to!(string[])(lambda), ", "), "]\n\n",
			"#with the corresponding π₀ values:     [", join(to!(string[])(pi0), ", "), "]\n\n",
			"#and spline-smoothed π₀ values:        [", join(to!(string[])(pi0Est), ", "), "]\n");
      paramFile.writeln("###R code to produce diagnostic plots and qvalue package estimate of π₀\n
plot.data <- data.frame(lambda = c(", join(to!(string[])(lambda), ", "), "),
                        pi0 = c(", join(to!(string[])(pi0), ", "), "),
                        pi0Est = c(", join(to!(string[])(pi0Est), ", "), "))\n");
      paramFile.writeln("qvalEst = smooth.spline(plot.data$lambda, plot.data$pi0, df = 3)$y ## replace 3 if different degrees of freedom is required\n
print(paste(c(\"Estimate of pi0 from qvalue package is:\", qvalEst[length(qvalEst)]), collapse = ' '))\n");
      paramFile.writeln("### Code to draw diagnostic plots with ggplot2\n");

      paramFile.writeln("library(ggplot2)
plot1 <- ggplot(data = plot.data, aes(x = lambda, y = pi0)) + geom_point() +
                                                              geom_line(aes(x = lambda, y = pi0Est)) +
                                                              geom_abline(slope = 0, intercept = plot.data$pi0Est[nrow(plot.data)], col = 'red') +
                                                              ylim(0, 1)
print(plot1)
");
    }
  return pi0Final;
}

pure nothrow double[] pValtoQ(in double[] pVal, ref size_t[] orderIndex, double pi0, bool robust){
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

void main(in string[] args){
  if (args.length==1)
    {
      writeln(helpString);
      exit(0);
    }

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
	enforce(splitLine.length > opts.col - 1,
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
  if (opts.issorted)
    {
      foreach(i, ref e; orderIndex)
	e = i;
    }
  else
    makeIndex!("a<b")(pVals, orderIndex);

  double pi0Final;

  if (opts.boot)
    pi0Final = getBootPi0(opts, pVals, orderIndex, paramFile);
  else if (opts.pi0.isNaN)
    pi0Final = getSmootherPi0(opts, pVals, orderIndex, paramFile);
  else
    {
      pi0Final = opts.pi0;
      if (opts.writeParam)
	paramFile.writeln("Using specified value of ", to!dchar(0x03C0), "0 = ", pi0Final);
    }

  double[] qVal = pValtoQ(pVals, orderIndex, pi0Final, opts.robust);

  reverse(orderIndex);

  if (qVal[orderIndex[0]] > 1)
    qVal[orderIndex[0]] = 1;

  foreach(ref e; zip(orderIndex[0 .. ($ - 1)], orderIndex[1 .. $]))
    {
      if (qVal[e[1]] > qVal[e[0]])
	qVal[e[1]] = qVal[e[0]];
      if (qVal[e[1]] > 1)
	qVal[e[1]] = 1;
    }

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
