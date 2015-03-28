import std.algorithm : makeIndex, min, reduce, reverse;
import std.array : array, join;
import std.math : exp, fabs, fmin, isNaN, log, pow;
import std.range : assumeSorted, chunks, indexed, iota, zip;
import std.stdio : File, stderr, stdin, stdout, tmpfile, write;
import std.string : chomp;
import std.utf;
import parse_arg;

extern (C)
{
    void bootSample(size_t* bootCount, double* probs, size_t total, size_t countSize,
        size_t seed);
}

pure nothrow extern (C)
{
  double gsl_cdf_beta_Pinv(double x, double a, double b);
}

extern (C)
{
    void splineFit(double* xs, double* ys, double* knot, int n, double dofoff, double* results);
}

double getBootPi0(in Opts opts, in double[] pVals, in size_t[] orderIndex, File paramFile)
{
    double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
    double[] pi0;
    size_t[] pi0Count;

    pi0Count ~= pVals.indexed(orderIndex).assumeSorted.lowerBound(lambda[0]).length;

    foreach (ref e; 1 .. lambda.length)
        pi0Count ~= pi0Count[$ - 1] + pVals.indexed(orderIndex[pi0Count[$ - 1] .. $]).assumeSorted.lowerBound(
            lambda[e]).length;

    foreach (i, ref e; pi0Count)
        pi0 ~= (pVals.length - e) / (1 - lambda[i]) / pVals.length;
    immutable double minP = pi0.reduce!min;

    double[] probs;
    probs ~= to!double(pi0Count[0]) / pVals.length;
    foreach (e; zip(pi0Count[0 .. ($ - 1)], pi0Count[1 .. $]))
        probs ~= to!double(e[1] - e[0]) / pVals.length;

    size_t[] bootCount = new size_t[lambda.length * 100];
    double[] bootPi0;

    bootSample(bootCount.ptr, probs.ptr, pVals.length, lambda.length, opts.seed);
    foreach (ref e; chunks(bootCount, lambda.length))
        foreach (i, ref f; e)
            bootPi0 ~= (pVals.length - f) / (1 - lambda[i]) / pVals.length;

    double[] mse = new double[lambda.length];
    mse[] = 0.0;

    foreach (ref e; chunks(bootPi0, lambda.length))
        foreach (i, ref f; e)
            mse[i] += pow(f - minP, 2);

    double minMSE = mse[0];
    double pi0Final = pi0[0];

    foreach (e; zip(mse[1 .. $], pi0[1 .. $]))
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
        paramFile.writeln("#λ values to calculate this were:         [", lambda.to!(string[]).join(", "), "]\n\n",
			  "#with the corresponding π₀ values:        [", pi0.to!(string[]).join(", "), "]\n\n",
			  "#and mean squared error estimates:        [", mse.to!(string[]).join(", "), "]\n");
        paramFile.writeln(
			  "###R code to produce diagnostic plots for bootstrap estimates of π₀:\n");
        paramFile.writeln("plot.pi0.data <- data.frame(x = rep(c(", lambda.to!(string[]).join(", "), "), 100),
                            y = c(", bootPi0.to!(string[]).join(", "), "))\n");
        paramFile.writeln("plot.data <- data.frame(x = c(", lambda.to!(string[]).join(", "), "),
                        y = c(", pi0.to!(string[]).join(", "), "),
                        mse = c(", mse.to!(string[]).join(", "), "),
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
                                                labs(x = expression(lambda), y = expression(pi[0])) +
                                                theme(axis.title = element_text(size = rel(2)))
print(plot1)
");
    }
    return pi0Final;
}

double getSmootherPi0(in Opts opts, in double[] pVals, in size_t[] orderIndex, File paramFile)
{
    import std.algorithm : map;

    double[] lambda = iota(opts.lambdaStart, opts.lambdaEnd, opts.lambdaStep).array;
    double[] pi0;
    size_t[] pi0Count;
    double[] pi0Est = new double[lambda.length];

    pi0Count ~= pVals.indexed(orderIndex).assumeSorted.lowerBound(lambda[0]).length;

    foreach (ref e; 1 .. lambda.length)
        pi0Count ~= pi0Count[$ - 1] + pVals.indexed(orderIndex[pi0Count[$ - 1] .. $]).assumeSorted.lowerBound(
            lambda[e]).length;

    foreach (i, ref e; pi0Count)
        pi0 ~= (pVals.length - e) / (1 - lambda[i]) / pVals.length;

    if (lambda.length != 1)
    {
        if (opts.logSmooth)
            foreach (ref e; pi0)
                e = log(e);
        double[] xs = lambda.map!(a => a / (lambda[$ - 1] - lambda[0])).array;
        double[] knot = [0.0, 0.0, 0.0] ~ xs ~ [1.0, 1.0, 1.0];
        splineFit(xs.ptr, pi0.ptr, knot.ptr, lambda.length.to!int, opts.df, pi0Est.ptr);

        if (opts.logSmooth)
        {
            foreach (ref e; pi0)
                e = exp(e);
            foreach (ref e; pi0Est)
                e = exp(e);
        }
    }
    else
        pi0Est[0] = pi0[0];

    double pi0Final = fmin(pi0Est[$ - 1], 1);

    try
    {
        enforce(pi0Final > 0, new InputException("Pi0 estimate is <= 0"));
    }
    catch (InputException e)
    {
        stderr.writeln(e.msg);
        exit(0);
    }

    if (opts.writeParam)
    {
        paramFile.writeln("#The estimated value of π₀ is:         ", pi0Final, "\n");
        paramFile.writeln("#λ values to calculate this were:      [", lambda.to!(string[]).join(", "), "]\n\n",
			  "#with the corresponding π₀ values:     [", pi0.to!(string[]).join(", "), "]\n\n",
			  "#and spline-smoothed π₀ values:        [", pi0Est.to!(string[]).join(", "), "]\n");
        paramFile.writeln("###R code to produce diagnostic plots for spline estimates of π₀:\n
plot.data <- data.frame(lambda = c(", lambda.to!(string[]).join(", "), "),
                        pi0 = c(", pi0.to!(string[]).join(", "), "),
                        pi0Est = c(", pi0Est.to!(string[]).join(", "), "))\n");
        paramFile.writeln("library(ggplot2)
plot1 <- ggplot(data = plot.data, aes(x = lambda, y = pi0)) + geom_point() +
                                                              geom_line(aes(x = lambda, y = pi0Est)) +
                                                              geom_abline(slope = 0, intercept = plot.data$pi0Est[nrow(plot.data)], col = 'red') +
                                                              labs(x = expression(lambda), y = expression(pi[0])) +
                                                              theme(axis.title = element_text(size = rel(2))) +
                                                              ylim(0, 1)
print(plot1)
");
    }
    return pi0Final;
}

pure nothrow double[] pValtoQ(in double[] pVal, ref size_t[] orderIndex, double pi0,
    bool robust)
{
    immutable size_t len = pVal.length;
    auto qVal = new double[len];
    size_t dupcount = 0;

    foreach (i, ref e; orderIndex)
    {
        dupcount++;
        if (i == (len - 1) || fabs(pVal[e]) < fabs(pVal[orderIndex[i + 1]]))
        {
            foreach (ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
                if (robust)
                    qVal[j] = pi0 * pVal[j] * len / ((i + 1) * (1 - pow(1 - pVal[
                        j
                    ], len)));
                else
                    qVal[j] = pi0 * pVal[j] * len / (i + 1);
            dupcount = 0;
        }
    }
    return qVal;
}

void main(in string[] args)
{
    if (args.length == 1)
    {
        writeln(helpString);
        exit(0);
    }

    Opts opts = new Opts(to!(string[])(args));

    File inFile;
    File outFile;
    File paramFile;
    File tmpFile;
    bool tmp = false;

    try
    {
        if (opts.input == "")
        {
            tmp = true;
            inFile = stdin;
            tmpFile = File.tmpfile();
        }
        else
            inFile = File(opts.input);
        if (opts.outF == "")
            outFile = stdout;
        else
            outFile = File(opts.outF, "w");
        if (opts.writeParam)
            paramFile = File(opts.param, "w");
    }
    catch (Exception e)
    {
        stderr.writeln(e.msg);
        exit(0);
    }

    double[] pVals;
    double pVal;

    if (opts.header)
        outFile.writeln(chomp(inFile.readln), "\tQvalue");

    foreach (ref line; inFile.byLine)
    {
        if (tmp)
            tmpFile.writeln(line);

        auto splitLine = line.split;
        try
        {
            enforce(splitLine.length > opts.col - 1,
                new InputException("Column with p values doesn't exist"));
            pVal = to!double(splitLine[opts.col - 1]);
            if (!pVal.isNaN)
                pVals ~= pVal;
        }
        catch (ConvException e)
        {
        }
        catch (InputException e)
        {
            stderr.writeln(e.msg);
            exit(0);
        }
    }

    try
    {
        enforce(pVals.length > 0, new InputException("No p values read"));
        foreach (ref e; pVals)
            enforce(e <= 1 && e >= 0, new InputException(
                "Some p values are outside [0, 1] interval"));
    }
    catch (InputException e)
    {
        stderr.writeln(e.msg);
        exit(0);
    }

    auto orderIndex = new size_t[pVals.length];
    if (opts.issorted)
    {
        foreach (i, ref e; orderIndex)
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
            paramFile.writeln("Using specified value of π₀ = ", pi0Final);
    }

    double[] qVal = pValtoQ(pVals, orderIndex, pi0Final, opts.robust);

    reverse(orderIndex);

    if (qVal[orderIndex[0]] > 1)
        qVal[orderIndex[0]] = 1;

    foreach (ref e; zip(orderIndex[0 .. ($ - 1)], orderIndex[1 .. $]))
        qVal[e[1]] = min(qVal[e[1]], qVal[e[0]], 1);

    import std.algorithm : countUntil;
    double nomThreshold;

    if (opts.fast)
      {
	size_t firstThreshold = orderIndex.countUntil!(a => qVal[a] < 0.05);
	if (firstThreshold==-1)
	  nomThreshold = 2;
	else if (firstThreshold==0)
	  nomThreshold = pVals[orderIndex[0]];
	else
	  nomThreshold = pVals[orderIndex[firstThreshold - 1]];
      }
    
    string sep = opts.sep;
    string noPvalue;
    string nanPvalue;
    if (opts.fast)
      {
	noPvalue = sep ~ "NA" ~ sep ~ "NA";
	nanPvalue = sep ~ "NaN" ~ sep ~ "NaN";
      }
    else
      {
	noPvalue = sep ~ "NA";
	nanPvalue = sep ~ "NaN";
      }
	
    if (tmp)
    {
        tmpFile.seek(0);
        size_t i = 0;
        foreach (ref line; tmpFile.byLine)
        {
            auto splitLine = line.split;
            try
            {
                pVal = to!double(splitLine[opts.col - 1]);
                if (!pVal.isNaN)
                {
                    outFile.write(line, sep, qVal[i]);
                    i++;
		    if (opts.fast)
		      {
			if (nomThreshold!=2)
			  outFile.write(sep, gsl_cdf_beta_Pinv(nomThreshold, splitLine[2].to!double, splitLine[3].to!double));
			else
			  outFile.write(sep, "NA");
		      }
		    outFile.writeln();
                }
                else
                    outFile.writeln(line, nanPvalue);
            }
            catch (ConvException e)
            {
                outFile.writeln(line, noPvalue);
            }
        }
    }
    else
    {
        inFile.seek(0);
        if (opts.header)
            inFile.readln;
        size_t i = 0;
        foreach (ref line; inFile.byLine)
        {
            auto splitLine = line.split;
            try
            {
                pVal = to!double(splitLine[opts.col - 1]);
                if (!pVal.isNaN)
                {
                    outFile.write(line, sep, qVal[i]);
                    i++;
		    if (opts.fast)
		      {
			if (nomThreshold!=2)
			  outFile.write(sep, gsl_cdf_beta_Pinv(nomThreshold, splitLine[2].to!double, splitLine[3].to!double));
			else
			  outFile.write(sep, "NA");
		      }
		    outFile.writeln();
                }
                else
                    outFile.writeln(line, nanPvalue);
            }
            catch (ConvException e)
            {
                outFile.writeln(line, noPvalue);
            }
        }
    }
}
