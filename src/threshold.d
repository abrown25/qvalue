import parse_arg : InputException, Opts;
import pi0_calc : getBootPi0, getSmootherPi0;

import core.stdc.stdlib : exit;
import std.array : array, split;
import std.algorithm : canFind, cumulativeFold, makeIndex, min, reverse;
import std.conv : ConvException, to;
import std.exception : enforce;
import std.math : fabs, isNaN, pow;
import std.range : zip;
import std.stdio : File, stderr, stdin, stdout, writeln;
import std.string : chomp;

struct SignificantLine
{
  string line;
  double pVal;
  double qVal;

  this(string resultsLine, double Pval)
  {
    line = resultsLine;
    pVal = Pval;
  }
}

void pValtoQ(ref SignificantLine[] pVals, ref size_t[] orderIndex, double pi0,
    bool robust, size_t len)
{
  size_t dupcount = 0;

  foreach (i, ref e; orderIndex)
  {
    dupcount++;
    if (i == (len - 1) || fabs(pVals[e].pVal) < fabs(pVals[orderIndex[i + 1]].pVal))
    {
      foreach (ref j; orderIndex[(i - dupcount + 1) .. (i + 1)])
      {
        pVals[j].qVal = robust ? pi0 * pVals[j].pVal * len / ((i + 1) * (1 - pow(1 - pVals[j].pVal,
            len))) : pi0 * pVals[j].pVal * len / (i + 1);
      }
      dupcount = 0;
    }
  }
}

void threshold_p_values(Opts opts)
{
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
    {
      inFile = File(opts.input);
    }
    if (opts.outF == "")
    {
      outFile = stdout;
    }
    else
    {
      outFile = File(opts.outF, "w");
    }
    if (opts.writeParam)
    {
      paramFile = File(opts.param, "w");
    }
  }
  catch (Exception e)
  {
    stderr.writeln(e.msg);
    exit(0);
  }

  SignificantLine[] storePVals;
  double pVal;

  if (opts.header && !opts.getPi)
  {
    outFile.writeln(chomp(inFile.readln), "\tQvalue");
  }
  size_t counter;

  size_t[] pValCounts = new size_t[](cast(int)((opts.lambdaEnd - opts.lambdaStart) / opts
      .lambdaStep));

  foreach (ref line; inFile.byLine)
  {
    auto splitLine = line.split;

    try
    {
      enforce(splitLine.length > opts.col,
          new InputException("Requested column " ~ to!string(opts.col + 1) ~ ", but row " ~ to!string(
            counter) ~ " has only " ~ splitLine.length.to!string ~ " columns."));
      pVal = to!double(splitLine[opts.col]);
      if (pVal >= 0 && pVal <= 1)
      {
        counter++;
        if (pVal < opts.lambdaStart)
        {
          pValCounts[0]++;
        }
        else if (pVal <= opts.lambdaEnd)
        {
          pValCounts[cast(int)((pVal - opts.lambdaStart) / opts.lambdaStep + 1)]++;
        }
        if (pVal < opts.threshold)
        {
          storePVals ~= SignificantLine(line.idup, pVal);
        }
      }
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

  pValCounts = pValCounts.cumulativeFold!((a, b) => a + b).array;

  auto orderIndex = new size_t[storePVals.length];

  double pi0Final;

  if (opts.boot)
  {
    pi0Final = getBootPi0(opts, pValCounts, counter, paramFile);
  }
  else if (opts.pi0.isNaN)
  {
    pi0Final = getSmootherPi0(opts, pValCounts, counter, paramFile);
  }
  else
  {
    pi0Final = opts.pi0;
    if (opts.writeParam)
    {
      paramFile.writeln("Using specified value of π₀ = ", pi0Final);
    }
  }

  makeIndex!((a, b) => a.pVal < b.pVal)(storePVals, orderIndex);

  pValtoQ(storePVals, orderIndex, pi0Final, opts.robust, counter);

  reverse(orderIndex);

  if (storePVals[orderIndex[0]].qVal > 1)
  {
    storePVals[orderIndex[0]].qVal = 1;
  }
  foreach (ref e; zip(orderIndex[0 .. ($ - 1)], orderIndex[1 .. $]))
  {
    storePVals[e[1]].qVal = min(storePVals[e[1]].qVal, storePVals[e[0]].qVal, 1);
  }

  if (storePVals.length > 0)
  {
    auto sep = storePVals[0].line.canFind("\t") ? "\t" : " ";
    foreach (ref e; storePVals)
    {
      outFile.writeln(e.line, sep, e.qVal.to!string);
    }
  }

}
