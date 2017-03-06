import core.stdc.stdlib : exit;
import std.algorithm : canFind, makeIndex, min, reduce, reverse;
import std.array : array, join, split;
import std.conv : ConvException, to;
import std.exception : enforce;
import std.math : exp, fabs, fmin, isNaN, log, pow;
import std.range : assumeSorted, chunks, indexed, iota, zip;
import std.stdio : File, stderr, stdin, stdout, tmpfile, write;
import std.string : chomp;

import parse_arg;
import pi0_calc : binPVals, getBootPi0, getSmootherPi0;

pure nothrow double[] pValtoQ(in double[] pVal, ref size_t[] orderIndex, double pi0, bool robust)
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
      {
        qVal[j] = robust ? pi0 * pVal[j] * len / ((i + 1) * (1 - pow(1 - pVal[j], len)))
          : pi0 * pVal[j] * len / (i + 1);
      }
      dupcount = 0;
    }
  }

  return qVal;
}

void all_p_values(Opts opts)
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

  double[] pVals;
  double pVal;

  if (opts.header && !opts.getPi)
  {
    outFile.writeln(chomp(inFile.readln), "\tQvalue");
  }

  size_t counter = 0;

  foreach (ref line; inFile.byLine)
  {
    counter++;

    if (tmp && !opts.getPi)
    {
      tmpFile.writeln(line);
    }

    auto splitLine = line.split;

    try
    {
      enforce(splitLine.length > opts.col,
          new InputException("Requested column " ~ to!string(opts.col + 1) ~ ", but row " ~ to!string(
            counter) ~ " has only " ~ splitLine.length.to!string ~ " columns."));
      pVal = to!double(splitLine[opts.col]);
      if (!pVal.isNaN)
      {
        pVals ~= pVal;
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

  try
  {
    enforce(pVals.length > 0, new InputException("No p values read"));
    foreach (ref e; pVals)
    {
      enforce(e <= 1 && e >= 0, new InputException("Some p values are outside [0, 1] interval"));
    }
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
    {
      e = i;
    }
  }
  else
  {
    makeIndex!("a<b")(pVals, orderIndex);
  }

  double pi0Final;

  auto pi0Count = binPVals(pVals, orderIndex, opts);

  if (opts.boot)
  {
    pi0Final = getBootPi0(opts, pi0Count, pVals.length, paramFile);
  }
  else if (opts.pi0.isNaN)
  {
    pi0Final = getSmootherPi0(opts, pi0Count, pVals.length, paramFile);
  }
  else
  {
    pi0Final = opts.pi0;
    if (opts.writeParam)
    {
      paramFile.writeln("Using specified value of π₀ = ", pi0Final);
    }
  }

  if (opts.getPi)
  {
    outFile.writeln(pi0Final);
  }
  else
  {
    double[] qVal = pValtoQ(pVals, orderIndex, pi0Final, opts.robust);

    reverse(orderIndex);

    if (qVal[orderIndex[0]] > 1)
    {
      qVal[orderIndex[0]] = 1;
    }
    foreach (ref e; zip(orderIndex[0 .. ($ - 1)], orderIndex[1 .. $]))
    {
      qVal[e[1]] = min(qVal[e[1]], qVal[e[0]], 1);
    }

    bool fastQTLData = opts.fast == 2.0 ? false : true;
    double nomThreshold;

    if (fastQTLData)
    {
      import std.algorithm : countUntil;

      size_t firstThreshold = orderIndex.countUntil!(a => qVal[a] < opts.fast);
      if (firstThreshold == -1)
      {
        nomThreshold = 2;
      }
      else if (firstThreshold == 0)
      {
        nomThreshold = pVals[orderIndex[0]];
      }
      else
      {
        nomThreshold = pVals[orderIndex[firstThreshold - 1]];
      }
    }

    string noPvalue;
    string nanPvalue;

    import std.mathspecial : betaIncompleteInverse;

    if (tmp)
    {
      tmpFile.seek(0);
      auto sep = tmpFile.readln.canFind("\t") ? "\t" : " ";
      if (fastQTLData)
      {
        noPvalue = sep ~ "NA" ~ sep ~ "NA";
        nanPvalue = sep ~ "NaN" ~ sep ~ "NaN";
      }
      else
      {
        noPvalue = sep ~ "NA";
        nanPvalue = sep ~ "NaN";
      }
      tmpFile.seek(0);
      size_t i = 0;
      foreach (ref line; tmpFile.byLine)
      {
        auto splitLine = line.split;
        try
        {
          pVal = to!double(splitLine[opts.col]);
          if (!pVal.isNaN)
          {
            outFile.write(line, sep, qVal[i]);
            i++;
            if (fastQTLData)
            {
              if (nomThreshold != 2)
              {
                outFile.write(sep, betaIncompleteInverse(splitLine[2].to!double,
                    splitLine[3].to!double, nomThreshold));
              }
              else
              {
                outFile.write(sep, "NA");
              }
            }
            outFile.writeln();
          }
          else
          {
            outFile.writeln(line, nanPvalue);
          }
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
      auto sep = inFile.readln.canFind("\t") ? "\t" : " ";
      if (fastQTLData)
      {
        noPvalue = sep ~ "NA" ~ sep ~ "NA";
        nanPvalue = sep ~ "NaN" ~ sep ~ "NaN";
      }
      else
      {
        noPvalue = sep ~ "NA";
        nanPvalue = sep ~ "NaN";
      }
      inFile.seek(0);

      if (opts.header)
      {
        inFile.readln;
      }
      size_t i = 0;
      foreach (ref line; inFile.byLine)
      {
        auto splitLine = line.split;
        try
        {
          pVal = to!double(splitLine[opts.col]);
          if (!pVal.isNaN)
          {
            outFile.write(line, sep, qVal[i]);
            i++;
            if (fastQTLData)
            {
              if (nomThreshold != 2)
              {
                outFile.write(sep, betaIncompleteInverse(splitLine[2].to!double,
                    splitLine[3].to!double, nomThreshold));
              }
              else
              {
                outFile.write(sep, "NA");
              }
            }
            outFile.writeln();
          }
          else
          {
            outFile.writeln(line, nanPvalue);
          }
        }
        catch (ConvException e)
        {
          outFile.writeln(line, noPvalue);
        }
      }
    }
  }
}
