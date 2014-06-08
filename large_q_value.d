import std.c.stdlib : exit;
import std.stdio;
import std.conv;
import std.array;
import std.algorithm : find, makeIndex, reduce;
import std.range;
import std.math : fabs;

enum double EPSILON = 0.00000001;

extern (C) {
  double spline_fit(double* lambda, double* pi0, size_t length);
}

pure size_t[] bestRank(in double[] rankArray, ref size_t[] orderIndex){
  import std.algorithm : makeIndex;

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
  if (args.length==1)
    exit(0);

  File inFile = File(args[1]);

  double[] pVals;

  inFile.readln;
  foreach(ref line; inFile.byLine)
    {
      auto splitLine = line.split;
      if (find(splitLine[$-1], "lcl_vqtl").empty)
	{
	  try{
	    pVals ~= to!double(splitLine[3]);
	  } catch (ConvException e){
	    }
	}
    }
  foreach(ref e; pVals)
    assert(e<=1 && e>=0);
  auto orderIndex = new size_t[pVals.length];
  makeIndex!("a<b")(pVals, orderIndex);

  size_t[] bestIndex =  bestRank(pVals, orderIndex);

  double[] lambda = iota(0, 0.95, 0.05).array;
  double[] pi0;

  foreach(ref e; lambda)
    pi0 ~= to!double(pVals.length - map!(a=>pVals[a])(orderIndex).assumeSorted
                                                                 .lowerBound(e)
		                                                 .length) / pVals.length / (1 - e);
  double pi0Est = spline_fit(lambda.ptr, pi0.ptr, lambda.length);

  double[] qVal;
  size_t i = 0;
  foreach(ref e; pVals)
    {
      qVal ~= pi0Est * e * pVals.length / bestIndex[i];
      ++i;
    }
  writeln(pi0Est);
  writeln(pVals[0..10]);
  writeln(qVal[0..10]);
  writeln(orderIndex[2], pVals[2], " ", bestIndex[1], pVals[1]);
  qVal[orderIndex[0]] = (qVal[orderIndex[0]] > 1) ? 1 : qVal[orderIndex[0]];
}
