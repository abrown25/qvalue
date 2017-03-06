import parse_arg;
import std.conv : to;

import all_p_vals : all_p_values;
import threshold : threshold_p_values;

void main(in string[] args)
{
  Opts opts = new Opts(to!(string[])(args));
  if (opts.threshold == 1)
  {
    all_p_values(opts);
  }
  else
  {
    threshold_p_values(opts);
  }
}
