// DO NOT EDIT: this file is automatically derived from usage/stfinvbase_description_usage.txt
#include "stfinvbase_description_usage.h"
char stfinvbase_description_usage[]=
{
  "Options and parameters in common for all procedures:\n"
  "  verbose       produce verbose output (if implemented)\n"
  "  DEBUG=l       produce debug output with level l\n"
  "  exp=k         apply offset dependent weights to signals\n"
  "\n"
  "Due to the amplitude decay with offset to the source signals from receivers at\n"
  "larger offset would contribute less to the optimization criterion for which\n"
  "the source wavelet correction filter is constructed. The option exp provides\n"
  "means to add a weight factor ((r/1m)**k) to each signal, where r is the\n"
  "receiver to source offset. This is used to compensate the decrease in signal\n"
  "amplitude. If the energy in the original signal decays with ((-r/1m)**(2*k))\n"
  "all traces will contribute at a similar level to the derived correction filter\n"
  "after application of the weight factors.\n"
};
