// DO NOT EDIT: this file is automatically derived from usage/stfinvidentity_description_usage.txt
#include "stfinvidentity_description_usage.h"
char stfinvidentity_description_usage[]=
{
  "Procedure: Scale with amplitude factor\n"
  "--------------------------------------\n"
  "This procedure does not modify the waveform of the synthetic data. It\n"
  "convolves the signales with a discrete delta pulse so to speak. Optionally the\n"
  "synthetics are scaled with an amplitude factor such that the weighted average\n"
  "signal energy of all traces in the scaled synthetics equals that of the\n"
  "recordings. Offset dependent weights are applied in this case. The\n"
  "appropriately scaled discrete delta pulse is returned as correction filter.\n"
  "\n"
  "Options and parameters:\n"
  "  scaleenergy   if flag is set: scale energy\n"
};
