#!/bin/bash

UNDERLINE="\033[4m"
BOLD="\033[1m"
RESET="\033[0m"

usage()
{
    echo
    echo -e "${UNDERLINE}Run a single SOFI FD modeling test${RESET}"
    echo
    echo -e "${BOLD}$(basename $0) test${RESET}"
    echo
    echo "where 'test' is:"
    echo "  1      :: acoustic wave equation"
    echo "  2      :: visco-acoustic wave equation"
    echo "  3      :: elastic wave equation"
    echo "  4      :: visco-elastic wave equation"
    echo "  5      :: elastic VTI wave equation"
    echo "  6      :: visco-elastic VTI wave equation"
    echo "  7      :: elastiv TTI wave equation"
    echo "  8      :: visco-elastic TTI wave equation"
    echo "  3_fy   :: elastic wave equation with y-force source"
    echo "  8_fy   :: visco-elastic TTI wave equation with y-force source"
    echo "  4_fdt4 :: visco-elastic wave equation with 4th order in time"
    echo
}

if [[ $# -eq 0 ]] ; then 
    usage
    exit 0
fi

if [[ $# -gt 1 ]] ; then 
    echo "Error: more than one input argument given." >&2
    exit 1
fi

weq="$1"

if [[ "${weq}" =~ ^("1"|"2"|"3"|"4"|"5"|"6"|"7"|"8"|"3_fy"|"8_fy"|"4_fdt4")$ ]]; then
    echo "Running test for weq${weq}..."
else
    echo "Error: argument not in allowed list of 1, 2, 3, 4, 5, 6, 7, 8, 3_fy, 8_fy or 4_fdt4."
    exit 1
fi

mpiexe="mpirun"

mpi_found=$(which ${mpiexe} &>/dev/null ; echo $?)
if [[ ${mpi_found} -ne 0 ]] ; then
    echo "Error: ${mpiexe} program not available."
    exit 1
fi

sofiexe="../bin/sofi2D"

sofi_found=$(ls ${sofiexe} &>/dev/null ; echo $?)
if [[ ${sofi_found} -ne 0 ]] ; then
    echo "Error: s${sofiexe} program not found."
    exit 1
fi

json="weq${weq}.json"

json_found=$(ls ${json} &>/dev/null ; echo $?)
if [[ ${json_found} -ne 0 ]] ; then
    echo "Error: json file ${json} not found."
    exit 1
fi

snapexe="../bin/snapmerge"

snap_found=$(ls ${snapexe} &>/dev/null ; echo $?)
if [[ ${snap_found} -ne 0 ]] ; then
    echo "Error: ${snapexe} program not found."
    exit 1
fi

outdir=$(echo ${weq} | sed 's/_fy//g' | sed 's/_fdt4//g')

${mpiexe} -np 4 ${sofiexe} ${json} && ${snapexe} ${json} && rm -f weq${outdir}/*.bin.{vx,vy,p,curl,div}.*.*
