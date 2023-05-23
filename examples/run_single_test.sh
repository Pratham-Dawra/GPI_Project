#!/bin/bash

UNDERLINE="\033[4m"
BOLD="\033[1m"
RESET="\033[0m"

usage()
{
    echo
    echo -e "${UNDERLINE}Run a single SOFI FD modeling test${RESET}"
    echo
    echo -e "${BOLD}$(basename $0) [-h] [-v] [-e exe] [-s snm] test${RESET}"
    echo
    echo "-h       :: show help message and exit"
    echo "-v       :: run test using 'valgrind' memory checker (attention: this is slow)"
    echo "-e exe   :: use SOFI executable 'exe' rather than default '../bin/sofi2D'"
    echo "-s snm   :: use SNAPMERGE executable 'snm' rather than default '../bin/snapmerge'"
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
    echo "  3_fdt4 :: elastic wave equation with 4th order in time"
    echo "  4_fdt4 :: visco-elastic wave equation with 4th order in time"
    echo
}

if [[ $# -eq 0 ]] ; then usage ; exit 0 ; fi

valgrind=0
snapexe="../bin/snapmerge"
sofiexe="../bin/sofi2D"

while getopts "hvs:e:" _options; do
    case $_options in
        h ) usage && exit 0 ;;
        v ) valgrind=1 ;;
        e ) sofiexe="$OPTARG" ;;
        s ) snapexe="$OPTARG" ;;
        \?) echo "Invalid option. Run $(basename $0) -h for help." && exit 1 ;;
        * ) ;;
    esac
done
shift $(($OPTIND-1))
weq="$*"

if [ -z "${weq}" ] ; then
    echo "No test specified. Run $(basename $0) -h for help."
    exit 1
fi

if [[ "${weq}" =~ ^("1"|"2"|"3"|"4"|"5"|"6"|"7"|"8"|"3_fy"|"8_fy"|"3_fdt4"|"4_fdt4")$ ]]; then
    echo "Running test for weq${weq}..."
else
    echo "Error: argument not in allowed list of 1, 2, 3, 4, 5, 6, 7, 8, 3_fy, 8_fy, 3_fdt4 or 4_fdt4."
    exit 1
fi

mpiexe="mpirun"

mpi_found=$(which ${mpiexe} &>/dev/null ; echo $?)
if [[ ${mpi_found} -ne 0 ]] ; then
    echo "Error: ${mpiexe} program not available."
    exit 1
fi

is_openmpi=$(${mpiexe} --version | grep -c -i "Open MPI")
if [[ ${is_openmpi} -gt 0 ]] ; then
    mpiexe="${mpiexe} --mca orte_base_help_aggregate 0 --mca btl_base_warn_component_unused 0"
fi

sofi_found=$(ls ${sofiexe} &>/dev/null ; echo $?)
if [[ ${sofi_found} -ne 0 ]] ; then
    echo "Error: ${sofiexe} program not found."
    exit 1
fi

json="weq${weq}.json"

json_found=$(ls ${json} &>/dev/null ; echo $?)
if [[ ${json_found} -ne 0 ]] ; then
    echo "Error: json file ${json} not found."
    exit 1
fi

snap_found=$(ls ${snapexe} &>/dev/null ; echo $?)
if [[ ${snap_found} -ne 0 ]] ; then
    echo "Error: ${snapexe} program not found."
    exit 1
fi

valexe="valgrind"
valcommand=""

if [[ ${valgrind} -eq 1 ]] ; then
    val_found=$(which ${valexe} &>/dev/null ; echo $?)
    if [[ ${val_found} -ne 0 ]] ; then
	echo "Error: ${valexe} program not available."
	exit 1
    else
	valcommand="${valexe} --tool=memcheck"
    fi
fi

outdir=$(echo ${weq} | sed 's/_fy//g' | sed 's/_fdt4//g')

echo "Command: ${mpiexe} -np 4 ${valcommand} ${sofiexe} ${json}"

${mpiexe} -np 4 ${valcommand} ${sofiexe} ${json} && \
${snapexe} ${json} && \
rm -f weq${outdir}/*.bin.{vx,vy,p,curl,div}.*.*
