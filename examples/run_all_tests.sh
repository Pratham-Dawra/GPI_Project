#!/bin/bash

# list="1 2 3 4 5 6 7 8 3_fy 8_fy"
list="3 4 5 6 7 8 3_fy 8_fy"

UNDERLINE="\033[4m"
BOLD="\033[1m"
RESET="\033[0m"

usage()
{
    echo
    echo -e "${UNDERLINE}Run all SOFI FD modeling test${RESET}"
    echo
    echo -e "${BOLD}$(basename $0) [-h] [-v] [-e exe] [-s snm]${RESET}"
    echo
    echo "-h       :: show help message and exit"
    echo "-v       :: run tests using 'valgrind' memory checker (attention: this is slow)"
    echo "-e exe   :: use SOFI executable 'exe' rather than default '../bin/sofi2D'"
    echo "-s snm   :: use SNAPMERGE executable 'snm' rather than default '../bin/snapmerge'"
    echo
}

valflag=0
sofiflag=0
snapflag=0
while getopts "hvs:e:" _options; do
    case $_options in
	h ) usage && exit 0 ;;
        v ) valflag=1 ;;
        e ) sofiflag=1; sofiexe="$OPTARG" ;;
        s ) snapflag=1; snapexe="$OPTARG" ;;
        \?) echo "Invalid option." && usage && exit 1 ;;
        * ) ;;
    esac
done
shift $(($OPTIND-1))
remain="$*"

if [ ! -z "${remain}" ] ; then
    echo "Unknown parameters specified: ${remain}"
    echo "Run $(basename $0) -h for help."
    exit 1
fi

if [[ ${valflag} -eq 1 ]] ; then
    valgrind="-v"
else
    valgrind=""
fi

if [[ ${sofiflag} -eq 1 ]] ; then
    sofi="-e ${sofiexe}"
else
    sofi=""
fi

if [[ ${snapflag} -eq 1 ]] ; then
    snap="-s ${snapexe}"
else
    snap=""
fi

for l in ${list} ; do
   ./run_single_test.sh ${valgrind} ${sofi} ${snap} ${l}
done
