#!/bin/bash

BOLD="\033[1m"
UNDERLINE="\033[4m"
RED="\033[0;31m"
BLUE="\033[0;34m"
RESET="\033[0m"
VERSION="-1"
MINREQVERSION="2.2.12"
USEOPTIONAL=1

################################################################################

version()
{
    echo "$@" | awk -F. '{ printf("%d%03d%03d%03d\n", $1,$2,$3,$4); }'
}

################################################################################

usage()
{
    echo 
    echo -e "${UNDERLINE}Code beautifier for C using indent ${VERSION}${RESET}"
    echo
    echo -e "${BOLD}$(basename $0) [-o] file(s)${RESET}"
    echo 
    echo "   -o      : overwrite original source code file(s)"
    echo "   file(s) : source code file(s) .c or .h to re-format"
    echo 
    echo "Each file is output to the same location but with an "
    echo "additional suffix '_new'. If the -o option is given, "
    echo "the original file is backed up as '_bck' and then "
    echo "overwritten with the generated '_new' version."
    echo
    echo "Current indent parameters:"
    echo "${options} ${options212} ${typedef}" | tr -s " " | fold -sw55
    echo
    if [ ${USEOPTIONAL} -eq 0 ] ; then
	echo "Version of indent <${MINREQVERSION}. Options '-sar -as -fnc'"
	echo "are disabled on this host."
	echo
    fi
    echo "Example:"
    echo "$(basename $0) -o ../src/*.c ../src/*.h"
    echo
}

################################################################################

HAVE_INDENT=$(which indent &>/dev/null ; echo $?)

if ! [ ${HAVE_INDENT} -eq 0 ] ; then
    echo "Program 'indent' required but not found. Please install."
    exit 1
fi

VERSION=$(indent --version | awk '{print $3}')

if [ $(version ${VERSION}) -lt $(version ${MINREQVERSION}) ] ; then
    USEOPTIONAL=0
fi

options="-nbad -bap -bbb -sob -sc -br -ce -cdw -cli2 -ss -npcs \
-ncs -saf -sai -saw -nbc -di1 -nbfda -nbfde -brs -blf \
-nut -i4 -ci0 -l120 -lp -ip0 -nbbo -bli0 -c33 -cbi0 -cd0 \
-cp33 -nlps -nprs -npsl -ts4"

typedef="-T AcqVar -T WEQTYPE -T GlobVar -T log_Level -T log_Program \
-T MemModel -T MemWavefield -T SUgather -T SUhead -T FILE -T Perform"

# the following options are only available for indent version >=2.2.12
if [ ${USEOPTIONAL} -eq 0 ] ; then
    options212=""
else
    options212="-sar -as -fnc"
fi

if [[ $# -eq 0 ]] ; then usage ; exit 0 ; fi

SUFFIX_NEW="_new"
SUFFIX_BCK="_bck"
OVERWRITE=0

while getopts "o" _options; do
    case $_options in
	o ) OVERWRITE=1 ;;
        \?) echo "Invalid option." && exit 1 ;;
        * ) ;;
    esac
done
shift $(($OPTIND-1))

ARGS="$*"
if [ -z "${ARGS}" ] ; then
    echo "No input file(s) given. Specify at least one .c or .h file." >&2
    exit 1
fi

for file in ${ARGS} ; do
    echo "Working on ${file}..."
    out="${file}${SUFFIX_NEW}"
    bck="${file}${SUFFIX_BCK}"
    echo -e -n "${RED}"
    indent ${options} ${options212} ${typedef} ${file} -o ${out}
    echo -e -n "${RESET}"
    if [ ${OVERWRITE} -eq 1 ] ; then
	cp -p ${file} ${bck}
	iret=$?
	if [ ${iret} -ne 0 ] ; then
	    echo "Warning: could not create backup file; skipping overwrite."
	else
	    echo "-> Moving ${out} back to ${file}"
	    mv ${out} ${file}
	    iret=$?
	    if [ ${iret} -ne 0 ] ; then
		echo "Warning: could not move newly formatted file back to original."
	    fi
	fi
    fi
done

