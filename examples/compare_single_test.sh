#!/bin/bash

UNDERLINE="\033[4m"
BOLD="\033[1m"
RESET="\033[0m"

usage()
{
    echo
    echo -e "${UNDERLINE}Compare a single SOFI FD modeling test with reference results${RESET}"
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
    echo "  3_fdt4 :: elastic wave equation with 4th order in time"
    echo "  4_fdt4 :: visco-elastic wave equation with 4th order in time"
    echo "  3_fwi  :: FWI of elastic wave equation"
    echo "  4_fwi  :: FWI of visco-elastic wave equation"
    echo "  8_fwi  :: FWI of visco-elastic TTI wave equation"
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

if [[ "${weq}" =~ ^("1"|"2"|"3"|"4"|"5"|"6"|"7"|"8"|"3_fy"|"8_fy"|"3_fdt4"|"4_fdt4"|"3_fwi"|"4_fwi"|"8_fwi")$ ]]; then
    echo "Running comparison weq${weq}."
else
    echo "Error: argument not in allowed list of 1, 2, 3, 4, 5, 6, 7, 8, 3_fy, 8_fy, 3_fdt4, 4_fdt4, 3_fwi, 4_fwi or 8_fwi."
    exit 1
fi

fy=$(echo ${weq} | grep -c "_fy")
if [[ ${fy} -eq 0 ]] ; then
  is_fy=false
else
  is_fy=true
fi

fdt=$(echo ${weq} | grep -c "_fdt4")
if [[ ${fdt} -eq 0 ]] ; then
  is_fdt=false
else
  is_fdt=true
fi

no="$(echo ${weq} | sed 's/_fy//g' | sed 's/_fdt4//g')"
testdir="weq${no}"
refdir="${testdir}_ref"

if ! test -d ${testdir} ; then
    echo "Error: directory ${testdir} not found."
    exit 1
else
    echo "Directory with test results: ${testdir}"
fi

if ! test -d ${refdir} ; then
    echo "Error: directory ${refdir} not found."
    exit 1
else
    echo "Directory with reference results: ${refdir}"
fi

if ${is_fy} ; then
    sufiles="seis_fy_curl.su seis_fy_div.su seis_fy_p.su seis_fy_vx.su seis_fy_vy.su"
    snapfiles="snap_fy.bin.curl snap_fy.bin.div snap_fy.bin.p snap_fy.bin.vx snap_fy.bin.vy"
    sigfiles="signal_out.shot1.su"
elif ${is_fdt} ; then
    sufiles="seis_fdt4_curl.su seis_fdt4_div.su seis_fdt4_p.su seis_fdt4_vx.su seis_fdt4_vy.su"
    snapfiles="snap_fdt4.bin.curl snap_fdt4.bin.div snap_fdt4.bin.p snap_fdt4.bin.vx snap_fdt4.bin.vy"
    sigfiles="signal_out.shot1.su"
else
    sufiles="seis_curl.su seis_div.su seis_p.su seis_vx.su seis_vy.su"
    snapfiles="snap.bin.curl snap.bin.div snap.bin.p snap.bin.vx snap.bin.vy"
    sigfiles="signal_out.shot1.su"
fi

for file in ${sufiles} ; do
    tfile="${testdir}/${file}"
    rfile="${refdir}/${file}"
    dfile="${testdir}/${file}_diff"
    if ! test -r ${tfile} ; then
	echo "Error: ${tfile} not found. Skipping comparison."
	continue
    fi
    if ! test -r ${rfile} ; then
	echo "Error: ${rfile} not found. Skipping comparison."
	continue
    fi
    echo "Comparing file ${file}..."
    sudiff ${tfile} ${rfile} > ${dfile}
    cat ${tfile} ${rfile} ${dfile} | sushw key=tracl a=1 b=1 j=0 | \
        suxwigb perc=99.9 key=tracl title="${file} (test result, ref. result, diff. : 3x8 traces)"   
done

for file in ${sigfiles} ; do
    tfile="${testdir}/${file}"
    rfile="${refdir}/${file}"
    dfile="${testdir}/${file}_diff"
    if ! test -r ${tfile} ; then
	echo "Error: ${tfile} not found. Skipping comparison."
	continue
    fi
    if ! test -r ${rfile} ; then
	echo "Error: ${rfile} not found. Skipping comparison."
	continue
    fi
    echo "Comparing file ${file}..."
    sudiff ${tfile} ${rfile} > ${dfile}
    cat ${tfile} ${rfile} ${dfile} | sushw key=tracl a=1 b=1 j=0 | \
        suxwigb perc=99.9 key=tracl title="${file} (test result, ref. result, diff. : 3x1 trace)"   
done

fwi=$(echo ${weq} | grep -c "_fwi")
if [[ ${fwi} -eq 0 ]] ; then
  nsamp=200
else
  nsamp=400
fi

for file in ${snapfiles} ; do
    tfile="${testdir}/${file}"
    rfile="${refdir}/${file}"
    dfile="${testdir}/${file}_diff"
    if ! test -r ${tfile} ; then
	echo "Error: ${tfile} not found. Skipping comparison."
	continue
    fi
    if ! test -r ${rfile} ; then
	echo "Error: ${rfile} not found. Skipping comparison."
	continue
    fi
    echo "Comparing file ${file}..."
    tmp1=$(mktemp -p .)
    suaddhead ns=${nsamp} < ${tfile} > ${tmp1}
    tmp2=$(mktemp -p .)
    suaddhead ns=${nsamp} < ${rfile} > ${tmp2}
    tmp3=$(mktemp -p .)
    sudiff ${tmp1} ${tmp2} > ${tmp3}
    cat ${tmp1} ${tmp2} ${tmp3} | sushw key=tracl,dt,d1 a=1,2,2 b=1,0,0 | \
        suximage perc=99.9 wbox=1800 hbox=200 title="${file} (test result, ref. result, diff. : 3x(3x200) traces)"
    echo "Maximum & minimum amplitude of reference data:"
    sumax < ${tmp1} mode=maxmin
    echo "Maximum & minimum amplitude of difference:"
    sumax < ${tmp3} mode=maxmin
    echo "RMS amplitude of reference data:"
    sumax < ${tmp1} mode=rms
    echo "RMS amplitude of difference:"
    sumax < ${tmp3} mode=rms
    rm -f ${tmp1} ${tmp2} ${tmp3} 
done
