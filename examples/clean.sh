#!/bin/bash

for i in 1 2 3 4 5 6 7 8 ; do
    dir="weq${i}"
    if test -d ${dir} ; then
	echo "Removing all files in directory ${dir}..."
        rm -f ./${dir}/*
    else
        echo "Directory ${dir} not found. Continuing."
    fi
done
