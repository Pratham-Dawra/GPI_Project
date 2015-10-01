#!/bin/bash
i=0
while [ "$i" -lt "90" ]
do
i=`expr $i + 1`
I=`echo $i | awk '{printf("%i.0",$1*2+38)}'`

echo running Shot $I $i

cp source_org.dat source.dat
replace 42.0 $I  source.dat

../bin/fdelast_ssg_4th < fdveps_ku1500.inp > fdveps_ku1500.$i.out

cp su/ku1500_div.su.0 su/ku1500_div_$i.su.0
cp su/ku1500_rot.su.0 su/ku1500_rot_$i.su.0
cp su/ku1500_p.su.0 su/ku1500_p_$i.su.0
cp su/ku1500_x.su.0 su/ku1500_x_$i.su.0
cp su/ku1500_y.su.0 su/ku1500_y_$i.su.0



done
