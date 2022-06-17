#----excecute with LAMMPI
#lamboot -v mpihosts
#lamboot -v

# mpirun -np 4 nice -19 ../bin/sofi2D ./in_and_out/sofi2D.json | tee ./in_and_out/sofi2D.jout
# mpirun -np 4 ../bin/sofi2D ./in_and_out/sofi2D.json | tee ./in_and_out/sofi2D.jout
mpirun -np 4 ../bin/sofi2D ./in_and_out/sofi2D.json > ./in_and_out/sofi2D.jout
../bin/snapmerge ./in_and_out/sofi2D.json

xmovie n1=200 n2=200  < ./snap/test.bin.vy loop=1 clip=3.0e-10 title=%g &

#----execute with OPENMPI2
#mpirun  --hostfile mpihosts -np 4 nice -19  ../bin/sofi2D ./in_and_out/sofi2D.json | tee ./in_and_out/sofi2D.jout


#merge snapshots
#../bin/snapmerge ./in_and_out/sofi2D.json

