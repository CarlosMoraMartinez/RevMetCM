cd /scratch/groups/assembly/shared/projects/cmora/results/mash1/

cd ill_sketch/
mash paste all_ill *.msh

cd ../ont_sketch/
mash paste all_ont *.msh

cd ../
mkdir all2all
mv ont_sketch/all_ont.msh ill_sketch/all_ill.msh all2all
cd all2all

mash dist all_ont.msh all_illumina.msh > on2ill.dist
mash dist all_ont.msh all_ont.msh > on2on.dist
mash dist all_illumina.msh >  ill2ill.dist


