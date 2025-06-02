#! /bin/bash


rho0=1.5

sigma=0.05
h=0.05
n_step=4000000
snap_dt=40000
D_t=0.01
seed=3001
IC=ordered

for T in 0.1
do
for gamma in 0.02
do
    job_name="L256_s${sigma}_T${T}_g${gamma}_h${h}"
    cat <<EOF > ${job_name}.sh
#!/bin/bash
#SBAtCH --job-name=${job_name}
#SBATCH --partition=localhost
#SBATCH -N 1 -n 8
#SBATCH --hint=nomultithread


mpirun --bind-to none ./a1_8.out $rho0 $T $D_t $gamma $sigma $h $n_step $snap_dt $seed $IC
EOF
sleep 0.25
sbatch ${job_name}.sh
done
done
