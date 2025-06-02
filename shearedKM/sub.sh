#! /bin/bash


rho0=1
T=0.5

sigma=0.
h=0.1
n_step=5000000
snap_dt=50000
seed=1000
IC=resume

for D_theta in 0
do
    job_name="L256_s${sigma}_D${D_theta}_h${h}"
    cat <<EOF > ${job_name}.sh
#!/bin/bash
#SBAtCH --job-name=${job_name}
#SBATCH --partition=localhost
#SBATCH -N 1 -n 1
#SBATCH --hint=nomultithread

srun ./a.out $rho0 $T $D_theta $sigma $h $n_step $snap_dt $seed $IC
EOF
sleep 0.25
sbatch ${job_name}.sh
done
