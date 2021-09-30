# Lenka's notes:

- 12 Nov: setting up for the "synoptic diffusion project" using the full dynamics model, plan:
    - held - suarez setup
    - held - suarez setup + damping of all waves (following walker and schneider 05)
    - held - suarez setup + damping of some wavenumbers

- 30 Nov: segmentation error when using vert diff up/down so switching to the full idealized Caltech FMS
    - 

- 14 Dec:  
    - add explicit truncation using the input namelist
    - the wavenumber selection by divergence damping etc will have to wait (hs_syn+_diff_tr_old) 

- 30 Jun:
    - running slab ocean wavenumber experiments
    - default environment messed up - use srun to trigger from compute node for now
    - num_fourier = (no four modes - 1); make sure \geq ntasks

- 5 Jul:
    - several analysis jobs failed - faulty node: hpc-23-12 (slurmd: execve error -> specify a nodelist,e.g. 
in slurm script, need to add `#SBATCH nodelist=hpc-22-15`)
    - the above fails are also suffering from analysis having the wrong permissions (chmod)




