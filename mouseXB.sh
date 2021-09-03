#!/bin/bash
# The interpreter used to execute the script

#â€œ#SBATCHâ€? directives that convey submission options:

#SBATCH --job-name=mouseXB1GAO
#SBATCH --mail-user=bmarzban@umich.edu
#SBATCH --mail-type=END
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=30
#SBATCH --mem-per-cpu=1200
#SBATCH --time=48:00:00
#SBATCH --account=beardda99
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j.log
#SBATCH --output=matlab_parfor.out
#SBATCH --error=matlab_parfor.err

cd /home/bmarzban/mouseXB3
module load matlab/R2018b
matlab -nodisplay -r -noFigureWindows -nosplash "driver_GAO" > mouseXBGAO.out

