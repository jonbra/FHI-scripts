# Future additions:
: <<'END'
- Ikke kopier pos og neg kontroll
- Bare kopiere de over 97% coverage
END

# NB! The script only works for Artic samples. Each sample should be in a separate directory named ArticXXX.
# To run the script you need to be inside a directory which contains sub-directories with fastq-files for each sample. 
# The fastq files must end with .fastq.gz
# E.g. For the Illumina Artic run 598 you should be inside the directory called Run598 which is inside Run598_Corona.

# Run the script with the command: bash ~/.fhiscripts/copy_and_rename_for_Elixir.sh -n "fastq_for_Elixir"

PPOSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
	-n|--newdir)
	newdir="$2"
	shift # past argument
	shift # past value
	;;	
esac
done

# Setting basedir (the directory which contains directories for each sample)
basedir=$(pwd)

# Newdir is the new location for the renamed files
fullnewdir="/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-ELIXIRsubmisjon/raw/${newdir}"

# The year for the run. Will be part of the filename
year="2021"

for dir in $(ls -d Artic*/)
do
        cd ${dir}
        R1=$(ls *_R1*.fastq.gz)
        R2=$(ls *_R2*.fastq.gz)
        
        # Extract unique number from file name
        unique="${R1:9:4}"

        # Copy and rename files
        cp $R1 "${fullnewdir}/${unique}_${year}_R1.fastq.gz"
        cp $R2 "${fullnewdir}/${unique}_${year}_R2.fastq.gz"

        # Go back and repeat for next sample
        cd "${basedir}"
done
