1. Legge scriptet på .FHI_scripts, deretter cd til ønsket mappe.
2. Ikke kopier pos og neg kontroll
3. Bare kopiere de over 97% coverage
4. Legge inn argument for hvilken mappe i newdir som skal brukes

594, 595, 597, 

# Basedir is the directory called e.g. Run598 which is inside Run598_Corona.
# Inside basedir are the directories for each sample, e.g. ArticXXX--
basedir=$(pwd)

# Newdir is the new location for the renamed files
newdir="/mnt/N/Virologi/NGS/1-NGS-Analyser/1-Rutine/2-Resultater/SARS-CoV-2/4-ELIXIRsubmisjon/raw/Run598"

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
        cp $R1 "${newdir}/${unique}_${year}_R1.fastq.gz"
        cp $R2 "${newdir}/${unique}_${year}_R2.fastq.gz"

        # Go back and repeat for next sample
        cd "${basedir}"
done
