# Basedir is the directory called e.g. Run598_Corona.
# Inside basedir is a single directory called e.g. Run598
basedir=$(pwd)

# Newdir is the new location for the renamed files
newdir="/home/jonbra/FHI/Prosjekter/FHI-scripts/renaming_test_files/For_Elixir"

# Runname will be the same as the directory in basedir. E.g. Run598
runname=${basedir##*/}

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


