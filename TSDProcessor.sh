## Changelog ##
#25 March 2021 Jon Bråte
# Use hard links instead of soft links
# Save error report for the md5 checks
# Run the md5 checks at the end

#26 March 2021. Nacho Garcia
# chown included after moving the files

## Usage ##
#./TSDProcessor.sh -r FHI
# Substitute FHI with MIK or another lab that generated the sample. 
# The name must match exactly what is used in the different filenames.

PPOSITIONAL=()
while [[ $# -gt 0 ]]
do
key="$1"

case $key in
    -r|--run)
    run_id="$2"
    shift # past argument
    shift # past value
    ;;
esac
done


if [ -z "${run_id}" ]
then
	echo Error! No ID provided. 
	echo Try again with arguments like -r FHI or -r MIK
	exit 1
	
fi

## Move files ##

cd /tsd/p1516/data/durable/file-import/p1516-member-group/

#Raw files
for files in *Project*${run_id}*.tar
do
	mv ${files} /tsd/p1516/data/durable/${run_id}_NSC/fastq/
	mv ${files}.md5 /tsd/p1516/data/durable/${run_id}_NSC/fastq/
	chown :p1516-member-group ${files}
	chown :p1516-member-group ${files}.md5
	ln /tsd/p1516/data/durable/${run_id}_NSC/fastq/${files} /tsd/p1516/data/durable/file-export/${files} 
done

#Consensus
for files in *${run_id}*_variants.tar
do
	mv ${files} /tsd/p1516/data/durable/${run_id}_NSC/variants/
	mv ${files}.md5 /tsd/p1516/data/durable/${run_id}_NSC/variants/
	chown :p1516-member-group ${files}
	chown :p1516-member-group ${files}.md5
	ln /tsd/p1516/data/durable/${run_id}_NSC/variants/${files} /tsd/p1516/data/durable/file-export/${files} 
done

#Analysis
for files in *${run_id}*.tar
do
	mv ${files} /tsd/p1516/data/durable/${run_id}_NSC/analysis/
	mv ${files}.md5 /tsd/p1516/data/durable/${run_id}_NSC/analysis/
	chown :p1516-member-group ${files}
	chown :p1516-member-group ${files}.md5
	ln /tsd/p1516/data/durable/${run_id}_NSC/analysis/${files} /tsd/p1516/data/durable/file-export/${files} 
done

## Checksums ##

#Analysis
cd /tsd/p1516/data/durable/${run_id}_NSC/analysis/
md5sum -c *.md5 >> /tsd/p1516/data/durable/${run_id}_NSC/md5/md5sum_check.log
mv *.md5 /tsd/p1516/data/durable/${run_id}_NSC/md5/

#Raw files
cd /tsd/p1516/data/durable/${run_id}_NSC/fastq/
md5sum -c *.md5 >> /tsd/p1516/data/durable/${run_id}_NSC/md5/md5sum_check.log
mv *.md5 /tsd/p1516/data/durable/${run_id}_NSC/md5/

#Consensus
cd /tsd/p1516/data/durable/${run_id}_NSC/variants/
md5sum -c *.md5 >> /tsd/p1516/data/durable/${run_id}_NSC/md5/md5sum_check.log
mv *.md5 /tsd/p1516/data/durable/${run_id}_NSC/md5/
