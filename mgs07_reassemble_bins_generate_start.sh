################################################
## 2020 11 08.
## For the meta-genome analysis
##
## mgs_04 binning by using Maxbin 
##  on super computer system
##
## https://sourceforge.net/projects/maxbin/files/
##
################################################

#!/bin/bash

## =================================== ##
## Bash que list
quick=( QUICK 4 72 )
small=( SMALL 12 30 )
apcq=( APC 30 40 )
apcmax=( APC 40 360 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}

contig="/aptmp/hiroakif/mgs02_genome_assembly"
abun="/aptmp/hiroakif/mgs01_fastp"

out="mgs10_binning_metawrap"
#out2="mgs05_check_bin_quality"


## --------------------------------------------------------  ##
mkdir "/aptmp/hiroakif/${out}" 
#mkdir "/aptmp/hiroakif/${out2}"
mkdir metagenome/${out}

if [ ! -d "/aptmp/hiroakif/${out}/fastas" ]; then 
  mkdir "/aptmp/hiroakif/${out}/fastas"
  cp /aptmp/hiroakif/mgs01_fastp/*fastq.gz /aptmp/hiroakif/mgs10_binning_metawrap/fastas/
  rename _R1.fastq.gz _1.fastq /aptmp/hiroakif/mgs10_binning_metawrap/fastas/*_R1.fastq.gz
  rename _R2.fastq.gz _2.fastq /aptmp/hiroakif/mgs10_binning_metawrap/fastas/*_R2.fastq.gz
fi

if [ ! -d "/aptmp/hiroakif/${out}/fasqsr" ]; then 
  mkdir "/aptmp/hiroakif/${out}/fasqsr"
  cp /aptmp/hiroakif/mgs01_fastp/*R2.fastq.gz /aptmp/hiroakif/mgs10_binning_metawrap/fasqsr/
  rename _R2.fastq.gz _2.fastq /aptmp/hiroakif/mgs10_binning_metawrap/fasqsr/*_R2.fastq.gz
fi

if [ ! -d "/aptmp/hiroakif/${out}/fasqsf" ]; then 
  mkdir "/aptmp/hiroakif/${out}/fasqsf"
  cp /aptmp/hiroakif/mgs01_fastp/*R1.fastq.gz /aptmp/hiroakif/mgs10_binning_metawrap/fasqsf/
  rename _R1.fastq.gz _1.fastq /aptmp/hiroakif/mgs10_binning_metawrap/fasqsf/*_R1.fastq.gz
fi

for i in `ls ${contig} | grep scaffold | cut -b 1-7 | uniq`
	do
	
	scaf=`ls ${contig} | grep ${i} | grep scaffold`

cat <<EOF > metagenome/${out}/${i}_${out}.qsub
#!/bin/bash
#PBS -q ${que}
#PBS -N ${i}_${out}
#PBS -l select=1:ncpus=${ncpu}:mem=${mem}gb
#PBS -e metagenome/${out}
#PBS -j eo
## -------------------------------------------------------- ##

source /etc/profile.d/modules.sh
module load metaWRAP

cd /aptmp/hiroakif/

assemblef=${contig}/${scaf}
## -------------------------------------------------------- ##

out=${out}
i=${i}
que=${small[0]}
ncpu=${small[1]}
mem=${small[2]}

## -------------------------------------------------------- ##

EOF

cat <<'EOF' >> metagenome/${out}/${i}_${out}.qsub
mkdir ${out}/${i}
mkdir ${out}/${i}

fasqs=` ls -d ${out}/fastas/* `
fasqsf=` ls -d ${out}/fasqsf/* `
fasqsr=` ls -d ${out}/fasqsr/* `

# metawrap binning -o ${out}/${i} -t ${ncpu} -m ${mem} -a ${assemblef} --metabat2 --maxbin2 --concoct $fasqs
# metawrap bin_refinement --quick -o ${out}/${i}/refinement -t ${ncpu} -m ${mem} -A ${out}/${i}/metabat2_bins/ -B ${out}/${i}/maxbin2_bins/ -C ${out}/${i}/concoct_bins/ -c 50 -x 5
# metawrap blobology -a ${assemblef} -t ${ncpu} -o ${out}/${i}/blobology --bins ${out}/${i}/refinement/metawrap_50_5_bins $fasqs
# metawrap quant_bins -b ${out}/${i}/refinement/metawrap_50_5_bins -o ${out}/${i}/quant_bins -a ${assemblef} $fasqs
metawrap reassemble_bins -o ${out}/${i}/reassemble -1 `ls $fasqs | grep _1.fastq | grep ${i}` -2 `ls $fasqs | grep _2.fastq | grep ${i}` -t ${ncpu} -m ${mem} -c 50 -x 5 -b ${out}/${i}/refinement/metawrap_50_5_bins


EOF

	chmod +x metagenome/${out}/${i}_${out}.qsub
	qsub metagenome/${out}/${i}_${out}.qsub

	echo " "
	echo "********* INFO: ${i}_${out}.qsub was generated and started! *********"

done
## -------------------------------------------------------- ##
