

#!/bin/bash
## =================================== ##
## Bash que list
quick=( QUICK 4 72 )
small=( SMALL 12 108 )
apcq=( APC 20 50 )
apcmax=( APC 40 360 )

## =================================== ##
que=${apcq[0]}
ncpu=${apcq[1]}
mem=${apcq[2]}


# Run on local pc

cd /Users/hiroakifujita/Desktop/Microbiome_TimeSeries/metagenome/mgs10_binning_metawrap/results/prokka_rep.bins

for i in `ls | grep .fa$`
do
	cd ${i}

	id=`echo ${i} |  awk -F ".fa$" '{print $1}' `

 	sed -e "s/>/>${id}_/g" *.faa > ../${id}_prokka.faa
 		
 	cd ../
 	 	
done
 
cat *_prokka.faa > ../prokka_summary_rep.bin.faa

echo done

cd /Users/hiroakifujita/Desktop/Microbiome_TimeSeries/metagenome/mgs10_binning_metawrap/results/prokka

for i in `ls | grep .fa$`
do
	cd ${i}
	echo $i
	id=`echo ${i} |  awk -F ".fa$" '{print $1}' `

 	sed -e "s/>/>${id}_/g" *.faa > ../${id}_prokka.faa
 		
 	cd ../
 	 	
done
 
cat *_prokka.faa > ../prokka_summary.faa