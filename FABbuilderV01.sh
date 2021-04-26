#!/bin/bash


###########################################################
#
#    this is the scrpt for automated FAB homology-modeling modified by YT wang
#
############################################################
#
rm -rf L
rm -rf H
rm -rf temp
rm -rf result


# author : yanhong hong  modify by YT wang

export fasta_file="$1"

# use the build_profile.py script to search templates for the target sequence and map the target sequence to the templates
# output files are .prf file and .ali file

process_num=10   # the number of processes that this script will create.
evalue_cutoff=0.01
non_template_text='non_template.txt' # the path to write the non-template gene_id



# # read the output of the rename.sh script and write it into an array
mapfile -t jobids < <( bash rename.sh ${fasta_file} )


# # # # # # multiple processes
for jobid in "${!jobids[@]}"
do
	ali_jobids[$jobid]=${jobids[$jobid]}"/"${jobids[$jobid]}".ali"
	echo ${ali_jobids[$jobid]}
done | xargs -n 1 -I {} -P ${process_num} bash -c "python build_profile.py {}"


# # # # # #wait

# # # # # # # write the alignment results in a .log file
for jobid in "${!jobids[@]}"
do
	log_jobids[$jobid]=${jobids[$jobid]}"/"${jobids[$jobid]}".log"
    echo -e "index\tpdbid\tidentity" > ${log_jobids[$jobid]}
	prf_jobids[$jobid]=${jobids[$jobid]}"/"${jobids[$jobid]}".prf"
	#!$3==0 delete the target itself 
        #modified by YT wang    
	#cat ${prf_jobids[$jobid]}|grep -v "^#"|awk '{print $1 "\t" $2 "\t" $11}'|awk '!$3==0'|sort  -k3nr  >> ${log_jobids[$jobid]} #OLD
        ###YT 
        A[$jobid]=${log_jobids[$jobid]}".i"
        B[$jobid]=${log_jobids[$jobid]}".ii"
        C[$jobid]=${log_jobids[$jobid]}".iii"
        D[$jobid]=${log_jobids[$jobid]}".iiii"
        E[$jobid]=${log_jobids[$jobid]}".v" 
        cat ${prf_jobids[$jobid]}|grep -v "^#"|awk '{print $1 "\t" $5 "\t" $11}'|awk '$3==0' > ${A[$jobid]} 
        cat ${A[$jobid]} | awk '{print $2}'  >${B[$jobid]}
        leng=$(<${B[$jobid]})
        cat ${prf_jobids[$jobid]}|grep -v "^#"|awk '{print $1 "\t" $2 "\t" $5 "\t" $11}'|awk '!$3==0' > ${C[$jobid]}
        cat ${C[$jobid]} |while read line; do echo $line $leng; done > ${D[$jobid]} 
        cat ${D[$jobid]}|awk '{ if ($5 >= $3 )print $1 "\t" $2 "\t" $4*$3/$5; else print $1 "\t" $2 "\t" $4*$5/$3}' > ${E[$jobid]}
        cat ${E[$jobid]}|sort  -k3nr  >> ${log_jobids[$jobid]}
        ###YT     
done

# # # # # # # # extract the highest identity pdbid and write it to the template.txt
echo "the non-templates gene id will write to "${non_template_text}
for jobid in "${!jobids[@]}"
do
	log_jobids[$jobid]=${jobids[$jobid]}"/"${jobids[$jobid]}".log"
	if [ x`cat ${log_jobids[$jobid]}|wc|awk '{print $1}'` = x1 ];
	then
		echo "sorry no templates found for "${log_jobids[$jobid]}" with evalue cutoff "${evalue_cutoff}".please try HHM search method or use ab initio modeling methods."
		#delete the job from the quene
		#https://stackoverflow.com/questions/16860877/remove-an-element-from-a-bash-array
		delete=${jobids[$jobid]}
		echo ${delete} >> ${non_template_text}
		for target in "${delete[@]}"; do
  			for i in "${!jobids[@]}"; do
    			if [[ ${jobids[i]} = $target ]]; then
      				unset 'jobids[i]'
    			fi
  			done
	done
	else
		cat ${log_jobids[$jobid]}|awk 'NR==2 {print $2}' > ${jobids[$jobid]}/${jobids[$jobid]}_template.txt
	fi
done

#no job needs to be processed.exit.

if [ ${#jobids[@]} -eq 0 ]; then
    exit 0
fi


# # # # # # using a python script to fetch the template  structure with specified chain from pdb database.

for jobid in "${!jobids[@]}"
do

	template_text=${jobids[$jobid]}"/"${jobids[$jobid]}"_template.txt"
	outdir=${jobids[$jobid]}"/"${jobids[$jobid]}"_template"
	#echo "-o "${outdir}" -p "${template_text}
##modified by YT
        mkdir ${outdir}
        python3 download_pdbchain.py -o "${outdir}" -p "${template_text}"
done
#done | xargs -n 1 -I {} -P ${process_num} bash -c "python download_pdbchain.py {}"


# # # # # # Align the target with the template.

for jobid in "${!jobids[@]}"
do
	templatepath=${jobids[$jobid]}"/"${jobids[$jobid]}"_template"
	targetname=${jobids[$jobid]}
	templatename=`ls ${templatepath}`
	echo "--templatepath "${templatepath}"/"${templatename}" --targetname "${targetname}
done | xargs -n 1 -I {} -P ${process_num} bash -c "python align2d.py {}"





# # #start build model.

for jobid in "${!jobids[@]}"
do
	cd ${jobids[$jobid]}
	templatepath=${jobids[$jobid]}"_template"
	templatename=`ls ${templatepath}|sed 's/.pdb//g'`
	cp ${templatepath}/${templatename}.pdb .
	alignmentname=${jobids[$jobid]}"_"${templatename}".ali"
	echo -e ${jobids[$jobid]}"\n--alignmentname "${alignmentname}
	cd ..
done  |xargs -n 2  -P ${process_num}  -d'\n'  bash -c 'cd $0; python ../model_single.py $1'
#############YT wang #############FAB assemble
mkdir temp
cd temp
rm -rf *.pdb                                                                                                                                                  
cp ../L/L_template.txt ref                                                                                                                                    
id=$(<ref)                                                                                                                                                    
pdbid=${id:0:4}                                                                                                                                               
echo $pdbid                                                                                                                                                   
wget http://www.rcsb.org/pdb/files/$pdbid.pdb.gz                                                                                                             
gzip -d $pdbid.pdb.gz                                                                                                                                        
wget http://opig.stats.ox.ac.uk/webapps/newsabdab/sabdab/pdb/$pdbid/                                                                                          
mv index.html ref.pdb                                                                                                                                         
cp ../L/*1.pdb L.pdb                                                                                                                                          
cp ../H/*1.pdb H.pdb                                                                                                                                          
cp ../TMalign .                                                                                                                                               
./TMalign L.pdb ref.pdb -o LL                                                                                                                                 
./TMalign H.pdb ref.pdb -o HH                                                                                                                                 
cat HH_all_atm.pdb |awk '{print $0}'|awk '$5=="A"'   >TH                                                                                                      
awk '{if ($5 == "A") gsub(/A/,"H",$5); printf "%4s%7.0f%5s%4s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'<TH >THH    
cat LL_all_atm.pdb |awk '{print $0}'|awk '$5=="A"'   >TL                                                                                                      
awk '{if ($5 == "A") gsub(/A/,"L",$5); printf "%4s%7.0f%5s%4s%2s%4.0f%12.3f%8.3f%8.3f%6.2f%6.2f\n", $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $11}'<TL >TLL    
rm -rf candidate.pdb                                                                                                                                          
cat THH >>candidate.pdb                                                                                                                                       
cat ../TER >>candidate.pdb                                                                                                                                    
cat TLL >>candidate.pdb                                                                                                                                       
cat ../TER >>candidate.pdb
mkdir -p ../result
cp candidate.pdb ../result/
cd ..
#rm -rf L
#rm -rf H
#rm -rf i*
#rm -rf v
echo "job is down by FABbuilder01"                                                                                                                                    
