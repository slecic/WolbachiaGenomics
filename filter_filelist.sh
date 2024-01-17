for sample in $(cat /mnt/raid10/cerasi/slecic/cerasi/cline/fst/vcflist.txt); do 
bash /mnt/raid10/cerasi/slecic/cerasi/cline/fst/flysnp_filterring.sh ${sample}; 
done
