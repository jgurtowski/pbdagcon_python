#!/bin/bash

# Simple pbdagcon workflow script.  Written for the benefit of running via
# smrtpipe so I can communicate pipe errors to the task.  We're overcoming
# the limitation of smrtpipe forcing tasks to run serially, enabling a new
# level of pipelining that's extremely efficient in an imperfect world ...

# generate pre-alignments
m4topre.py $mym4 $allm4 $subreads ${bestn-24} | \
# pipe it to consensus and generate fasta
pbdagcon -a -j ${nproc-15} | tee ${fasta-"corrected.fa"} | \
# generate a fastq
awk '{if($0~/>/){sub(/>/,"@",$0);print;}else{l=length($0);q="";while(l--){q=q "9"}printf("%s\n+\n%s\n",$0,q)}}' > ${fastq-"corrected.fq"}

# check the status of each piped command and exit non-zero if found
for exitval in ${PIPESTATUS[*]}
do
    if [ $exitval -gt 0 ]
    then
        exit $exitval
    fi
done

exit 0;
