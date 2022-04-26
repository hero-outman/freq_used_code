# merge sample reps bam
samtools merge -@ 20 A8301.merged.bam A8301_a.bam A8301_b.bam
samtools merge -@ 20 B14_2.merged.bam B14_2_a.bam B14_2_b.bam
samtools merge -@ 20 control.merged.bam control_a.bam control_b.bam
samtools merge -@ 20 Eli_2.merged.bam Eli_2_a.bam Eli_2_b.bam
samtools merge -@ 20 TGF.merged.bam TGF_a.bam TGF_b.bam

suffix=".noDup.nochrM.sortpos.bam"
for f in  ./*$suffix
do
    root=`basename $f`

    if [[ $root == *R1.fq* ]]
    then # PE
	bf=`echo $root | sed -r 's/_R1.fq//g' | sed -r s/.gz//g` # control_a
        tt=`echo $bf_cut_R*.fq.gz` # control_a.bam
        p1=`echo $f` # control_a_cut_R1.fq.gz
	p2=`echo $f | sed 's/\_R1/_R2/g'` # change control_a_cut_R1.fq.gz to control_a_cut_R2.fq.gz
        if ! [ -f $tt ] # -f check a file if exist; if output.bam don.t exist, do align; if script won't work,possible target bam already exist(like 0 size)
        then
            echo PE.p33 ... $tt $bf
            s_dt=$(date)
            echo $s_dt $bf' start......'
            # $bf for sample name and out put name
            # $p1 for pair1 $p2 for pair2
            # in do_bt2.pe.hg38.sh script use $1 $2 $2 to catch arguments
            bash ./trimmomatic.sh $bf $p1 $p2
            f_dt=$(date)
            echo $f_dt $bf' finished......'
            sleep 1
        fi
    fi
done




suffix=".noDup.nochrM.sortpos.bam"
for f in  ./*$suffix
do
    root=`basename $f`

    echo $root
done  