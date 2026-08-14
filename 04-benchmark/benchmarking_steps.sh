## CORRE Guillaume @ GENETHON
## 2026-08-14


############################
# Benchmarking of GUIDE-seq v1.0.2,v2, iGUIDE-seq v1.2.0 and GENETHOFF
############################

## All pipelines have been installed according to their corresponding documentation on github.

## To reproduce this code, please install sra-tool, seqkit and cutadapt (in a conda environment for example)

############################
# Libraries used in the benchmark
############################

# we used the libraries from the iGUIDE-seq paper with gRNA VEGFAs2, VEGFAs3, TRAC, B2M and mock

# library metadata are available at : https://seqout.org/p/PRJNA506241 

    # B2M libraries IDs
    B2M="SRX5035779 SRX5035780 SRX5035781 SRX5035782 SRX5035783 SRX5035784 SRX5035785 SRX5035786 SRX5035787 SRX5035788 SRX5035912 SRX5035913 SRX5035914 SRX5035915 SRX5035916 SRX5035917 SRX5036060 SRX5036061 SRX5036062 SRX5036063 SRX5036066 SRX5036067 SRX5036069 SRX5036070 SRX5036092 SRX5036093 SRX5036094 SRX5036095 SRX5036096 SRX5036097 SRX5036100 SRX5036101"
    # TRAC libraries IDs
    TRAC="SRX5035821 SRX5035823 SRX5035824 SRX5035825 SRX5035826 SRX5035827 SRX5035828 SRX5035833 SRX5035834 SRX5035884 SRX5035885 SRX5035886 SRX5035887 SRX5035930 SRX5035931 SRX5035932 SRX5035933 SRX5035934 SRX5035935 SRX5035936 SRX5035937 SRX5035939 SRX5035940" 
    # VEGFAs2 libraries IDs
    VEGFAs2="SRX5035748  SRX5035749 SRX5035750  SRX5035751 SRX5035752 SRX5035754 SRX5035755 SRX5035756 SRX5035757 SRX5035758 SRX5035857 SRX5035858"
    # VEGFAs3 libraries IDs
    VEGFAs3="SRX5035944  SRX5035945 SRX5035946  SRX5035947 SRX5035948 SRX5035949 SRX5035952 SRX5035953 SRX5036086 SRX5036087 SRX5036088 SRX5036089"
    # Mock samples IDs
    mock="SRX5035805 SRX5035876 SRX5035877  SRX5035964 SRX5035965 SRX5035806 SRX5035809 SRX5035810 SRX5035890 SRX5035891 SRX5035892 SRX5035893 SRX5035894 SRX5035895 SRX5035896 SRX5035897 SRX5035997 SRX5035998 SRX5035999 SRX5036000 SRX5036001 SRX5036002 SRX5036003 SRX5036004 SRX5035900"
    
# raw sequences are downloaded from SRA using sra-tool prefetch and fasterq-dump commands 
    
# SRA data are already demultiplexed and have only the R1,R2 and I2 reads files. 
# For the pipelines to work end to end, we need an I1 file. We will make a fake I1 read for each gRNA using the folowing barcodes:

    declare -A barcode

    barcode[B2M]="AAAAAAAA"
    barcode[TRAC]="TTTTTTTT"
    barcode[VEGFAs2]="CCCCCCCC"
    barcode[VEGFAs3]="GGGGGGGG"
    barcode[mock]="CGCGCGCG"
        
        
    for grna in B2M TRAC VEGFAs2 VEGFAs3 mock; do
        libraries="${!grna}"
        I1=${barcode[$grna]}
        
        ### -->>> activate the environment containing sra-tool if necessary
        conda activate sra-tool 
        
        prefetch $grna
        
        # here we add the comment with read name as some pipeline do require a 2 parts read ID to work
        fasterq-dump --include-technical --seq-defline '@$sn 1:N:0:0' --qual-defline '+' SR*/
        

        ### -->>> activate the environment containing seqkit if necessary
        conda activate GENETHOFF
        
        ## Extract list of unique I2 barcodes (leading 8nt) per library and per gRNA
        
        find . -name "*_1.fastq" | while read f; do
            seqkit head -n 100 "$f" \
            | seqkit fx2tab \
            | awk -F'\t' -v file="$f" '{print file"\t"substr($2,1,8)}'
        done \
        | sort \
        | uniq -c \
        | awk '{print $2"\t"$3"\t"$1}' > $grna.I2.barcodes
        
        # merge all libraries for this gRNA
        cat *_1.fastq > $grna.I2.fastq 
        cat *_2.fastq > $grna.R1.fastq 
        cat *_3.fastq > $grna.R2.fastq 

        # make the fake I1 reads (replace R1 sequence with barcode)
        ./add_pattern_to_sequence.sh $I1 replace $grna.R1.fastq > $grna.I1.fastq

        # compress fastq files
        pigz $grna*.fastq

        # clean sra files
        rm -r SRR*

        # add barcodes to the barcode list for manifest or sample data sheet later
        awk -v bc="$I1" -v target="$grna" 'BEGIN{OFS="\t"; FS="\t"} {print $1,target,bc,$2,$3,$4}' $grna.I2.barcodes >> barcodes.txt
    done


## get unique barcode combination per target to prepare the different manifests
cut -f 2,3,4 barcodes.txt | uniq > unique.barcodes.txt


## build the "undetermined" files

cat *R1.fastq.gz* > Undetermined_R1.fastq.gz
cat *R2.fastq.gz* > Undetermined_R2_temp.fastq.gz
cat *I1.fastq.gz* > Undetermined_I1.fastq.gz
cat *I2.fastq.gz* > Undetermined_I2_temp.fastq.gz

### -->>> activate the environment containing cutadapt if necessary
# this step renames the R2 and I2 reads name to repalce 1:N:0:0 with 2:N:0:0
cutadapt -j 12 --rename '{id} 2:N:0:0' -o Undetermined_I2.fastq.gz Undetermined_I2_temp.fastq.gz
cutadapt -j 12 --rename '{id} 2:N:0:0' -o Undetermined_R2.fastq.gz Undetermined_R2_temp.fastq.gz

rm *temp.fastq.gz

### -->>> return to base environment
conda activate base







############################
# Run iGUIDE-seq v1.2.0
############################
cd /media/Data/common/iGUIDE_v1.2.0

# add paperReview.config.yml manifest file in configs/
# add paperReview.supp.csv and paperReview.sampleInfo.csv to sampleInfo/

conda activate iGuideSeq_v1.2.0

# configure the analysis
iguide setup configs/paperReview.config.yml
time iguide run configs/paperReview.config.yml --cores=24 -k -n








############################
# Run GUIDE-seq Tsai v1.0.2
############################

cd /media/Data/common/guideseq_tsai_v1.0.2/

conda activate guideseq_v1.0.2

# add the manifest file in paper_revision/
# update path to undetermined files

time python guideseq/guideseq.py all -m paper_revision/manifest_v1.0.2.yaml









############################
# Run GUIDE-seq Tsai v2
############################
# Tsai V2, we need to provide a control for each sample, which we don't have here.

# we run the demultiplexing step alone first and duplicate the control sample as background for all other samples

cd /media/Data/common/guideseq_tsai_v2/guideseq/paper_revision/

# add the manifest file to current folder
# add undetermined files to current folder

time docker run --rm  \
    -v /media/Data/common/guideseq_tsai_v2/guideseq/:/app \
    -v /media/Data/common/guideseq_tsai_v2/guideseq/paper_revision:/app/data \
    -v /media/References/Human/ensembl/GRCh38/Sequences/:/app/genome \
    -w /app/data liyc1989/tsailabsj \
    sh -c "python /app/guideseq/guideseq.py main -step demultiplex -m manifest_v2.yaml"

# it takes about 2 minutes
# real    2m27,798s
# user    0m0,033s
# sys     0m0,022s




# we just duplicated that control for each sample we have : 

sudo chown "$USER":"$USER" demultiplexed/

cd demultiplexed

ctrl=$(ls Control_*r1.fastq | sed 's/.r1.fastq//')

for sample in $(ls *.r1.fastq | sed 's/.r1.fastq//'); do

    [[ "$sample" =~ $ctrl ]] && continue

    
    if [[ ! -f "Control_${sample}.r1.fastq" ]]; then
        echo $sample    
        cp $ctrl.r1.fastq "Control_${sample}.r1.fastq"
        cp $ctrl.r2.fastq "Control_${sample}.r2.fastq"
    fi
done
 
 cd ..
 
 # we will start directly from the demultiplexed sample, so skiping the demultiplexing step
 
time docker run --rm  \
    -v /media/Data/common/guideseq_tsai_v2/guideseq/:/app \
    -v /media/Data/common/guideseq_tsai_v2/guideseq/paper_revision:/app/data \
    -v /media/References/Human/ensembl/GRCh38/Sequences/:/app/genome \
    -w /app/data liyc1989/tsailabsj \
    sh -c "python /app/guideseq/guideseq.py main --step  align+identify+visualize -m manifest_v2.yaml"
    
    


############################
# GENETHOFF
############################
conda activate base
cd /media/Data/gcorre/GENETHOFF/paper_review/

# add "undetermined" fastq.gz files (R1,R2,I1,I2)
# add sampleInfo.csv and configuration.yml to same folder
# update path to undetermined files in configuration file

time snakemake -s ../00-pipeline/genethOFF.snakemake -k -j 24 --use-conda
