#!/bin/bash

#### meta pipeline ###
# params
rawfolder="/srv/scratch/ananihu/PREVHAP_REVELO_FINAL/raw_reads"
workfile="/srv/scratch/ananihu/PREVHAP_REVELO_FINAL/work"

# 0) Create sampleID
cd ${rawfolder}
samples=$(ls *_R[12].fastq.gz 2>/dev/null | sed -E 's/_R[12]\.fastq\.gz//' | sort -u)
printf '%s\n' "$samples" > "${workfile}/sampleID.txt"

cd ../work/

# 1) Interleave
mkdir -p ${workfile}/raw_interleaved
while IFS= read -r sampleId; do
    output_file="${workfile}/raw_interleaved/${sampleId}.fastq.gz"
    if [[ ! -f "$output_file" ]]; then
        singularity exec --bind /srv/scratch /srv/scratch/ananihu/denovo.sif \
        reformat.sh threads=40 in1="${rawfolder}/${sampleId}_R1.fastq.gz" in2="${rawfolder}/${sampleId}_R2.fastq.gz" out="$output_file" qin=auto qout=33
    else
        echo "Skipping interleaving for ${sampleId}, already done."
    fi
done < "sampleID.txt"

# 2) Trimming
mkdir -p trimmed
mkdir -p reports

while IFS= read -r sample_id; do
    input_fastq="raw_interleaved/${sample_id}.fastq.gz"
    output_fastq="trimmed/${sample_id}.fastq.gz"
    report_html="reports/${sample_id}_fastp.html"
	report_json="reports/${sample_id}_fastp.json"

    if [[ -f "$output_fastq" ]]; then
        echo "Skipping trimming for ${sample_id}, already done."
        continue
    fi

    singularity exec --bind /srv/scratch:/srv/scratch /srv/scratch/ananihu/meta.sif \
fastp \
 -i "${raw_dir}/${sample_id}.fastq.gz" \
 -o "${trimmed_dir}/${sample_id}.fastq.gz" \
 --length_required 30 \
 -w 40 \
 --html "${reports_dir}/${sample_id}_fastp.html" \
 --json "${reports_dir}/${sample_id}_fastp.json"


    if [ $? -eq 0 ]; then
        echo "Trimming complete for sample: ${sample_id}"
    else
        echo "Error: Trimming failed for sample: ${sample_id}"
    fi
done < "sampleID.txt"

# 3) Summarize QC results
singularity exec --bind /srv/scratch /srv/scratch/ananihu/meta.sif bash -c '
mkdir -p reports
out_file="reports/fastp_summary.tsv"

# Header
echo -e "SampleID\tTotalReads\tFilteredReads\tReadsRetained(%)\tQ20(%)\tQ30(%)\tGC(%)\tDuplication(%)" > "$out_file"

# Loop over samples
while IFS= read -r sample_id; do
    json_file="reports/${sample_id}_fastp.json"
    if [[ -f "$json_file" ]]; then
        total_reads=$(jq -r ".summary.before_filtering.total_reads" "$json_file")
        filtered_reads=$(jq -r ".summary.after_filtering.total_reads" "$json_file")
        retained=$(awk -v t="$total_reads" -v f="$filtered_reads" 'BEGIN { if (t>0) printf "%.2f", (f/t)*100; else print "NA" }')
        q20=$(jq -r ".summary.after_filtering.q20_rate * 100" "$json_file" | awk '{printf "%.2f", $1}')
        q30=$(jq -r ".summary.after_filtering.q30_rate * 100" "$json_file" | awk '{printf "%.2f", $1}')
        gc=$(jq -r ".summary.after_filtering.gc_content * 100" "$json_file" | awk '{printf "%.2f", $1}')
        dup=$(jq -r ".duplication.rate * 100" "$json_file" | awk '{printf "%.2f", $1}')

        echo -e "${sample_id}\t${total_reads}\t${filtered_reads}\t${retained}\t${q20}\t${q30}\t${gc}\t${dup}" >> "$out_file"
    else
        echo "⚠️  Warning: JSON report missing for ${sample_id}" >&2
    fi
done < "sampleID.txt"

echo
echo "✅ QC summary generated: $out_file"
'

# 3) Dehost
mkdir -p dehosted
export TMPDIR=/srv/scratch/ananihu/tmp
singularity exec --bind /srv/scratch /srv/scratch/ananihu/srahumanscrubber.sif bash -c '
while IFS= read -r sampleId || [[ -n "$sampleId" ]]; do
    input_file="trimmed/${sampleId}.fastq.gz"
    output_file="dehosted/${sampleId}.fastq.gz"

    if [[ -f "$output_file" ]]; then
        echo "Skipping dehosting for ${sampleId}, already done."
        continue
    fi

    if [[ ! -f "$input_file" ]]; then
        echo "Error: Input file $input_file not found. Skipping ${sampleId}."
        continue
    fi

    if unpigz -p 16 -c "$input_file" | \
       /opt/scrubber/scripts/scrub.sh -p 16 -x -d "/opt/scrubber/data/human_filter.db" -o - -i - | \
       repair.sh threads=16 in=stdin out=stdout | \
       pigz -6 -p 16 -c > "$output_file"; then
        echo "Successfully processed ${sampleId}."
    else
        echo "Error processing ${sampleId}. Check logs for details."
    fi
done < sampleID.txt
'

# 4) Kraken2
mkdir -p paired
while IFS= read -r sampleId; do
    output_R1="paired/${sampleId}_R1.fastq.gz"
    output_R2="paired/${sampleId}_R2.fastq.gz"
    if [[ ! -f "$output_R1" || ! -f "$output_R2" ]]; then
        singularity exec --bind /srv/scratch /srv/scratch/ananihu/denovo.sif \
        reformat.sh threads=114 in="dehosted/${sampleId}.fastq.gz" out1="$output_R1" out2="$output_R2" int=t qin=auto qout=33
    else
        echo "Skipping Kraken2 reformatting for ${sampleId}, already done."
    fi
done < "sampleID.txt"

mkdir -p krak2
input_directory="${workfile}/paired"
kraken2_db="/srv/scratch/kraken2_db_files_iai/kraken_db_2"
while IFS= read -r sample_id; do
    output_report="krak2/${sample_id}.tsv"
    if [[ ! -f "$output_report" ]]; then
        singularity exec --bind /srv/ /srv/scratch/ananihu/meta.sif \
        kraken2 --db "$kraken2_db" \
                --memory-mapping \
                --confidence 0.1 \
                --minimum-hit-groups 2 \
                --classified-out "krak2/${sample_id}#.fastq" \
                --output "krak2/${sample_id}.out" \
                --report-minimizer-data \
                --report "$output_report" \
                --paired "$input_directory/${sample_id}_R1.fastq.gz" "$input_directory/${sample_id}_R2.fastq.gz" \
                --threads 114
    else
        echo "Skipping Kraken2 for ${sample_id}, already done."
    fi
done < "sampleID.txt"

# 5) Extract Eukaryotic Reads
mkdir -p krak2/krak2_filter
while IFS= read -r sample_id; do
    output_fastq1="krak2/krak2_filter/${sample_id}_R1.fastq"
    output_fastq2="krak2/krak2_filter/${sample_id}_R2.fastq"
    if [[ ! -f "$output_fastq1" || ! -f "$output_fastq2" ]]; then
        mv "krak2/${sample_id}_1.fastq" "krak2/${sample_id}_R1.fastq"
        mv "krak2/${sample_id}_2.fastq" "krak2/${sample_id}_R2.fastq"
        singularity exec --bind /srv/scratch /srv/scratch/ananihu/krakentools.sif \
        extract_kraken_reads.py --max 1000000000 -t 2759 -k "krak2/${sample_id}.out" -s1 "krak2/${sample_id}_R1.fastq" -s2 "krak2/${sample_id}_R2.fastq" -o "$output_fastq1" -o2 "$output_fastq2" -r "krak2/${sample_id}.tsv" --include-children --exclude --fastq-output
        echo "Extraction complete for sample: ${sample_id}"
    else
        echo "Skipping extraction for ${sample_id}, already done."
    fi
done < "sampleID.txt"

# 6) New Step: Interleave and Deinterleave Before Spades Assembly

mkdir -p krak2/krak2_filter/inter

# Define paths
INPUT_DIR="krak2/krak2_filter"
INTERLEAVED_DIR="${INPUT_DIR}/inter"
PAIR_DIR="${INTERLEAVED_DIR}/pair"

# Ensure directories exist
mkdir -p "$INTERLEAVED_DIR" "$PAIR_DIR"

# Singularity container path
SINGULARITY_IMG="/srv/scratch/ananihu/denovo.sif"

# **Interleave**
while IFS= read -r sampleId; do
    output_interleaved="${INTERLEAVED_DIR}/${sampleId}.fastq"
    if [[ ! -f "$output_interleaved" ]]; then
        echo "Interleaving: $sampleId"
        singularity exec --bind /srv/scratch "$SINGULARITY_IMG" \
            reformat.sh threads=70 \
            in1="${INPUT_DIR}/${sampleId}_R1.fastq" \
            in2="${INPUT_DIR}/${sampleId}_R2.fastq" \
            out="$output_interleaved" \
            qin=auto qout=33
    else
        echo "Skipping interleaving for ${sampleId}, already done."
    fi
done < "sampleID.txt"

# **Deinterleave**
while IFS= read -r sampleId; do
    output_R1="${PAIR_DIR}/${sampleId}_R1.fastq"
    output_R2="${PAIR_DIR}/${sampleId}_R2.fastq"
    if [[ ! -f "$output_R1" || ! -f "$output_R2" ]]; then
        echo "Deinterleaving: $sampleId"
        singularity exec --bind /srv/scratch "$SINGULARITY_IMG" \
            reformat.sh threads=70 \
            in="${INTERLEAVED_DIR}/${sampleId}.fastq" \
            out1="$output_R1" \
            out2="$output_R2"
    else
        echo "Skipping deinterleaving for ${sampleId}, already done."
    fi
done < "sampleID.txt"

echo "Processing complete!"

# 7) De novo Assembly (Spades)
mkdir -p denovo/spades
while IFS= read -r sampleId; do
    output_dir="denovo/spades/${sampleId}"
    if [[ ! -d "$output_dir" ]]; then
        singularity exec --bind /srv/scratch /srv/scratch/ananihu/denovo.sif \
        spades.py -t 114 -1 "${PAIR_DIR}/${sampleId}_R1.fastq" -2 "${PAIR_DIR}/${sampleId}_R2.fastq" -o "$output_dir" --meta
    else
        echo "Skipping de novo assembly for ${sampleId}, already done."
    fi
done < "sampleID.txt"

mkdir -p denovo/contigs
while IFS= read -r sampleId; do
    input_contig="denovo/spades/${sampleId}/contigs.fasta"
    output_contig="denovo/contigs/${sampleId}_contigs.fasta"
    if [[ ! -f "$output_contig" && -f "$input_contig" ]]; then
        cp "$input_contig" "$output_contig"
    else
        echo "Skipping copy of contigs for ${sampleId}, already done."
    fi
done < "sampleID.txt"

# 8) Dereplicate Contigs
mkdir -p denovo/cdhit
cat denovo/contigs/*.fasta > denovo/cdhit/allcontigs.fasta
singularity exec --bind /srv/scratch /srv/scratch/ananihu/denovo.sif \
reformat.sh minlength=$ in=denovo/cdhit/allcontigs.fasta out=denovo/cdhit/allcontigs_filtered.fasta

singularity exec --bind /srv/scratch /srv/scratch/ananihu/cd-hit.sif \
cd-hit-est -i denovo/cdhit/allcontigs_filtered.fasta -o denovo/cdhit/allcontigs_filtered_500_dereplicated.fasta -c 0.95 -G 0 -aS 0.95 -g 1 -r 1 -M 0 -d 0 -T 40


# Phabox2

conda activate phabox2
phabox2 --task end_to_end --dbdir /srv/scratch/ananihu/phabox_db_v2_1 --outpth  ${workfile}/phabox --contigs denovo/cdhit/allcontigs_filtered_500_dereplicated.fasta --threads 40 --len 500


# diamond for bacteria

db2="refseq_bacteria"
/srv/scratch/ananihu/diamond blastx --query "allcontigs_filtered_1000_dereplicated.fasta" -p 118 -e 1e-3 --id 85 --top 1 --db $db2 --taxonlist 2 --outfmt 102 --out "diamond.txt"

tax2lin="ncbitax2lin_final.txt"
cat diamond.txt | awk -F'\t' '{print $2 "\t" $1}' > LCA_OTU_sorted.txt
sort LCA_OTU_sorted.txt >> LCA_OTU_sorted2.txt
join -a1 -o 2.2,1.2 -t$'\t' LCA_OTU_sorted2.txt $tax2lin | grep "Bacteria" > bacterial_contigs.tsv
cat allcontigs_filtered_1000_dereplicated.fasta |awk '/^>/ {if(N>0) printf("\n"); printf("%s ",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' > allcontigs_lin.fasta
cat viral_contigs.tsv | awk -F'\t' '{print $2}'|uniq > bacterial_contigs_list.txt
grep -f bacterial_contigs_list.txt allcontigs_lin.fasta | sed 's/ /\n/'g > all_bacterial_contigs.fasta


//sortmerna
for i in *.fastq* ; do  sortmerna --ref /srv/scratch/ananihu/silva/SILVA_138.2_NR99_LSU_SSU_Ref.fasta  --reads $i --aligned /srv/scratch/ananihu/Big_IBIS/aligned/$i --fastx --other /srv/scratch/ananihu/Big_IBIS/nonaligned/$i -v --workdir /srv/scratch/ananihu/Big_IBIS/ --threads 118 -m 80000 
rm -r /srv/scratch/ananihu/Big_IBIS/kvdb ; done

