### install tools as indicated in the README

### prepare databases


1. General Requirements

System dependencies

sudo apt-get update
sudo apt-get install -y wget curl tar gzip unzip python3


2.1 Install Kraken2

git clone https://github.com/DerrickWood/kraken2.git
cd kraken2
./install_kraken2.sh .
export PATH=$PWD:$PATH

2.2 Download Standard Kraken2 DB (recommended)
kraken2-build --standard --db kraken2_db

This includes:

bacteria
archaea
viral genomes
human genome (depending on build)
plasmids

2.3 Custom RefSeq / NR-style database

Download RefSeq genomes
kraken2-build --download-library bacteria --db kraken2_db
kraken2-build --download-library viral --db kraken2_db
kraken2-build --download-library archaea --db kraken2_db
Build database
kraken2-build --build --db kraken2_db

3. NCBI RefSeq / NR Database (for DIAMOND / HUMAnN / BLAST)

Option A: DIAMOND NR-style protein DB
wget ftp://ftp.ncbi.nlm.nih.gov/blast/db/FASTA/nr.gz
gunzip nr.gz

Build DIAMOND index:

diamond makedb --in nr --db nr
Option B: RefSeq protein database
wget ftp://ftp.ncbi.nlm.nih.gov/refseq/release/protein/protein.*.protein.faa.gz

Concatenate:

cat protein.*.faa.gz > refseq_protein.faa.gz
gunzip refseq_protein.faa.gz

4. SILVA rRNA Database (16S / 18S)

Download SILVA
wget https://www.arb-silva.de/fileadmin/silva_databases/release_138_2/Exports/SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
gunzip SILVA_138.2_SSURef_NR99_tax_silva.fasta.gz
Optional (QIIME2-ready format)
wget https://data.qiime2.org/2023.9/common/silva-138-99-seqs.qza
wget https://data.qiime2.org/2023.9/common/silva-138-99-tax.qza

5. MetaPhlAn Database (bioBakery)

Install MetaPhlAn
conda install -c bioconda metaphlan
Download database
metaphlan --install

or manually:

metaphlan --install --bowtie2db metaphlan_db
Set database path
metaphlan --bowtie2db /path/to/metaphlan_db

6. HUMAnN (bioBakery) Databases

Install HUMAnN
conda install -c bioconda humann
Download ChocoPhlAn (nucleotide DB)
humann_databases --download chocophlan full ./humann_db
Download UniRef protein database
humann_databases --download uniref90_diamond full ./humann_db
Configure HUMAnN
humann_config --update database_folders nucleotide ./humann_db/chocophlan
humann_config --update database_folders protein ./humann_db/uniref

7. PhaBOX (Phage analysis database)

Install PhaBOX
git clone https://github.com/KennthShang/PhaBOX.git
cd PhaBOX
Setup dependencies
conda env create -f environment.yml
conda activate phabox
Databases

PhaBOX uses:

viral reference databases
host prediction models
protein annotation databases

Download instructions depend on module:

python download_db.py

(or from official server if provided in updates)

8. Suggested Directory Structure

databases/
│
├── kraken2/
├── silva/
├── metaphlan/
├── humann/
├── refseq_nr/
├── phabox/

9. Notes & Best Practices

Always use absolute paths in config files
Kraken2 databases can exceed 50–200 GB
HUMAnN requires both:
ChocoPhlAn (DNA)
UniRef (protein)
SILVA is mainly for 16S/18S amplicon studies
MetaPhlAn is marker-gene based (faster, lighter than Kraken2)

10. Summary

Tool	Purpose	Database
Kraken2	Taxonomic classification	RefSeq / custom k-mer DB
SILVA	rRNA taxonomy	16S/18S reference
MetaPhlAn	microbial profiling	marker genes
HUMAnN	functional profiling	UniRef + ChocoPhlAn
PhaBOX	phage analysis	viral + host models
NCBI NR	protein annotation	DIAMOND / BLAST

#### Pipeline start

This workflow is not automated and requires manual execution of each module independently.

#!/bin/bash

#### meta pipeline ###
# params
rawfolder="/raw_reads"
workfile="/work"

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
        singularity exec images_singularity/denovo.sif \
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

    singularity exec images_singularity/meta.sif \
fastp \
 -i "${raw_dir}/${sample_id}.fastq.gz" \
 -o "${trimmed_dir}/${sample_id}.fastq.gz" \
 --length_required 30  --detect_adapter_for_pe\
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
singularity exec images_singularity/meta.sif bash -c '
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
export TMPDIR=/tmp
singularity exec images_singularity/srahumanscrubber.sif bash -c '
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
        singularity exec images_singularity/denovo.sif \
        reformat.sh threads=114 in="dehosted/${sampleId}.fastq.gz" out1="$output_R1" out2="$output_R2" int=t qin=auto qout=33
    else
        echo "Skipping Kraken2 reformatting for ${sampleId}, already done."
    fi
done < "sampleID.txt"

mkdir -p krak2
input_directory="${workfile}/paired"
kraken2_db="kraken_db"
while IFS= read -r sample_id; do
    output_report="krak2/${sample_id}.tsv"
    if [[ ! -f "$output_report" ]]; then
        singularity exec images_singularity/meta.sif \
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
        singularity exec images_singularity/krakentools.sif \
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
SINGULARITY_IMG="images_singularity/denovo.sif"

# **Interleave**
while IFS= read -r sampleId; do
    output_interleaved="${INTERLEAVED_DIR}/${sampleId}.fastq"
    if [[ ! -f "$output_interleaved" ]]; then
        echo "Interleaving: $sampleId"
        singularity exec "$SINGULARITY_IMG" \
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
        singularity exec "$SINGULARITY_IMG" \
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
        singularity exec images_singularity/denovo.sif \
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
singularity exec images_singularity/denovo.sif \
reformat.sh minlength=$ in=denovo/cdhit/allcontigs.fasta out=denovo/cdhit/allcontigs_filtered.fasta

singularity exec images_singularity/cd-hit.sif \
cd-hit-est -i denovo/cdhit/allcontigs_filtered.fasta -o denovo/cdhit/allcontigs_filtered_500_dereplicated.fasta -c 0.95 -G 0 -aS 0.95 -g 1 -r 1 -M 0 -d 0 -T 40


# Phabox2

conda activate phabox2
phabox2 --task end_to_end --dbdir /phabox_db_v2_1 --outpth  ${workfile}/phabox --contigs denovo/cdhit/allcontigs_filtered_500_dereplicated.fasta --threads 40 --len 500


# diamond for bacteria

db2="refseq_bacteria"
diamond blastx --query "allcontigs_filtered_1000_dereplicated.fasta" -p 118 -e 1e-3 --id 85 --top 1 --db $db2 --taxonlist 2 --outfmt 102 --out "diamond.txt"

tax2lin="ncbitax2lin_final.txt"
cat diamond.txt | awk -F'\t' '{print $2 "\t" $1}' > LCA_OTU_sorted.txt
sort LCA_OTU_sorted.txt >> LCA_OTU_sorted2.txt
join -a1 -o 2.2,1.2 -t$'\t' LCA_OTU_sorted2.txt $tax2lin | grep "Bacteria" > bacterial_contigs.tsv
cat allcontigs_filtered_1000_dereplicated.fasta |awk '/^>/ {if(N>0) printf("\n"); printf("%s ",$0);++N;next;} { printf("%s",$0);} END {printf("\n");}' > allcontigs_lin.fasta
cat viral_contigs.tsv | awk -F'\t' '{print $2}'|uniq > bacterial_contigs_list.txt
grep -f bacterial_contigs_list.txt allcontigs_lin.fasta | sed 's/ /\n/'g > all_bacterial_contigs.fasta


//sortmerna
for i in *.fastq* ; do  sortmerna --ref silva/SILVA_138.2_NR99_LSU_SSU_Ref.fasta  --reads $i --aligned /aligned/$i --fastx --other /nonaligned/$i -v --workdir /dir/ --threads 118 -m 80000 
rm -r /kvdb ; done

