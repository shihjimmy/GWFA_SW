# Run under ./GWFA folder
CASES=(
  "out_sequence_err5%_10k:10240"
  "out_sequence_err5%_20k:20480"
  "out_sequence_err5%_40k:40480"
  "out_sequence_err5%_60k:60480"
  "out_sequence_err5%_80k:80480"
  "out_sequence_err5%_100k:100480"
  "out_trim_err5%_10k:0"
  "out_trim_err5%_20k:1"
  "out_trim_err5%_40k:2"
  "out_trim_err5%_60k:3"
  "out_trim_err5%_80k:4"
  "out_trim_err5%_100k:5"
  "pbsim3_trim_err5%_10k:6"
  "pbsim3_trim_err5%_20k:7"
  "pbsim3_trim_err5%_40k:8"
  "pbsim3_trim_err5%_60k:9"
  "pbsim3_trim_err5%_80k:10"
  "pbsim3_trim_err5%_100k:11"
)


for item in "${CASES[@]}"; do
  IFS=: read -r DIR LEN <<<"$item"
  echo "Ensuring directory: $DIR"
  mkdir -p "$DIR"
done



for dir in ./out_sequence_err5%*; do
    # Check if it's a directory
    if [ -d "$dir" ]; then
        echo "Clearing contents of directory: $dir"
        rm -rf "$dir"/*  # Remove all files inside the directory but keep the directory
    fi
done



# 10240 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 10240 --prefix ./out_sequence_err5%_10k/chr1_pbsim3_10k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 10240 --prefix ./out_sequence_err5%_10k/chr22_pbsim3_10k

gzip -d ./out_sequence_err5%_10k/chr1_pbsim3_10k_0001.fq.gz
gzip -d ./out_sequence_err5%_10k/chr1_pbsim3_10k_0001.maf.gz

gzip -d ./out_sequence_err5%_10k/chr22_pbsim3_10k_0001.fq.gz
gzip -d ./out_sequence_err5%_10k/chr22_pbsim3_10k_0001.maf.gz


# 20480 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 20480 --prefix ./out_sequence_err5%_20k/chr1_pbsim3_20k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 20480 --prefix ./out_sequence_err5%_20k/chr22_pbsim3_20k

gzip -d ./out_sequence_err5%_20k/chr1_pbsim3_20k_0001.fq.gz
gzip -d ./out_sequence_err5%_20k/chr1_pbsim3_20k_0001.maf.gz

gzip -d ./out_sequence_err5%_20k/chr22_pbsim3_20k_0001.fq.gz
gzip -d ./out_sequence_err5%_20k/chr22_pbsim3_20k_0001.maf.gz



# 40480 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 40480 --prefix ./out_sequence_err5%_40k/chr1_pbsim3_40k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 40480 --prefix ./out_sequence_err5%_40k/chr22_pbsim3_40k

gzip -d ./out_sequence_err5%_40k/chr1_pbsim3_40k_0001.fq.gz
gzip -d ./out_sequence_err5%_40k/chr1_pbsim3_40k_0001.maf.gz

gzip -d ./out_sequence_err5%_40k/chr22_pbsim3_40k_0001.fq.gz
gzip -d ./out_sequence_err5%_40k/chr22_pbsim3_40k_0001.maf.gz




# 60480 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 60480 --prefix ./out_sequence_err5%_60k/chr1_pbsim3_60k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 60480 --prefix ./out_sequence_err5%_60k/chr22_pbsim3_60k

gzip -d ./out_sequence_err5%_60k/chr1_pbsim3_60k_0001.fq.gz
gzip -d ./out_sequence_err5%_60k/chr1_pbsim3_60k_0001.maf.gz

gzip -d ./out_sequence_err5%_60k/chr22_pbsim3_60k_0001.fq.gz
gzip -d ./out_sequence_err5%_60k/chr22_pbsim3_60k_0001.maf.gz





# 80480 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 80480 --prefix ./out_sequence_err5%_80k/chr1_pbsim3_80k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 80480 --prefix ./out_sequence_err5%_80k/chr22_pbsim3_80k

gzip -d ./out_sequence_err5%_80k/chr1_pbsim3_80k_0001.fq.gz
gzip -d ./out_sequence_err5%_80k/chr1_pbsim3_80k_0001.maf.gz

gzip -d ./out_sequence_err5%_80k/chr22_pbsim3_80k_0001.fq.gz
gzip -d ./out_sequence_err5%_80k/chr22_pbsim3_80k_0001.maf.gz






# 102400 and error rate = 5%
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr1.fasta   --length-sd 0 --accuracy-mean 0.95 --length-mean 100480 --prefix ./out_sequence_err5%_100k/chr1_pbsim3_100k
../pbsim3/src/pbsim --strategy wgs --method qshmm --qshmm ../pbsim3/data/QSHMM-RSII.model --depth 1 --genome ./Genome_data/GRCh38_chr22.fasta  --length-sd 0 --accuracy-mean 0.95 --length-mean 100480 --prefix ./out_sequence_err5%_100k/chr22_pbsim3_100k

gzip -d ./out_sequence_err5%_100k/chr1_pbsim3_100k_0001.fq.gz
gzip -d ./out_sequence_err5%_100k/chr1_pbsim3_100k_0001.maf.gz

gzip -d ./out_sequence_err5%_100k/chr22_pbsim3_100k_0001.fq.gz
gzip -d ./out_sequence_err5%_100k/chr22_pbsim3_100k_0001.maf.gz



