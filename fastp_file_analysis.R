library("pipelineTools")

# Path for fastp executable 
fastp.path <- "anaconda3/pkgs/fastp-0.20.0-hd9629dc_0/bin/fastp"

# Read in the fastq file paths to lists and then create sample names lists and trimmed read file paths
reads.path <- "Documents/GitHub/ComputationalEcology/data_analysis_files/edna_data/fastq_raw_reads"
# Trimmed reads directory
trimmed.reads.dir <- "Documents/GitHub/ComputationalEcology/data_analysis_files/edna_data/trimmed_reads"
# FastP results directory
fastp.results.dir <- "Documents/GitHub/ComputationalEcology/data_analysis_files/edna_data/fastpQC"

# Set-up pattern for reads file
reads.patt.1 <- "*_1.fastq.gz$"
reads.patt.2 <- "*_2.fastq.gz$"
sample.dataframe <- prepare_samples(reads.path, c(reads.patt.1,reads.patt.2),trimmed.reads.dir)

# Create lists of filepaths for forward reads
mate1 <- as.character(sample.dataframe$reads.path.1)
mate1.trim <- as.character(sample.dataframe$trimmed.reads.path.1)
# Create lists of filepaths for reverse reads
mate2 <- as.character(sample.dataframe$reads.path.2)
mate2.trim <- as.character(sample.dataframe$trimmed.reads.path.2)

sample.names <- as.character(sample.dataframe$sample.names)

# MiSeq adapter sequences
adapter1 <- "ACACTGACGACATGGTTCTACA"
adapter2 <- "TACGGTAGCAGAGACTTGGTCT"

# Run fast p
fastp.cmds <- run_fastp(mate1 = mate1,
                        mate2 = mate2,
                        mate1.out = mate1.trim,
                        mate2.out = mate2.trim,
                        adapter1 = adapter1,
                        adapter2 = adapter2,
                        sample.name =  sample.names,
                        out.dir = fastp.results.dir,
                        fastp = fastp.path)

write.table(fastp.cmds,"FastP_commands.sh", quote = FALSE, row.names = FALSE, col.names = FALSE)

