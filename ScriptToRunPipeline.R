#-------------------------------------------------------------------------------
# Cross-tissue macrophage profiling in Sus scrofa (Using PigGTEx data)
# Run from RStudio Desktop (Windows), executing via WSL
#-------------------------------------------------------------------------------

library(Seurat)

# ---- CONFIG ------------------------------------------------------------------
# these are paths INSIDE WSL, not Windows paths.

config <- list(
  genome_dir = "~/pig_genome",
  index_dir  = "~/pig_genome/index",
  fastq_dir  = "~/data/fastq",
  qc_dir     = "~/qc_reports",
  whitelist  = "~/pig_genome/10x_whitelist.txt",
  threads    = 8,
  samples    = c("SRR24723931")
)

# run a command inside WSL from Windows R
wsl_run <- function(cmd) {
  full_cmd <- paste("wsl bash -c", shQuote(cmd))
  system(full_cmd)
}

# ---- FUNCTIONS ---------------------------------------------------------------

run_fastqc <- function(sample, cfg) {
  r1 <- file.path(cfg$fastq_dir, paste0(sample, "_1.fastq.gz"))
  r2 <- file.path(cfg$fastq_dir, paste0(sample, "_2.fastq.gz"))
  wsl_run(paste("fastqc", r1, r2, "-o", cfg$qc_dir))
}

run_trim <- function(sample, cfg) {
  r2      <- file.path(cfg$fastq_dir, paste0(sample, "_2.fastq.gz"))
  r2_trim <- file.path(cfg$fastq_dir, paste0(sample, "_2.trimmed.fastq.gz"))
  html    <- file.path(cfg$qc_dir, paste0(sample, "_fastp.html"))
  wsl_run(paste("fastp -i", r2, "-o", r2_trim, "--html", html))
}

build_genome_index <- function(cfg) {
  check_cmd <- paste0("test -d ", cfg$index_dir, " && echo EXISTS")
  exists <- system(paste("wsl bash -c", shQuote(check_cmd)), intern = TRUE)
  if (length(exists) > 0 && exists == "EXISTS") {
    message("Genome index already exists, skipping.")
    return(invisible())
  }
  wsl_run(paste("mkdir -p", cfg$index_dir))
  wsl_run(paste(
    "STAR --runMode genomeGenerate",
    "--genomeDir", cfg$index_dir,
    "--genomeFastaFiles", file.path(cfg$genome_dir, "Sus_scrofa.Sscrofa11.1.dna.toplevel.fa"),
    "--sjdbGTFfile", file.path(cfg$genome_dir, "Sus_scrofa.Sscrofa11.1.111.gtf"),
    "--sjdbOverhang 89",
    "--runThreadN", cfg$threads
  ))
}

run_starsolo <- function(sample, cfg) {
  r1      <- file.path(cfg$fastq_dir, paste0(sample, "_1.fastq.gz"))
  r2_trim <- file.path(cfg$fastq_dir, paste0(sample, "_2.trimmed.fastq.gz"))
  wsl_run(paste(
    "STAR --runMode alignReads",
    "--genomeDir", cfg$index_dir,
    "--readFilesIn", r2_trim, r1,
    "--readFilesCommand zcat",
    "--soloType CB_UMI_Simple",
    "--soloCBwhitelist", cfg$whitelist,
    "--soloCBstart 1 --soloCBlen 16 --soloUMIstart 17 --soloUMIlen 10",
    "--outSAMtype BAM SortedByCoordinate",
    "--runThreadN", cfg$threads,
    "--outFileNamePrefix", paste0(sample, "_")
  ))
}

# Read10X needs a Windows-visible path, since it's Read10X() 

load_into_seurat <- function(sample, wsl_distro = "Ubuntu") {
  wsl_path <- sprintf("//wsl$/%s/home/%s/%s_Solo.out/Gene/filtered",
                      wsl_distro, Sys.getenv("USER"), sample)
  counts <- Read10X(data.dir = wsl_path)
  CreateSeuratObject(counts = counts, project = sample)
}

# ---- MAIN PIPELINE -----------------------------------------------------------

main <- function(cfg) {
  wsl_run(paste("mkdir -p", cfg$qc_dir))
  build_genome_index(cfg)
  
  seurat_objects <- list()
  for (s in cfg$samples) {
    run_fastqc(s, cfg)
    run_trim(s, cfg)
    run_starsolo(s, cfg)
    seurat_objects[[s]] <- load_into_seurat(s)
  }
  seurat_objects
}

# ---- RUN ---------------------------------------------------------------------
seurat_objects <- main(config)
saveRDS(seurat_objects, "seurat_objects_raw.rds")