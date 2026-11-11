install.packages("BiocManager")
BiocManager::install(c("Biostrings", "msa"))
install.packages("base64enc")
install.packages(c("rphylopic", "rentrez"))

library(rphylopic)
library(rentrez)
library(Biostrings)
library(msa)

fetch_sequences <- function(species_list, db, search_term, image_dir, fasta_file) {
  dir.create(image_dir, showWarnings = FALSE)
  file.create(fasta_file)
  
  for (specie in species_list) {
    safe_name <- gsub(" ", "_", specie)
    
    tryCatch({
      uuid <- get_uuid(name = specie)
      img <- get_phylopic(uuid)
      save_phylopic(img, path = paste0(image_dir, "/", safe_name, ".png"))
      cat("Saved image for:", specie, "\n")
    }, error = function(e) {
      genus <- strsplit(specie, " ")[[1]][1]
      tryCatch({
        uuid <- get_uuid(name = genus)
        img <- get_phylopic(uuid)
        save_phylopic(img, path = paste0(image_dir, "/", safe_name, ".png"))
        cat("Used genus-level image for:", specie, "\n")
      }, error = function(e2) {
        cat("No PhyloPic image for:", specie, "\n")
      })
    })
    
    term <- gsub("\\{species\\}", specie, search_term)
    ids <- entrez_search(db = db, term = term, retmax = 1)$ids
    
    if (length(ids) > 0) {
      seqs <- entrez_fetch(db = db, id = ids, rettype = "fasta")
      seqs <- sub("^>.*?\n", paste0(">", specie, "\n"), seqs)
      if (!grepl("\n$", seqs)) seqs <- paste0(seqs, "\n")
      cat(seqs, file = fasta_file, append = TRUE)
      cat("Retrieved sequence for:", specie, "\n")
    } else {
      cat("No sequence found for:", specie, "\n")
    }
  }
}

align_sequences <- function(fasta_file, aligned_file, seq_type = "protein") {
  if (seq_type == "protein") {
    seqs <- readAAStringSet(fasta_file)
    aln <- msa(seqs, type = "protein")
  } else {
    seqs <- readDNAStringSet(fasta_file)
    aln <- msa(seqs, type = "dna")
  }
  aln_seqs <- as(aln, ifelse(seq_type == "protein", "AAStringSet", "DNAStringSet"))
  writeXStringSet(aln_seqs, filepath = aligned_file)
  cat("Alignment saved to:", aligned_file, "\n")
  return(aln)
}

# --- Primate BRCA1 ---
primates <- c(
  "Pan troglodytes",       # Chimpanzee
  "Gorilla gorilla",       # Western Gorilla
  "Pongo pygmaeus",        # Bornean Orangutan
  "Hylobates lar",         # White-handed Gibbon
  "Macaca mulatta",        # Rhesus Macaque
  "Papio anubis",          # Olive Baboon
  "Saimiri sciureus",      # Common Squirrel Monkey
  "Callithrix jacchus",    # Common Marmoset
  "Homo sapiens"           # Human
)

fetch_sequences(
  species_list = primates,
  db = "protein",
  search_term = "BRCA1[Gene] AND {species}[Organism]",
  image_dir = "primate_images",
  fasta_file = "primate_BRCA1.fasta"
)

primate_aln <- align_sequences("primate_BRCA1.fasta", "primate_BRCA1_aligned.fasta", seq_type = "protein")

# --- Mammal MT genomes ---
mammals <- c(
  "Homo sapiens",             # Human
  "Macaca mulatta",           # Rhesus Macaque
  "Ornithorhynchus anatinus", # Platypus
  "Macropus giganteus",       # Eastern Grey Kangaroo
  "Loxodonta africana",       # African Elephant
  "Ursus maritimus",          # Polar Bear
  "Pteropus vampyrus",        # Large Flying Fox (Fruit Bat)
  "Canis lupus",              # Wolf
  "Panthera onca",            # Tiger
  "Balaenoptera musculus",    # Blue Whale
  "Tursiops truncatus",       # Bottlenose Dolphin
  "Mus musculus"              # Mouse
)

fetch_sequences(
  species_list = mammals,
  db = "nucleotide",
  search_term = "{species}[Organism] AND complete genome[Title] AND mitochondrion[Filter]",
  image_dir = "mammal_images",
  fasta_file = "mammal_mt.fasta"
)

mammal_aln <- align_sequences("mammal_mt.fasta", "mammal_mt_aligned.fasta", seq_type = "dna")


# --- Universal SSU rRNA (Woese et al. 1990, Figure 1) ---
ssu_species <- c(
  "Homo sapiens",                  # Animals
  "Drosophila melanogaster",      # Animals
  "Euglena gracilis",             # Euglenozoans
  "Arabidopsis thaliana",         # Green plants
  "Chlamydomonas reinhardtii",    # Green algae
  "Saccharomyces cerevisiae",     # Fungi
  "Neurospora crassa",            # Fungi
  "Tetrahymena thermophila",      # Ciliates
  "Paramecium tetraurelia",       # Ciliates
  "Dictyostelium discoideum",     # Amoebozoans
  "Giardia intestinalis"          # Metamonads
)

fetch_sequences(
  species_list = ssu_species,
  db = "nucleotide",
  search_term = "{species}[Organism] AND (18S ribosomal RNA[Title] OR small subunit ribosomal RNA[Title]) AND 1500:2100[SLEN]",
  image_dir = "ssu_images",
  fasta_file = "ssu_rRNA.fasta"
)

ssu_aln <- align_sequences("ssu_rRNA.fasta", "ssu_rRNA_aligned.fasta", seq_type = "dna")
