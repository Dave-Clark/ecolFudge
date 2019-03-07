ehei <- function(otuTable, taxonomyCol, sampleCols){
  # print message telling user how many samples are being used
  message("Calculating Ecological Hydrocarbon exposure index for ",
    length(sampleCols), " sample(s)")

  # list of target genera
  hcDegraders <- c("Arthrobacter",
    "Dietzia",
    "Gordonia",
    "Microbacterium",
    "Micrococcus",
    "Mycobacterium",
    "Nocardia",
    "Nocardioides",
    "Prauserella",
    "Rhodococcus",
    "Streptomyces",
    "Terrabacter",
    "Cytophaga",
    "Flavobacterium",
    "Pedobacter",
    "Bacillus",
    "Paenibacillus",
    "Planomicrobium",
    "Kordiimonas",
    "Novosphingobium",
    "Ochrobactrum",
    "Sphingobium",
    "Sphingomonas",
    "Sphingopyxis",
    "Roseobacter_clade",
    "Roseovarius",
    "Jannaschia",
    "Silicibacter",
    "Sulfitobacter",
    "Thalassospira",
    "Tranquillimonas",
    "Tropicibacter",
    "Tropicimonas",
    "Acidovorax",
    "Alcaligenes",
    "Burkholderia",
    "Comamonas",
    "Delftia",
    "Polaromonas",
    "Ralstonia",
    "Desulfatibacillum",
    "Desulfatiferula",
    "Desulfobacterium",
    "Desulfococcus",
    "Desulfoglaeba",
    "Desulfothermus",
    "Alcanivorax",
    "Alkanindiges",
    "Cycloclasticus",
    "Oleibacter",
    "Oleiphilus",
    "Oleispira",
    "Thalassolituus",
    "Acinetobacter",
    "Alteromonas",
    "Halomonas",
    "Marinobacter",
    "Microbulbifer",
    "Neptunomonas",
    "Pseudoalteromonas",
    "Pseudomonas",
    "Shewanella",
    "Vibrio")

  # detect row indices of target genera in our OTU table
  hcDegOtus <- grep(paste(hcDegraders, collapse = "|"),
    as.character(otuTable[, taxonomyCol]), ignore.case = T)

  # calculate total abundance of all target OTUs in each sample
  hcDegAbunds <- data.frame(sample = sampleCols,
    hcAbund = colSums(otuTable[hcDegOtus, sampleCols]),
    librarySize = colSums(otuTable[, sampleCols]))

  # calculate total rel abundance of all hc degraders
  hcDegAbunds$exposure <- hcDegAbunds$hcAbund/hcDegAbunds$librarySize

  # delete unnecessary cols
  hcDegAbunds[, c("hcAbund", "librarySize")] <- NULL

  # return data.frame with sample names and exposure indices
  return(hcDegAbunds)
}
