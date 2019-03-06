ehei <- function(otuTable, taxonomyCol, sampleCols){
  library(data.table)
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

  # coerce x to data.table class
  otuTable <- as.data.table(otuTable)

  # detect row indices of target genera in our OTU table
  hcDegOtus <- grep(paste(hcDegraders, collapse = "|"),
    as.character(otuTable[[taxonomyCol]]), ignore.case = T)

  # calculate total abundance of all target OTUs in each sample
  hcDegAbunds <- as.data.table(otuTable[
    hcDegOtus, colSums(.SD), .SDcols = sampleCols],
    keep.rownames = T)

  # get total library sizes
  sampSizes <- as.data.table(otuTable[, colSums(.SD), .SDcols = sampleCols],
    keep.rownames = T)

  # merge total abundances with total library sizes
  hcDegAbunds[sampSizes, librarySize := i.V2, on = c(V1 = "V1")]

  # adjust names of cols
  names(hcDegAbunds)[1:2] <- c("sample", "hdAbund")

  # calculate total rel abundance of all hc degraders
  hcDegAbunds[, exposure := hdAbund/librarySize]

  # delete unnecessary cols
  hcDegAbunds[, c("hdAbund", "librarySize") := NULL]

  # return data.table with sample names and exposure indices
  return(hcDegAbunds)
}
