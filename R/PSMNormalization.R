PSMnormalization <- function(targetPeptide_name, targetProtein_name) {

  library(data.table)
  library(stringi)

  #Carico file TargetPeptideSpectrumMatch
  targetPeptide <- read.csv(targetPeptide_name , sep = "" )

  #Carico TargetProtein
  targetProtein <- read.csv(targetProtein_name, sep = "" )

  colnames(targetProtein) <- c("Checked", "Master", "Accession", "Description", "Contaminant", "D",
                               "Coverage", "# PSMs", "Score Sequest HT", "# Peptides", "# Unique Peptides",
                               "# Protein Groups", "# AAs", "MW [kDa]", "calc. pI", "Modifications",
                               "# Peptides Sequest HT")

  #Calcolo PSMtot
  PSMtot <- nrow(targetPeptide)

  columns <- c("Accession", "Description", "D", "Coverage", "# Peptides",
               "# PSMs", "# Unique Peptides",
               "# Protein Groups", "# AAs","MW [kDa]", "calc. pI", 	"Score Sequest HT",
               "PSM/AA"	, "/PSM tot", "PSM tot", "MS/MS",	"Somma")

  #Creo il file per lo sheet Master
  Master <- data.frame(matrix(nrow = nrow(targetProtein), ncol = length(columns)))
  colnames(Master)<- columns
  Master[, 1:which(columns=='Score Sequest HT')] <- targetProtein[, columns[1:which(columns=='Score Sequest HT')]]

  #Gestisco le Immunoglobulin lambda-like polypeptide: se presenti,
  #vado a vedere quante volte compare la sequenza VTVLGQPK nel file TargetPeptide
  #e sottraggo quel valore al #PSMs della lambda-like

  prot <- c("Immunoglobulin lambda-like")
  indices_lambda <- unlist(lapply(prot, function(x) grep(x, Master$Description)))
  check_lamb_like <- 0
  if(length(indices_lambda) != 0){
    IDXpeptides <- unlist(lapply("VTVLGQPK", function(x) grep(x, targetPeptide$Annotated.Sequence)))
    check_lamb_like <- 1
  }
  Master[indices_lambda, '# PSMs'] <- Master[indices_lambda, '# PSMs'] - length(IDXpeptides)
  Master[indices_lambda, '# AAs'] <- 106

  #Calcolo PSM/AA
  Master[, "PSM/AA"] <- round(((Master$`# PSMs`)/(Master$`# AAs`)*100), 3)

  #Calcolo PSM/AA/PSM tot
  PSM_norm <- (Master$`PSM/AA`)/(PSMtot)*100
  Master[, "/PSM tot"] <- round(PSM_norm, 3)

  #Metto in ordine per psm normalizzato
  Master <- setorderv(Master, cols = "/PSM tot", order = -1)
  Master[1, "PSM tot"] <- PSMtot

  #Somma le proteine dell'Amyloid signature
  #-----------------------------------------------------------------------------
  #Apolipoprotein E OS=Homo sapiens OX=9606 GN=APOE PE=1 SV=1
  #Apolipoprotein A-IV OS=Homo sapiens OX=9606 GN=APOA4 PE=1 SV=3
  #Serum amyloid P-component OS=Homo sapiens OX=9606 GN=APCS PE=1 SV=2
  #-----------------------------------------------------------------------------
  AmySig <- c("Apolipoprotein E",
              "Apolipoprotein A-IV",
              "Serum amyloid P")

  AmySigIndex <- unlist(lapply(AmySig, function(x) grep(x, Master$Description)))

  Master[1, "Somma"] <- sum(Master[AmySigIndex, "/PSM tot"])

  #Seleziono proteine per foglio amyloid
  amyPROT <- c("Apolipoprotein A", "Apolipoprotein C", "Apolipoprotein E", "Beta-2-microglobulin",
               "Fibrinogen", "Gelsolin", "Immunoglobulin", "Lysozyme", "Serum albumin", "Serum amyloid",
               "Transthyretin", "Vitronectin")
  idx <- unlist(lapply(amyPROT, function(x) grep(x, Master$Description)))
  Amyloid <- Master[idx, -c(which(colnames(Master)=='D'), which(colnames(Master)=='Coverage'))]
  Amyloid <- setorderv(Amyloid, cols = "/PSM tot", order = -1)
  Master <- Master[, -which(colnames(Master)=="Somma")]

  ##################################################################################
  #----------------------Gestisco file xlsx in uscita--------------------------#
  library(xlsx)

  tmp <- stri_split_boundaries(targetPeptide_name)[[1]][3] #prendo la terza sottostringa
  tmp <- stri_replace_all_fixed(tmp, " ", "")
  nomeFile <- stri_c(tmp, ".xlsx")
  write.xlsx(Master, nomeFile, sheetName = "Master", showNA = FALSE, row.names = FALSE)
  write.xlsx(Amyloid, nomeFile, sheetName = "Amyloid", showNA = FALSE, row.names = FALSE, append = TRUE)

  ####Crea la formattazione in automatico per lo sheet  Amyloid ####
  file <- loadWorkbook(nomeFile)

  #Seleziono il foglio del file excel (Amyloid, sheet numero 2)
  AllList <- getSheets(file)
  sheet <- AllList[[2]]
  #-----------------------------------------------------------------------------
  # Metto il testo della colonna #PSMs in rosso

  index <- which(colnames(Amyloid) == '# PSMs')

  cells <- getCells(getRows(sheet), colIndex = index)

  fill_style_PSM <- CellStyle(file) + Font(file, color = 'red')

  #Applying color to the column
  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_style_PSM)
  }
  #-----------------------------------------------------------------------------
  ####Crea la formattazione in automatico per le proteine dell'Amyloid signature nel foglio Amyloid####

  #Trovo le righe di interesse
  AmySig <- c("Apolipoprotein E",
              "Apolipoprotein A-IV",
              "Serum amyloid P")

  AmySigIndex <- unlist(lapply(AmySig, function(x) grep(x, Amyloid$Description)))
  row <- getRows(sheet, AmySigIndex+1)
  cells <- getCells(row, 1:which(colnames(Amyloid)=='/PSM tot'))

  #Specifico il colore delle righe
  fill_style <- CellStyle(file) +

    Fill(backgroundColor = "yellow",

         foregroundColor = "yellow",

         pattern         = "SOLID_FOREGROUND")

  #Applico il colore alle righe
  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_style)
  }
  #-----------------------------------------------------------------------------
  #Tengo il font della colonna #PSMs rosso
  cells <- getCells(row, index)
  fill_new <- fill_style + Font(file, color = 'red')
  fill_new$font <- Font(file, color = 'red')

  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_new)
  }
  #-----------------------------------------------------------------------------
  ####------Crea la formattazione in automatico per le altre proteine------####

  #### IMMUNOGLOBULINE
  #-----------------------------------------------------------------------------
  # Immunoglobulin lambda constant
  # Immunoglobulin kappa constant
  # Immunoglobulin lambda-like polypeptide
  #-----------------------------------------------------------------------------

  prot_lambda <- c("Immunoglobulin lambda constant", "Immunoglobulin lambda-like")
  prot_kappa <- "Immunoglobulin kappa constant"

  indices_lambda <- unlist(lapply(prot_lambda, function(x) grep(x, Amyloid$Description)))

  # Find the indices of lambda with the maximum value of Amyloid$#PSMs
  max_index <- which.max(Amyloid$`# PSMs`[indices_lambda])

  # Select only the index with the maximum value of Amyloid$PSMs
  indices_lambda_max <- indices_lambda[max_index]

  indices_kappa <- unlist(lapply(prot_kappa, function(x) grep(x, Amyloid$Description)))

  indices <- c(indices_kappa, indices_lambda_max)

  row <- getRows(sheet, indices+1)
  cells <- getCells(row, 1:which(colnames(Amyloid)=='/PSM tot'))

  #Specifico il colore delle righe
  fill_style <- CellStyle(file) +

    Fill(backgroundColor = "lightblue",

         foregroundColor = "lightblue",

         pattern         = "SOLID_FOREGROUND")

  #Applying color to the row
  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_style)
  }
  #-----------------------------------------------------------------------------
  #Tengo il font della colonna #PSMs rosso
  cells <- getCells(row, index)
  fill_new <- fill_style + Font(file, color = 'red')
  fill_new$font <- Font(file, color = 'red')

  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_new)
  }
  #-----------------------------------------------------------------------------

  # Se presente la sequenza VTVLGQPK  colora il testo della riga lambda like di rosso

  if (check_lamb_like == 1){

    idx_l_like <- unlist(lapply("Immunoglobulin lambda-like", function(x) grep(x, Amyloid$Description)))
    row <- getRows(sheet, idx_l_like+1)
    cells <- getCells(row, 1:which(colnames(Amyloid)=='/PSM tot'))

    if (idx_l_like == indices_lambda_max){
      #se lambda-like è la proteina con # PSMs maggiore

      fill_style <- CellStyle(file) +
        Fill(backgroundColor = "lightblue",
             foregroundColor = "lightblue",
             pattern         = "SOLID_FOREGROUND") +
        Font(file, color = 'red')

    } else {
      fill_style <- CellStyle(file) +
        Font(file, color = 'red')
    }

    #Applying color to the rows
    for (i in names(cells)) {
      setCellStyle(cells[[i]], fill_style)
    }
  }
  #### TTR e SAA
  #-----------------------------------------------------------------------------
  #Transthyretin OS=Homo sapiens OX=9606 GN=TTR PE=1 SV=1
  #Serum Amyloid A (1,2,4) protein
  #delle SAA coloro solo la riga della prot con #PSMs max
  #-----------------------------------------------------------------------------
  prot_TTR <- c("Transthyretin")
  prot_SAA <- c("Serum amyloid A")

  indices_TTR <- unlist(lapply(prot_TTR, function(x) grep(x, Amyloid$Description)))

  indices_SAA <- unlist(lapply(prot_SAA, function(x) grep(x, Amyloid$Description)))

  # Find the indices of SAA with the maximum value of Amyloid$#PSMs
  max_index <- which.max(Amyloid$`# PSMs`[indices_SAA])

  # Select only the index with the maximum value of Amyloid$PSMs
  indices_SAA_max <- indices_SAA[max_index]


  indices <- c(indices_TTR, indices_SAA_max)

  row <- getRows(sheet, indices+1)
  cells <- getCells(row, 1:which(colnames(Amyloid)=='/PSM tot'))
  #Specifico il colore delle righe
  fill_style <- CellStyle(file) +

    Fill(backgroundColor = "pink3",

         foregroundColor = "pink3",

         pattern         = "SOLID_FOREGROUND")

  #Applying color to the row
  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_style)
  }
  #-----------------------------------------------------------------------------
  #Tengo il font della colonna #PSMs rosso
  cells <- getCells(row, index)
  fill_new <- fill_style + Font(file, color = 'red')
  fill_new$font <- Font(file, color = 'red')

  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_new)
  }
  #-----------------------------------------------------------------------------
  ##### Fibrinogeno
  #-----------------------------------------------------------------------------
  #Fibrinogen alpha chain OS=Homo sapiens OX=9606 GN=FGA PE=1 SV=2
  #Fibrinogen beta chain OS=Homo sapiens OX=9606 GN=FGB PE=1 SV=2
  #Fibrinogen gamma chain OS=Homo sapiens OX=9606 GN=FGG PE=1 SV=3
  #-----------------------------------------------------------------------------
  prot <- c("Fibrinogen")
  indices <- unlist(lapply(prot, function(x) grep(x, Amyloid$Description)))
  row <- getRows(sheet, indices+1)

  cells <- getCells(row, 1:which(colnames(Amyloid)=='/PSM tot'))
  #Specifico il colore delle righe
  fill_style <- CellStyle(file) +

    Fill(backgroundColor = "peachpuff",

         foregroundColor = "peachpuff",

         pattern         = "SOLID_FOREGROUND")

  #Applying color to the row
  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_style)
  }
  #-----------------------------------------------------------------------------
  #Tengo il font della colonna #PSMs rosso
  cells <- getCells(row, index)
  fill_new <- fill_style + Font(file, color = 'red')
  fill_new$font <- Font(file, color = 'red')

  for (i in names(cells)) {
    setCellStyle(cells[[i]], fill_new)
  }
  #-----------------------------------------------------------------------------


  autoSizeColumn(getSheets(file)[[1]], colIndex = 1:20)
  autoSizeColumn(getSheets(file)[[2]], colIndex = 1:20)
  saveWorkbook(file, nomeFile)

}
