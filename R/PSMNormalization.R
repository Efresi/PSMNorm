PSMnormalization <- function(targetPeptide_name, targetProtein_name, out_folder) {

  library(stringi)
  library(openxlsx)
  library(readr)

  #Carico file TargetPeptideSpectrumMatch
  targetPeptide <- read_delim(targetPeptide_name , delim = "\t", show_col_types = FALSE, progress = FALSE)

  #Carico TargetProtein
  targetProtein <- read_delim(targetProtein_name, delim = "\t", show_col_types = FALSE, progress = FALSE)


  #Calcolo PSMtot
  PSMtot <- nrow(targetPeptide)

  #Colonne dello sheet Master
  columns <- c("Accession", "Description", "D", "Coverage", "# Peptides",
               "# PSMs", "# Unique Peptides",
               "# Protein Groups", "# AAs","MW [kDa]", "calc. pI", 	"Score Sequest HT",
               "PSM/AA"	, "/PSM tot", "PSM tot", "MS/MS",	"Somma")

  #Creo il file per lo sheet Master
  if(!('D' %in% colnames(targetProtein))){ columns <- columns[-which(columns == 'D')]}
  Master <- data.frame(matrix(nrow = nrow(targetProtein), ncol = length(columns)))
  colnames(Master)<- columns
  Master[, 1:which(columns=='Score Sequest HT')] <- targetProtein[, columns[1:which(columns=='Score Sequest HT')]]

  #-----------------------------------------------------------------------------
  # Immunoglobulin lambda-like polypeptide
  # Gestisco le Immunoglobulin lambda-like polypeptide: se presenti,
  # vado a vedere quante volte compare la sequenza VTVLGQPK nel file TargetPeptide
  # e sottraggo quel valore al #PSMs della lambda-like

  prot <- c("Immunoglobulin lambda-like")
  indices_lambda <- unlist(lapply(prot, function(x) grep(x, Master$Description)))
  check_lamb_like <- 0

  if(length(indices_lambda) != 0){ #se presente la lambda-like

    IDXpeptides <- unlist(lapply("VTVLGQPK", function(x) grep(x, targetPeptide$'Annotated Sequence'))) #cerco la sequenza
    check_lamb_like <- ifelse(length(IDXpeptides) != 0, 1, 0) #se presente la sequenza metti check = 1 per colorare la riga di rosso

    # Correggo psm della lamda like sottraendo il numero di volte che trovo la sequenza VTVLGQPK
    # Correggo #AAs che deve esssere 106 nella lambda-like

    Master[indices_lambda, '# PSMs'] <- Master[indices_lambda, '# PSMs'] - length(IDXpeptides)
    Master[indices_lambda, '# AAs'] <- 106
  }

  #-----------------------------------------------------------------------------
  #Calcolo PSM/AA
  Master[, "PSM/AA"] <- ((Master$`# PSMs`)/(Master$`# AAs`)*100)

  #Calcolo PSM/AA/PSM tot
  PSM_norm <- (Master$`PSM/AA`)/(PSMtot)*100
  Master[, "/PSM tot"] <- PSM_norm

  #Metto in ordine per psm normalizzato
  Master <- Master[order(- Master$`/PSM tot`), ]
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
               "Transthyretin", "Vitronectin", "heavy chain_V")

  idx <- unlist(lapply(amyPROT, function(x) grep(x, Master$Description)))
  Amyloid <- Master[idx, -c(which(colnames(Master)=='D'), which(colnames(Master) == 'Coverage'))]
  Amyloid <- Amyloid[order(- Amyloid$`/PSM tot`), ]
  Amyloid[1, "PSM tot"] <- PSMtot
  Amyloid[1, "Somma"] <- sum(Master[AmySigIndex, "/PSM tot"])

  ##################################################################################
  #----------------------Gestisco file xlsx in uscita--------------------------#

  tmp <- stri_split_boundaries(targetPeptide_name)[[1]][3] #prendo la terza sottostringa
  tmp <- stri_replace_all_fixed(tmp, " ", "")
  nomeFile <- stri_c(out_folder, tmp, ".xlsx")
  list_sheets <- list("Master" = Master, "Amyloid" = Amyloid)
  write.xlsx(list_sheets, nomeFile)

  ####Crea la formattazione in automatico per lo sheet  Amyloid ####
  file <- loadWorkbook(nomeFile)

  # Select the Amyloid sheet
  sheet <- "Amyloid"

  #Metto l'intestazione di colonna in bold
  fill_style <- createStyle(textDecoration = 'bold', border = 'bottom')
  addStyle(file, sheet = 1, rows = 1, cols = 1:ncol(Master), style = fill_style, gridExpand = TRUE)
  addStyle(file, sheet = 2, rows = 1, cols = 1:ncol(Amyloid), style = fill_style, gridExpand = TRUE)

  # Metto il testo della colonna #PSMs in rosso
  index <- which(colnames(Amyloid) == '# PSMs')
  fill_style_PSM <- createStyle(fontColour = "red")
  addStyle(file, sheet = sheet, rows = 1:(nrow(Amyloid)+1), cols = index, style = fill_style_PSM, stack = TRUE)

  #-----------------------------------------------------------------------------
  ####Crea la formattazione in automatico per le proteine dell'Amyloid signature nel foglio Amyloid####

  #Trovo le righe di interesse
  AmySig <- c("Apolipoprotein E",
              "Apolipoprotein A-IV",
              "Serum amyloid P")

  AmySigIndex <- unlist(lapply(AmySig, function(x) grep(x, Amyloid$Description)))
  fill_style <- createStyle(fgFill = rgb(255,255,0, maxColorValue = 255)) #yellow
  addStyle(file, sheet = sheet, rows = AmySigIndex + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

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

  fill_style <- createStyle(fgFill = rgb(218, 238, 243, maxColorValue = 255)) #azzurro
  addStyle(file, sheet = sheet, rows = indices + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)


  #-----------------------------------------------------------------------------

  # Se presente la sequenza VTVLGQPK  colora il testo della riga lambda like di rosso

  if (check_lamb_like == 1){

    idx_l_like <- unlist(lapply("Immunoglobulin lambda-like", function(x) grep(x, Amyloid$Description)))
    fill_style <- createStyle(fontColour = "red")
    addStyle(file, sheet = sheet, rows = idx_l_like + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

  }
  #### TTR e SAA
  #-----------------------------------------------------------------------------
  #Transthyretin OS=Homo sapiens OX=9606 GN=TTR PE=1 SV=1
  #Serum Amyloid A (1,2,4) protein
  #delle SAA evidenzio solo la riga della prot con #PSMs max
  #-----------------------------------------------------------------------------
  prot_TTR <- c("Transthyretin")
  prot_SAA <- c("Serum amyloid A")
  prot <- c("Apolipoprotein A-I ", 'Gelsolin')

  indices_TTR <- unlist(lapply(prot_TTR, function(x) grep(x, Amyloid$Description)))

  indices_SAA <- unlist(lapply(prot_SAA, function(x) grep(x, Amyloid$Description)))

  indices <- unlist(lapply(prot, function(x) grep(x, Amyloid$Description)))

  # Find the indices of SAA with the maximum value of Amyloid$#PSMs
  max_index <- which.max(Amyloid$`# PSMs`[indices_SAA])

  # Select only the index with the maximum value of Amyloid$PSMs
  indices_SAA_max <- indices_SAA[max_index]

  indices <- c(indices_TTR, indices_SAA_max, indices)

  fill_style <- createStyle(fgFill = rgb(242, 220, 219, maxColorValue = 255)) #rosso tenue
  addStyle(file, sheet = sheet, rows = indices + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)


  #-----------------------------------------------------------------------------
  ##### Fibrinogeno
  #-----------------------------------------------------------------------------
  #Fibrinogen alpha chain OS=Homo sapiens OX=9606 GN=FGA PE=1 SV=2
  #Fibrinogen beta chain OS=Homo sapiens OX=9606 GN=FGB PE=1 SV=2
  #Fibrinogen gamma chain OS=Homo sapiens OX=9606 GN=FGG PE=1 SV=3
  #-----------------------------------------------------------------------------
  prot <- c("Fibrinogen")
  indices <- unlist(lapply(prot, function(x) grep(x, Amyloid$Description)))

  fill_style <- createStyle(fgFill = rgb(228, 223, 236, maxColorValue = 255)) #lilla
  addStyle(file, sheet = sheet, rows = indices + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

  #-----------------------------------------------------------------------------
  ##### Serum Albumin and Vitronectin
  #-----------------------------------------------------------------------------
  #Serum albumin OS=Homo sapiens GN=ALB PE=1 SV=2
  #Vitronectin OS=Homo sapiens OX=9606 GN=VTN PE=1 SV=1
  #-----------------------------------------------------------------------------
  prot <- c("Serum albumin", 'Vitronectin')
  indices <- unlist(lapply(prot, function(x) grep(x, Amyloid$Description)))

  fill_style <- createStyle(fgFill = rgb(253, 233, 217, maxColorValue = 255)) #arancio
  addStyle(file, sheet = sheet, rows = indices + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

  #-----------------------------------------------------------------------------
  #### heavy chain_V
  #-----------------------------------------------------------------------------
  # Quando è presente heavy chain_V, check se presente la sequenza GLEWVSAISGSGGGTYYADSVK
  # nel file targetpeptide

  prot <- c("heavy chain_V")
  indices <- unlist(lapply(prot, function(x) grep(x, Amyloid$Description)))
  check_heavy_chainV <- 0

  if(length(indices) != 0){
    IDXpeptides <- unlist(lapply("GLEWVSAISGSGGGTYYADSVK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
    check_heavy_chainV <- ifelse(length(IDXpeptides) != 0, 1, 0)
  }

  # Se presente la sequenza GLEWVSAISGSGGGTYYADSVK  colora il testo della riga di verde

  if (check_heavy_chainV == 1){

    fill_style <- createStyle(fgFill = rgb(235,241,222, maxColorValue = 255))#verde
    addStyle(file, sheet = sheet, rows = indices + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

  }

  #-----------------------------------------------------------------------------
  #### Apolipoprotein A-IV
  #-----------------------------------------------------------------------------
  # Gestisco la Apolipoprotein A-IV, check se presente almeno una delle seguenti sequenze
  # - ARAEVSADQVATVMWDYFSQLSNNAK
  # - RAEVSADQVATVMWDYFSQLSNNAK
  # - AEVSADQVATVMWDYFSQLSNNAK

  check_apo_A_IV <- 0

  if(length(indices) != 0){
    IDXpeptides_1 <- unlist(lapply("ARAEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
    IDXpeptides_2 <- unlist(lapply("RAEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
    IDXpeptides_3 <- unlist(lapply("AEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
    cond = length(IDXpeptides_1) != 0 | length(IDXpeptides_2) != 0 | length(IDXpeptides_3) != 0
    check_apo_A_IV <- ifelse(cond, 1, 0)
  }

  # Se presente check==1colora il testo della riga lambda like di rosso

  if (check_apo_A_IV == 1){

    idx <- unlist(lapply("Apolipoprotein A-IV", function(x) grep(x, Amyloid$Description)))
    fill_style <- createStyle(fontColour = "red")
    addStyle(file, sheet = sheet, rows = idx + 1, cols = 1:which(colnames(Amyloid) == '/PSM tot'), style = fill_style, stack = TRUE, gridExpand = TRUE)

  }

  #-----------------------------------------------------------------------------

  setColWidths(file, 'Master', widths = 'auto', cols = 1:20)
  setColWidths(file, 'Amyloid', widths = 'auto', cols = 1:20)
  saveWorkbook(file, nomeFile, overwrite = TRUE)

}
