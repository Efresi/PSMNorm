appendTabellona <- function(targetPeptide_name, targetProtein_name, Tabellona_PATH, id_pz_escl = NA) {
  library(stringi)
  library(openxlsx)
  library(readr)

  #Carico file TargetPeptideSpectrumMatch
  targetPeptide <- read_delim(targetPeptide_name , delim = "\t", show_col_types = FALSE, progress = FALSE)

  #Carico TargetProtein
  targetProtein <- read_delim(targetProtein_name, delim = "\t", show_col_types = FALSE, progress = FALSE)

  #Rimuovo le proteine Bos taurus dal file targetProtein
  #non  deve essere preso in considerazione nei calcoli

  indices_BosTaurus <- unlist(lapply("Bos taurus", function(x) grep(x, targetProtein$Description)))
  if(length(indices_BosTaurus)!=0){
    targetProtein <- targetProtein[-indices_BosTaurus, ]
  }


  #prendo numero paziente, codice e file dal nome file
  #prendo PSM tot dal numero righe di TargetPeptide

  tmp <- stri_split_fixed(targetPeptide_name, ' ')[[1]] #attenzione

  n_paziente <- as.numeric(tmp[1]) #attenzione
  codice_pz  <- tmp[3] #attenzione

  #Alcuni file hanno il trattino basso, altri il trattino alto, altri hanno lo spazio tra la data
  #e il nome del file. Per prendere il nome del file corretto partiamo da dx verso sx '<-' prendendo
  #la terzultima(-2) e penultima(-1) sottostringa

  file_pz    <- stri_split_fixed(tmp[length(tmp)], "_")[[1]]
  last_item  <- length(file_pz)
  file_pz    <- stri_c(file_pz[(last_item-2)], "_", file_pz[(last_item-1)])

  PSMtot <- nrow(targetPeptide)

  #prendo coverage, #psm e score sequest ht dal file TargetProtein
  amyPROT <- c("Apolipoprotein A-I", "Apolipoprotein C", "Apolipoprotein E", "Beta-2-microglobulin",
               "Fibrinogen", "Gelsolin",
               "Immunoglobulin heavy constant alpha", "Immunoglobulin heavy constant gamma",
               "Immunoglobulin heavy constant mu", "Immunoglobulin kappa constant",
               "Immunoglobulin lambda constant", "Immunoglobulin lambda-like",
               "Lysozyme", "Serum albumin", "Serum amyloid",
               "Transthyretin", "Vitronectin")

  columns <- c("Description","# AAs", "Coverage", "# PSMs", "Score Sequest HT")

  idx <- unlist(lapply(amyPROT, function(x) grep(x, targetProtein$Description)))

  #Check delle colonne mancanti: se al file targetProtein manca qualche colonna, ne creo una fittizia di tutti missing

  for (i in 1:length(columns[1:which(columns=='Score Sequest HT')])){
    if(!(columns[i] %in% colnames(targetProtein))){

      tmp <- rep('ND', nrow(targetProtein))
      targetProtein <- cbind(targetProtein, tmp)
      colnames(targetProtein)[ncol(targetProtein)] <- columns[i]

    }
  }

  #Correggo #AAs e #PSMs per Immunoglobulin lambda-like, se presente

  idx_ll <- unlist(lapply("Immunoglobulin lambda-like", function(x) grep(x, targetProtein$Description)))
  if(length(idx_ll) != 0){

    #Cerco per la sequenza VTVLGQPK
    IDXpeptides <- unlist(lapply("VTVLGQPK", function(x) grep(x, targetPeptide$'Annotated Sequence'))) #cerco la sequenza

    targetProtein$`# PSMs`[idx_ll] <- targetProtein$`# PSMs`[idx_ll] - length(IDXpeptides)
    targetProtein$`# AAs`[idx_ll] <- 106

  }

  #Calcolo PSM/AA
  PSM_AA <- ((targetProtein$`# PSMs`[idx])/(targetProtein$`# AAs`[idx])*100)

  #Calcolo PSM/AA/PSM tot
  PSM_norm <- (PSM_AA)/(PSMtot)*100
  PSM_normfinal <- PSM_norm

  #Creazione del dataframe dato_new_pz

  dato_new_pz <- targetProtein[idx, columns]
  dato_new_pz <- cbind(dato_new_pz, PSM_AA, PSM_normfinal)

  colnames(dato_new_pz) <- c("Description","# AAs", "Coverage", "# PSMs", "Score Sequest HT", "PSM/AA", "/PSM tot")


  #Lettura di Tabellona
  Tabellona <- read.xlsx(Tabellona_PATH, sheet = 1, sep.names = " ")

  idx_prot <- c('Apolipoprotein A-I ' , 'Apolipoprotein A-II ' , 'Apolipoprotein A-IV ' ,
                'Apolipoprotein C-II ' , 'Apolipoprotein C-III ' , 'Apolipoprotein E' ,
                'Beta-2-microglobulin' , 'Fibrinogen alpha' , 'Fibrinogen beta' ,
                'Fibrinogen gamma' , 'Gelsolin' ,'Immunoglobulin heavy constant alpha' ,
                'Immunoglobulin heavy constant gamma' , 'Immunoglobulin heavy constant mu' ,
                'Immunoglobulin kappa constant' ,
                'Lysozyme' , 'Serum albumin' ,'Serum amyloid A' ,
                'Serum amyloid P-component' , 'Transthyretin' , 'Vitronectin')
  colonne <- colnames(Tabellona)

  #Se non ho pazienti da escludere (default=NA) o il paziente non è nella lista

  cond1 = length(id_pz_escl) == 1 && is.na(id_pz_escl)#nessuna lista passata alla funzione, caso default metto tutti i pz in tabellona
  cond2 = !(n_paziente %in% id_pz_escl) #paziente non è incluso nella lista degli esclusi

  if (cond1 || cond2){

    if(!(n_paziente %in% Tabellona$N.)){ #se il paziente non è già registrato in Tabellona

      Tabellona <- rbind(Tabellona, vector(mode = 'numeric', length = ncol(Tabellona)))
      colnames(Tabellona) <- colonne

      for(i in 1:length(idx_prot)){

        prot <- idx_prot[i]

        #trova la riga relativa alla proteina i-esima nel file targetProtein
        idx_new_dato <- grep(prot, dato_new_pz$Description)

        #trovo le colonne di Tabellona relative alla proteina i-esima
        idx_colonna <- unlist(lapply(prot, function(x) grep(x, colnames(Tabellona))))

        if (length(idx_new_dato) == 0) {#se non trova la proteina metto tutti i valori a 0

          Tabellona[nrow(Tabellona), idx_colonna] <- vector(mode = "numeric", length = length(idx_colonna))

        } else if (length(idx_new_dato) == 1){ #se trova la proteina i-esima

          Tabellona[nrow(Tabellona), idx_colonna] <- dato_new_pz[idx_new_dato, c(-1)]

        } else { #se trova più occorrenze della stessa proteina (es. + frammenti)

          #trovo la proteina con #PSMs maggiore
          max_psm_row <- which.max(dato_new_pz[idx_new_dato, "# PSMs"])

          #prendo la riga relativa alla proteina con #PSMs maggiore
          tmp <- dato_new_pz[idx_new_dato[max_psm_row], c(-1)]

          Tabellona[nrow(Tabellona), idx_colonna] <- tmp

        }
      }
      Tabellona[nrow(Tabellona), 'N.']      <- n_paziente
      Tabellona[nrow(Tabellona), 'Codice']  <- codice_pz
      Tabellona[nrow(Tabellona), 'File']    <- file_pz
      Tabellona[nrow(Tabellona), 'PSM tot'] <- PSMtot

      #-----------------------------------------------------------------------------
      #### Apolipoprotein A-IV
      #-----------------------------------------------------------------------------
      # Gestisco la Apolipoprotein A-IV, check se presente almeno una delle seguenti sequenze
      # - ARAEVSADQVATVMWDYFSQLSNNAK
      # - RAEVSADQVATVMWDYFSQLSNNAK
      # - AEVSADQVATVMWDYFSQLSNNAK

      prot <- c("Apolipoprotein A-IV")
      indices <- unlist(lapply(prot, function(x) grep(x, targetProtein$Description)))
      check_apo_A_IV <- 0

      if(length(indices) != 0){
        IDXpeptides_1 <- unlist(lapply("ARAEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
        IDXpeptides_2 <- unlist(lapply("RAEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
        IDXpeptides_3 <- unlist(lapply("AEVSADQVATVMWDYFSQLSNNAK", function(x) grep(x, targetPeptide$'Annotated Sequence')))
        cond = length(IDXpeptides_1) != 0 | length(IDXpeptides_2) != 0 | length(IDXpeptides_3) != 0
        check_apo_A_IV <- ifelse(cond, 1, 0)
      }

      # Se check==1 metto colonna APOA-IV a 1

      Tabellona[nrow(Tabellona), 'ApoA-IV'] <- ifelse(check_apo_A_IV==1, 1, 0)

      #### IMMUNOGLOBULINA Lambda-Like
      #-----------------------------------------------------------------------------
      # se presente lambda-like (idx_ll != vuoto)
      # metto colonna Lambda-Like a 1
      Tabellona[nrow(Tabellona), 'Lambda-Like'] <- ifelse(length(idx_ll)!=0, 1, 0)

      #### Daratumumab -new-
      #-----------------------------------------------------------------------------
      # se presente Daratumumab in TargetPeptideSpectrumMatch (idx_dara != vuoto)
      # metto colonna D a 1

      IDXdara <- unlist(lapply("GLEWVSAISGSGGGTYYADSVK", function(x) grep(x, targetPeptide$'Annotated Sequence')))

      Tabellona[nrow(Tabellona), 'D'] <- ifelse(length(IDXdara)!=0, 1, 0)

      #### Lambda 1 -new-
      #-----------------------------------------------------------------------------
      # se presente lambda 1 in TargetPeptideSpectrumMatch (idx_lambda != vuoto)
      # metto colonna D a 1

      IDXprot <- unlist(lapply("ANPTVTLFPPSSEELQANK", function(x) grep(x, targetPeptide$'Annotated Sequence')))

      Tabellona[nrow(Tabellona), 'Lambda 1'] <- ifelse(length(IDXprot)!=0, 1, 0)

      #### SAA -new-
      #-----------------------------------------------------------------------------
      # se presente lambda 1 in TargetPeptideSpectrumMatch (idx_lambda != vuoto)
      # metto colonna D a 1

      IDXprot <- unlist(lapply("ANPTVTLFPPSSEELQANK", function(x) grep(x, targetPeptide$'Annotated Sequence')))

      Tabellona[nrow(Tabellona), 'Lambda 1'] <- ifelse(length(IDXprot)!=0, 1, 0)





      ##### IMMUNOGLOBULINE Lambda (prendo prot a PSM maggiore)
      #-----------------------------------------------------------------------------
      # Immunoglobulin lambda constant
      # Immunoglobulin lambda-like polypeptide
      #-----------------------------------------------------------------------------

      prot_lambda <- c("Immunoglobulin lambda constant", "Immunoglobulin lambda-like")

      #trova la riga relativa alla proteina i-esima nel file targetProtein
      idx_new_dato <- unlist(lapply(prot_lambda, function(x) grep(x, dato_new_pz$Description)))

      #trovo le colonne di Tabellona relative alla proteina i-esima
      idx_colonna <- unlist(lapply(prot_lambda, function(x) grep(x, colnames(Tabellona))))

      if(length(idx_new_dato != 0)){
        # Find the indices of lambda with the maximum value of #PSMs
        max_index <- which.max(dato_new_pz[idx_new_dato, "# PSMs"])

        # Select only the index with the maximum value of #PSMs
        ind_lambda_max <- idx_new_dato[max_index]

        Tabellona[nrow(Tabellona), idx_colonna] <- dato_new_pz[ind_lambda_max, c(-1)]

      } else{ #se non trova nessuna delle lambda mette tutto a zero
        Tabellona[nrow(Tabellona), idx_colonna] <- vector(mode = "numeric", length = length(idx_colonna))
      }


      #-----------------------------------------------------------------------------
      #-----------------------------------------------------------------------------

      #Metto in ordine per N. paziente
      Tabellona <- Tabellona[order(Tabellona$'N.'), ]

      #Sovrascrivo il file
      write.xlsx(Tabellona, Tabellona_PATH)

      #Creo tutta la formattazione del workbook
      file <- loadWorkbook(Tabellona_PATH)

      #Metto i bordi alle celle
      fill_style <- createStyle(fontColour = "red", textDecoration = 'bold', border = 'bottom')
      addStyle(file, sheet = 1, rows = 1, cols = 1:ncol(Tabellona), style = fill_style, gridExpand = TRUE)

      index <- which(colnames(Tabellona) == 'Peptidi')
      fill_style <- createStyle(fgFill = "gray90")
      addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = 1:index, style = fill_style, gridExpand = TRUE, stack = TRUE)
      fill_style <- createStyle(border = 'right')
      addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = index, style = fill_style, gridExpand = TRUE, stack = TRUE)

      #metto insieme anche la immunoglobulina lambda
      idx_prot_colors =   idx_prot <- c('Apolipoprotein A-I ' , 'Apolipoprotein A-II ' , 'Apolipoprotein A-IV ' ,
                                        'Apolipoprotein C-II ' , 'Apolipoprotein C-III ' , 'Apolipoprotein E' ,
                                        'Beta-2-microglobulin' , 'Fibrinogen alpha' , 'Fibrinogen beta' ,
                                        'Fibrinogen gamma' , 'Gelsolin' ,'Immunoglobulin heavy constant alpha' ,
                                        'Immunoglobulin heavy constant gamma' , 'Immunoglobulin heavy constant mu' ,
                                        'Immunoglobulin kappa constant' , "Immunoglobulin lambda constant",
                                        'Lysozyme' , 'Serum albumin' ,'Serum amyloid A' ,
                                        'Serum amyloid P-component' , 'Transthyretin' , 'Vitronectin')

      for (i in 1:length(idx_prot_colors)){
        if(i%%2 != 0){
          index <- unlist(lapply(idx_prot_colors[i], function(x) grep(x, colnames(Tabellona))))

          fill_style <- createStyle(border = 'right')
          addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = index[length(index)], style = fill_style, gridExpand = TRUE, stack = TRUE)

          fill_style <- createStyle(fgFill = "peachpuff")
          addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = index[1]:index[length(index)], style = fill_style, gridExpand = TRUE, stack = TRUE)

        } else {
          index <- unlist(lapply(idx_prot_colors[i], function(x) grep(x, colnames(Tabellona))))

          fill_style <- createStyle(border = 'right')
          addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = index[length(index)], style = fill_style, gridExpand = TRUE, stack = TRUE)

          fill_style <- createStyle(fgFill = "thistle2")
          addStyle(file, sheet = 1, rows = 1:(nrow(Tabellona)+1), cols = index[1]:index[length(index)], style = fill_style, gridExpand = TRUE, stack = TRUE)

        }

      }
      setColWidths(file, 1, widths = 'auto', cols = 1:200)
      saveWorkbook(file, Tabellona_PATH, overwrite = TRUE)

    } else {
      cat('\n\n !!! ATTENZIONE: PAZIENTE', n_paziente, 'GIÀ ESISTENTE IN TABELLONA !!!')
    }
  } else { #è stata passata una lista di pazienti da escludere dalla Tabellona in cui n_paziente è presente
    cat('\n\n !!! ATTENZIONE: IL PAZIENTE NON PUÒ ESSERE INSERITO IN TABELLONA !!!')
    cat("\n  Paziente N°", n_paziente, "risulta incluso in quelli da non inserire in Tabellona!")
  }

}
