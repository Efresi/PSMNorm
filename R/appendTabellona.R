appendTabellona <- function(targetPeptide_name, targetProtein_name, Tabellona_name) {

  #funzione che aggiunge i dati relativi alle proteine di un paziente alla tabellona
  library(stringi)

  #Carico file TargetPeptideSpectrumMatch
  targetPeptide <- read.csv(targetPeptide_name , sep = "" )

  #Carico TargetProtein
  targetProtein <- read.csv(targetProtein_name, sep = "" )

  colnames(targetProtein) <- c("Checked", "Master", "Accession", "Description", "Contaminant", "D",
                               "Coverage", "# PSMs", "Score Sequest HT", "# Peptides", "# Unique Peptides",
                               "# Protein Groups", "# AAs", "MW [kDa]", "calc. pI", "Modifications",
                               "# Peptides Sequest HT")

  #prendo numero paziente, codice e file dal nome file
  #prendo PSM tot dal numero righe di TargetPeptide

  tmp <- stri_split_boundaries(targetPeptide_name)[[1]] #attenzione

  n_paziente <- as.numeric(stri_replace_all_fixed(tmp[1], " ", "")) #attenzione
  codice_pz  <- stri_replace_all_fixed(tmp[3], " ", "")#attenzione
  file_pz    <- stri_split_fixed(tmp[length(tmp)], "_")[[1]]
  file_pz    <- stri_c(file_pz[2], "_", file_pz[3])
  PSMtot <- nrow(targetPeptide)

  #prendo coverage, #psm e score sequest ht dal file TargetProtein
  amyPROT <- c("Apolipoprotein A-I", "Apolipoprotein C", "Apolipoprotein E", "Beta-2-microglobulin",
               "Fibrinogen", "Gelsolin",
               "Immunoglobulin heavy constant alpha", "Immunoglobulin heavy constant gamma",
               "Immunoglobulin heavy constant mu", "Immunoglobulin kappa constant",
               "Immunoglobulin lambda constant",
               "Lysozyme", "Serum albumin", "Serum amyloid",
               "Transthyretin", "Vitronectin")

  columns <- c("Description","# AAs", "Coverage", "# PSMs", "Score Sequest HT")

  idx <- unlist(lapply(amyPROT, function(x) grep(x, targetProtein$Description)))
  dato_new_pz <- targetProtein[idx, columns]

  library(openxlsx)

  Tabellona <- openxlsx::read.xlsx(Tabellona_name, sheet = 1, sep.names = " ")

  idx_prot <- c('Apolipoprotein A-I ' , 'Apolipoprotein A-II ' , 'Apolipoprotein A-IV ' ,
                'Apolipoprotein C-II ' , 'Apolipoprotein C-III ' , 'Apolipoprotein E' ,
                'Beta-2-microglobulin' , 'Fibrinogen alpha' , 'Fibrinogen beta' ,
                'Fibrinogen gamma' , 'Gelsolin' ,'Immunoglobulin heavy constant alpha' ,
                'Immunoglobulin heavy constant gamma' , 'Immunoglobulin heavy constant mu' ,
                'Immunoglobulin kappa constant' ,  'Immunoglobulin lambda constant' ,
                'Lysozyme' , 'Serum albumin' ,'Serum amyloid A' ,
                'Serum amyloid P-component' , 'Transthyretin' , 'Vitronectin')

  if(!(n_paziente %in% Tabellona$N.)){ #Se il paziente non è già registrato in Tabellona

    Tabellona <- rbind(Tabellona, vector(mode = 'numeric', length = ncol(Tabellona)))

    for(i in 1:length(idx_prot)){

      prot <- idx_prot[i]

      #trova la riga relativa alla proteina i-esima nel file targetProtein
      idx_new_dato <- grep(prot, dato_new_pz$Description)

      #trovo le colonne di Tabellona relative alla proteina i-esima
      idx_colonna <- unlist(lapply(prot, function(x) grep(x, colnames(Tabellona))))

      if (length(idx_new_dato) == 0) {#se non trova la proteina metto tutti i valori a 0

        Tabellona[nrow(Tabellona), idx_colonna] <- vector(mode = "numeric", length = 4)

      } else if (length(idx_new_dato) == 1){ #se trova la proteina i-esima

        Tabellona[nrow(Tabellona), idx_colonna] <- dato_new_pz[idx_new_dato, columns[-1]]

      } else { #se trova più occorrenze della stessa proteina (es. + frammenti)
        #prendo il valore con #PSMs maggiore
        max_psm_row <- which.max(dato_new_pz[idx_new_dato, "# PSMs"])

        Tabellona[nrow(Tabellona), idx_colonna] <- dato_new_pz[idx_new_dato[max_psm_row], columns[-1]]

      }
    }
    Tabellona[nrow(Tabellona), 'N.'] <- n_paziente
    Tabellona[nrow(Tabellona), 'Codice'] <- codice_pz
    Tabellona[nrow(Tabellona), 'File'] <- file_pz
    Tabellona[nrow(Tabellona), 'PSM tot'] <- PSMtot

    #Metto in ordine per N. paziente
    Tabellona <- Tabellona[order(Tabellona$'N.'), ]
    #Sovrascrivo il file

    #creo tutta la formattazione del workbook


    openxlsx::write.xlsx(Tabellona, Tabellona_name)
  } else {
    cat('\n\n !!! ATTENZIONE: PAZIENTE GIÀ ESISTENTE IN TABELLONA !!!')
  }
}

