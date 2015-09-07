#' \code{PQI_to_MA} converts protein QI data output
#' @rdname QI_to_MA
#' @export
PQI_to_MA <- function(directory = getwd(),
                      abundance = "Normalised",
                      facNames  = NULL,
                      compoundID = "shortDescription") {
    oldwd <- getwd()
    on.exit(setwd(oldwd))    
    setwd(directory)
    match.arg(abundance, c("Normalised", "Raw"))
    match.arg(compoundID, c("Accession", "Description", "shortDescription"))
    
    rawfiles <- list.files(pattern = ".csv", full.names = TRUE, include.dirs = TRUE)
    
    for (j in 1:length(rawfiles)) {
        in.tab <- read.csv(rawfiles[j], stringsAsFactors = FALSE, header = FALSE)
        tab <- in.tab
        
        if (compoundID == "Accession") {
            idcmpd <- pmatch("Acces", tab[3, ])
        } else {
            idcmpd <- pmatch("Descr", tab[3, ])
        }
        
        noempty2 <- which(tab[2, ] != "")
        
        # Filling the sample names
        for (i in 1:(length(noempty2) - 1)) {
            tab[2, noempty2[i]:(noempty2[i + 1] - 1)] <- tab[2, noempty2[i]]
        }
        
        # Cutting the table, moving protein ID column into the beginning
        normstart <- pmatch("Norm", in.tab[1, ])
        rawstart  <- pmatch("Raw", in.tab[1, ])
        specstart <- pmatch("Spec", in.tab[1, ])
        if (abundance == "Normalised") {
            tab <- tab[c(idcmpd, normstart:(rawstart - 1))]
        } else {
            tab <- tab[c(idcmpd, rawstart:(specstart - 1))]
        }
        
        if (compoundID == "shortDescription") {
            for (i in 4:nrow(tab)) tab[i,1] <- strsplit(tab[i,1], " OS=")[[1]][1]
        }
        
        tab[1,]  <- tab[3,]
        tab[1,1] <- "Sample"
        
        # replace spaces with "_" in compound names
        tab[1] <- sapply(tab[1], gsub, pattern = " ", replacement = "_")
        
        sp <- strsplit(as.character(tab[2,-1]), "_")
        facNum <- unique(sapply(sp, length))
        if (length(facNum) != 1) {
            warning("Different number of factors in group labels!
                    Empty factor values are filled with NAs");
            facNum <- max(facNum);
        }
        
        if (is.null(facNames)) {
            if (facNum == 1) {
                facNames <- "Group"
            } else {
                facNames <- paste0("Factor", 1:facNum);
            }
        } else {
            if (length(facNames) < facNum) {
                warning("The number of factor names is higher than the number of detected factors.
                        Extra names are removed");
                facNames <- facNames[1:facNum];
            } else if (length(facNames) > facNum){
                warning("The number of factor names is lower than the number of detected factors.
                        Extra factor names are generated automatically");
                facNames <- c(facNames, paste0("Factor", length(facNames):facNum));
            }
        }
        
        groupData <- data.frame(cbind(facNames, sapply(sp, '[', 1:max(facNum))))
        names(groupData) <- names(tab)
        tab <- rbind(tab[1, ], groupData, tab[-(1:3), ])
        
        # Write the output file
        write.table (tab,
                     file = paste0(strsplit(rawfiles[1], ".csv")[[1]], "_for_MA.csv"),
                     col.names = FALSE, row.names = FALSE, sep = ",")
        }
    return("DONE")
    }