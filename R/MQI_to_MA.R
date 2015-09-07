#' Convert QI to MetaboAnalyst
#'
#' Function for converting QI output .csv files into format convenient for further analysis with MetaboAnalyst.
#' For proper work, please, make complete QI output files with all the options marked.\cr
#' This function:
#' \enumerate{
#'   \item Cuts out all the extra data, leaving only columns with abundance values 
#' (either raw or normalised) and a column with compound's IDs.
#'   \item Moves chosen compound ID column to the first position.
#'   \item Expands sample group labels.
#'   \item Separate group labels into grouping factors and creates rows with factors.
#'   \item (for \code{MQI_to_MA}) Merges neg and pos files with matching names.
#'   \item Writes the resulting table into new .csv file.
#' }
#' 
#' @param directory  Path to the directory containing .csv QI output files.
#' All the files in all subdirectories are proceeded too.
#' Remember that in file paths you should use either "/" or double "\", not just "\".
#' @param abundance If \code{"Normalised"}, then normalised data are used; if \code{"Raw"}, then raw data are used
#' @param facNames Character vector of grouping factor names.
#' @param compoundID The name of the column in QI output table used as the set of compound IDs.\cr
#' For \code{MQI_to_MA} one of: \code{"Compound"}, \code{"Accepted Compound ID"}, \code{"Formula"}.\cr
#' For \code{PQI_to_MA} one of: \code{"Accession"}, \code{"Description"}, \code{"shortDescription"}.\cr
#' \code{"shortDescription"} means that the data from "Description" column are used, 
#' but they are cut starting with "OS=".
#'
#' @examples
#' PQI_to_MA("//mpimp-golm/user/Homes/Zubkov/QI_data")
#' PQI_to_MA("\\\\mpimp-golm\\user\\Homes\\Zubkov\\QI_data")
#' PQI_to_MA("C:/QI_data", abundance = "Raw", compoundID = "Accession")
#' MQI_to_MA("C:/QI_data", facNames = c("Phenotype", "Treatment", "Time"), unite_neg_pos = FALSE)
#' @name QI_to_MA
NULL

#' \code{MQI_to_MA} converts lipid and metabolite QI data output
#' @param unite_neg_pos If \code{TRUE}, creates combined file from pos and neg file pairs with matching names;
#' if \code{FALSE}, just skips this step
#' @rdname QI_to_MA
#' @export

MQI_to_MA <- function(directory = getwd(),
                      abundance = "Normalised",
					  compoundID = "Accepted Compound",
                      facNames = NULL,
                      unite_neg_pos = TRUE) {
    oldwd <- getwd()
    on.exit(setwd(oldwd))    
    setwd(directory)
    match.arg(abundance, c("Normalised", "Raw"))
    
    rawfiles <- list.files(pattern = "raw.csv", full.names = TRUE, include.dirs = TRUE)
    
    for (j in 1:length(rawfiles)) {
        in.tab <- read.csv(rawfiles[j], stringsAsFactors = FALSE, header = FALSE)
        tab <- in.tab
        
        idcmpd <- pmatch("Accepted Compound", tab[3, ])
        noempty2 <- which(tab[2, ] != "")
        noempty2 <- c(noempty2, idcmpd)
        
        # Filling the sample names
        for (i in 1:(length(noempty2) - 1)) {
            tab[2, noempty2[i]:(noempty2[i + 1] - 1)] <- tab[2, noempty2[i]]
        }
        
        # Cutting the table, removing "Accepted Compound ID" into the beginning
        normstart <- pmatch("Norm", in.tab[1, ])
        rawstart  <- pmatch("Raw", in.tab[1, ])
        if (abundance == "Normalised") {
            tab <- tab[c(idcmpd, normstart:(rawstart - 1))]
        } else {
            tab <- tab[c(idcmpd, rawstart:(idcmpd - 1))]
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
    
    # Combine all the pairs of neg and pos files
    if (unite_neg_pos) {
        MAfiles  <- list.files(pattern = "for_MA.csv")
        posfiles <- MAfiles[grep("_pos_", MAfiles)]
        if (length(posfiles) == 0 || length(negfiles) == 0) {
            setwd(oldwd)
            return("DONE")
        }
        negfiles <- MAfiles[grep(paste0(unlist(strsplit(posfiles, "_pos_raw_for_MA.csv")), "_neg"), MAfiles)]
        paired   <- sapply(negfiles, length) == 1
        posfiles <- posfiles[paired]
        negfiles <- negfiles[paired]
        if (length(posfiles) == 0 || length(negfiles) == 0) {
            setwd(oldwd)
            return("DONE")
        }
        for (j in 1:length(posfiles)) {
            tab.pos <- read.csv(posfiles[j], stringsAsFactors = FALSE, header = FALSE)
            tab.neg <- read.csv(negfiles[j], stringsAsFactors = FALSE, header = FALSE)
            if (emptyrow) {
                tab <- rbind(tab.pos, tab.neg[-c(1:3), ])
                tab[3, grep("neg", tab)] <- sapply(tab[3, grep("neg", tab)], strsplit, "_neg")
            } else {
                tab <- rbind(tab.pos, tab.neg[-c(1:2), ])
                tab[2, grep("neg", tab)] <- sapply(tab[2, grep("neg", tab)], strsplit, "_neg")
            }
            write.table (tab, 
                         file = paste0(strsplit(posfiles[1], "_pos_raw_for_MA.csv")[[1]], "_for_MA.csv"),
                         col.names = FALSE, row.names = FALSE, sep = ",")
        }
    }
    
    return("DONE")
    }