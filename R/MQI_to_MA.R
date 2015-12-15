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
#' @param path Character vector. There could be two types of elements:\cr
#' 1) file path(s) to the proceeded .csv QI output file(s).\cr
#' 2) path(s) to the directory(ies) containing .csv files.\cr
#' In the second case all the files in all subdirectories are proceeded.\cr
#' Remember that in file paths you should use either "/" or double "\", not just "\".
#' @param abundance If \code{"Normalised"}, then normalised data are used; if \code{"Raw"}, then raw data are used
#' @param autoFac If \code{"TRUE"} then group names are automatically converted to grouping factors' values.
#' @param facNames Character vector of grouping factor names.
#' @param compoundID The name of the column in QI output table used as the set of compound IDs.\cr
#' For \code{MQI_to_MA} one of: \code{"Compound"}, \code{"Accepted Compound ID"}, \code{"Formula"}.\cr
#' For \code{PQI_to_MA} one of: \code{"Accession"}, \code{"Description"}, \code{"shortDescription"}.\cr
#' \code{"shortDescription"} means that the data from "Description" column are used, 
#' but they are cut starting with "OS=".
#'
#' @examples
#' # PQI_to_MA("//mpimp-golm/user/Homes/Zubkov/QI_data")
#' # PQI_to_MA("\\\\mpimp-golm\\user\\Homes\\Zubkov\\QI_data")
#' # PQI_to_MA("C:/QI_data", abundance = "Raw", compoundID = "Accession")
#' # MQI_to_MA("C:/QI_data", facNames = c("Phenotype", "Treatment", "Time"))
#' @name QI_to_MA
NULL

#' \code{MQI_to_MA} converts lipid and metabolite QI data output
#' @rdname QI_to_MA
#' @export

MQI_to_MA <- function(file = getwd(),
                      abundance = "Raw",
					  compoundID = "Accepted Compound ID",
                      autoFac = FALSE,
                      facNames = NULL) {

    match.arg(abundance, c("Normalised", "Raw"))
    match.arg(compoundID, c("Compound", "Formula", "Accepted Compound ID"))
	
	files <- grepl(".csv", path)
	dirs  <- !files
	rawfiles <- c(path[files], unlist(sapply(path[dirs], list.files, pattern = ".csv", full.names = TRUE, include.dirs = TRUE)))
	#facNums <- integer()
    MAfiles <- character()
	
    for (j in 1:length(rawfiles)) {
        format.check <- readLines(rawfiles[j], n = 10)
        
        if (all(grepl(";", format.check))) {
            sep <- ";"
        } else if (all(grepl(",", format.check))) {
            sep <- ","
        } else {
            sep <- ""
        }
        
        in.tab <- read.table(rawfiles[j],
                             sep = sep, dec = dec,
                             stringsAsFactors = FALSE,
                             header = FALSE,
                             colClasses = "character",
                             comment.char = "")
        tab <- in.tab
        
        idcmpd <- pmatch(compoundID, tab[3, ])
        final <- min(pmatch("Accepted Compound ID", tab[3, ]), ncol(tab) + 1, na.rm = T)
        noempty2 <- which(tab[2, ] != "")
        noempty2 <- c(noempty2, ncol(tab) + 1)
        
        # Filling the sample names
        for (i in 1:(length(noempty2) - 1)) {
            tab[2, noempty2[i]:(noempty2[i + 1] - 1)] <- tab[2, noempty2[i]]
        }
        
        # Cutting the table, removing CompoundID column into the beginning
        normstart <- pmatch("Norm", in.tab[1, ])
        rawstart  <- pmatch("Raw", in.tab[1, ])
        if (abundance == "Normalised") {
            tab <- tab[c(idcmpd, normstart:(rawstart - 1))]
        } else {
            tab <- tab[c(idcmpd, rawstart:(final - 1))]
        }
        
        tab[1, ]  <- tab[3, ]
        tab[1, 1] <- "Sample"
		# remove commas and replace spaces with "_" in compound and sample names
        tab[[1]] <- gsub(tab[[1]], pattern = " ", replacement = "_")
		tab[[1]] <- gsub(tab[[1]], pattern = ",", replacement = "_")
		tab[1, ] <- gsub(tab[1, ], pattern = " ", replacement = "_")
		tab[1, ] <- gsub(tab[1, ], pattern = ",", replacement = "_")
        
        # replace commas in intensity values with dots
		# tab <- apply(tab, 2, gsub, pattern = ",", replacement = ".")
        
        if (!autoFac) {
            tab <- tab[-(2:3), ]
            
        } else {
            sp <- strsplit(as.character(tab[2,-1]), "_")
            facNum <- unique(sapply(sp, length))
            if (length(facNum) != 1) {
                warning("Different number of factors in group labels! Empty factor values are filled with NAs");
                facNum <- max(facNum);
            }
            
            if (is.null(facNames)) {
                if (facNum == 1) {
                    facNames <- "Group"
                } else {
                    facNames <- paste0("Factor", 1:facNum);
                }
            } else {
                if (length(facNames) > facNum) {
                    warning("The number of factor names is higher than the number of detected factors.
                        Extra names are removed");
                    facNames <- facNames[1:facNum];
                } else if (length(facNames) < facNum){
                    warning("The number of factor names is lower than the number of detected factors.
                        Extra factor names are generated automatically");
                    facNames <- c(facNames, paste0("Factor", length(facNames):facNum));
                }
            }
            
            groupData <- data.frame(cbind(facNames, sapply(sp, '[', 1:max(facNum))), stringsAsFactors = FALSE)
            names(groupData) <- names(tab)
            tab <- rbind(tab[1, ], groupData, tab[-(1:3), ])
        }
		
        
		# Write the output file
		MAfiles[j] <- paste0(strsplit(rawfiles[j], ".csv")[[1]], "_for_MA.csv")
		#facNums[j] <- facNum
        write.table (tab,
                     file = MAfiles[j],
                     col.names = FALSE, row.names = FALSE, sep = ";", dec = ",")
        }
    
#     # Combine all the pairs of neg and pos files
#     if (unite_neg_pos) {
#         posfiles <- MAfiles[grepl("pos", MAfiles)]
# 		negfiles <- MAfiles[grepl("neg", MAfiles)]
#         if (length(posfiles) == 0 || length(negfiles) == 0) {
#             return("DONE")
#         }
# 		
# 		posnames <- lapply(strsplit(posfiles, split="pos"), '[', 1)
# 		negnames <- lapply(strsplit(negfiles , split="neg"), '[', 1)		
#         paired <- match(negnames, posnames, nomatch = 0)
#         posfiles <- posfiles[paired > 0]
#         negfiles <- negfiles[paired]
# 		facNums <- facNums[grepl("pos", MAfiles)][paired > 0]
# 		
#         if (length(posfiles) == 0 || length(negfiles) == 0) {
#             return("DONE")
#         }
# 		
#         for (j in 1:length(posfiles)) {
# 			facNum <- facNums[j]
#             tab.pos <- read.csv(posfiles[j], stringsAsFactors = FALSE, header = FALSE)
#             tab.neg <- read.csv(negfiles[j], stringsAsFactors = FALSE, header = FALSE)
# 			tab.pos[1, grep("pos", tab.pos)] <- sapply(tab.pos[1, grep("pos", tab.pos)], strsplit, "_pos")
# 			tab.neg[1, grep("neg", tab.neg)] <- sapply(tab.neg[1, grep("neg", tab.neg)], strsplit, "_neg")
#             if (!identical(tab.pos[1:(facNum + 1), ], tab.neg[1:(facNum + 1), ])) {
# 				warning("Sample names or grouping data does not match in paired files ", 
# 						posfiles[j], " and ", negfiles[j], "\nThese files will not be merged.")
# 			} else {
# 				tab <- rbind(tab.pos, tab.neg[-(1:(facNum + 1)), ])
# 			}
# 			
#             write.table (tab,
#                          file = paste0(strsplit(posfiles[j], "pos")[[1]][1], "_for_MA.csv"),
#                          col.names = FALSE, row.names = FALSE, sep = ",")
#         }
#     }
    
    return("DONE")
    }
