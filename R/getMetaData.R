#' Constructor for metadata
#'
#' This function allows constructor of phenotype information necessary for 
#' downstream analysis. See format of required fields. Function will also 
#' rearrange the rows such that trios are grouped together - with proband first,
#' mother second, and father third. 
#' @param path The path pointing to the file that contains information on each
#' subject in the dataset.
#' @param id The column name that identifies each sample. Defaults to 'subj_id'.
#' @param familyId The column name that identifies which family the sample
#' belongs to. Defaults to 'family_id'.
#' @param fatherId The column name that identifies the id of the father. 
#' Defaults to 'father_id'.
#' @param motherId The column name that identifies the id of the mother.
#' Defaults to 'mother_id'.
#' @param bamPath The column name that identifies where to find the bam file
#' for each subject. Defaults to 'bam_path'.
#' @examples 
#'	meta <- getMetaData(
#'	'https://raw.githubusercontent.com/JMF47/MDTSData/master/data/pD.ped')
#' @export
#' @return Returns a \code{data.frame} of required sample information for 
#' running MDTS.
getMetaData <- function(path, id="subj_id", familyId="family_id", 
                        fatherId="father_id", motherId="mother_id", 
                        bamPath="bam_path"){
	tab <- utils::read.table(path, header=TRUE, colClasses = "character")
	hdrs <- c(id, familyId, fatherId, motherId, bamPath)
	if(sum(hdrs %in% colnames(tab))<length(hdrs))
	      stop("One of the required variables is missing.")
	col_ind <- match(hdrs, colnames(tab))
	tab_sel <- tab[,col_ind]
	colnames(tab_sel) <- c("subj_id", "family_id", "father_id", "mother_id",
	                       "bam_path")
	proband_id <- sort(tab_sel$subj_id[
	      tab_sel[[fatherId]] %in% tab_sel[[id]] & 
		tab_sel[[motherId]] %in% tab_sel[[id]]])
	
	tab_proband <- tab_sel[match(proband_id, tab_sel[[id]]),]
	tab_mother <- tab_sel[match(tab_proband$mother_id, tab_sel[[id]]),]
	tab_father <- tab_sel[match(tab_proband$father_id, tab_sel[[id]]),]

	tab_out <- tab[,col_ind]
	tab_out[seq(1, length(proband_id)*3, by=3),] <- tab_proband
	tab_out[seq(2, length(proband_id)*3, by=3),] <- tab_mother
	tab_out[seq(3, length(proband_id)*3, by=3),] <- tab_father
	return(tab_out)
}