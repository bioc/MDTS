#' Constructor for pData
#'
#' This function allows constructor of phenotype information necessary for 
#' downstream analysis. See format of required fields. Function will also 
#' rearrange the rows such that trios are grouped together - with proband first,
#' mother second, and father third. 
#' @param path The path to tab-delimited file storing the phenotype information.
#' @examples 
#'	system.file('extdata', package='MDTS')
#'	pD = pData('https://raw.githubusercontent.com/JMF47/MDTSData/master/data/pD.ped')
#' @export
#' @return Returns a \code{data.frame} of required sample information for 
#' running MDTS.
pData <- function(path){
	tab = utils::read.table(path, header=TRUE, colClasses = "character")
	hdrs = c("subj_id", "family_id", "father_id", 
	         "mother_id", "gender", "bam_path")
	col_ind = match(hdrs, colnames(tab))
	tab_sel = tab[,col_ind]
	proband_id = sort(tab_sel$subj_id[
	      tab_sel$father_id %in% tab_sel$subj_id & 
		tab_sel$mother_id %in% tab_sel$mother_id])
	
	tab_proband = tab_sel[match(proband_id, tab_sel$subj_id),]
	tab_mother = tab_sel[match(tab_proband$mother_id, tab_sel$subj_id),]
	tab_father = tab_sel[match(tab_proband$father_id, tab_sel$subj_id),]

	tab_out = tab[1:length(proband_id)*3, col_ind]
	tab_out[seq(1, length(proband_id)*3, by=3),] = tab_proband
	tab_out[seq(2, length(proband_id)*3, by=3),] = tab_mother
	tab_out[seq(3, length(proband_id)*3, by=3),] = tab_father
	return(tab_out)
}