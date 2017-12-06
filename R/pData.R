#' Constructor for pData
#'
#' This function allows constructor of phenotype information necessary for downstream analysis. See format of required fields. Function will also rearrange the rows such that trios are grouped together - with proband first, mother second, and father third. 
#' @param path The path to tab-delimited file storing the phenotype information.
#' @export
pData <- function(path){
	tab = read.table(path, header=T, colClasses = "character")
	col_ind = match(c("subj_id", "family_id", "father_id", "mother_id", "gender", "bam_path"), colnames(tab))
	tab_selected = tab[,col_ind]
	proband_id = sort(tab_selected$subj_id[tab_selected$father_id %in% tab_selected$subj_id & 
		tab_selected$mother_id %in% tab_selected$mother_id])
	
	tab_proband = tab_selected[match(proband_id, tab_selected$subj_id),]
	tab_mother = tab_selected[match(tab_proband$mother_id, tab_selected$subj_id),]
	tab_father = tab_selected[match(tab_proband$father_id, tab_selected$subj_id),]

	tab_out = tab[1:length(proband_id)*3, col_ind]
	tab_out[seq(1, length(proband_id)*3, by=3),] = tab_proband
	tab_out[seq(2, length(proband_id)*3, by=3),] = tab_mother
	tab_out[seq(3, length(proband_id)*3, by=3),] = tab_father
	return(tab_out)
}