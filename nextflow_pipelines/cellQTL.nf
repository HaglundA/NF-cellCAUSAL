
nextflow.enable.dsl=2

params.ncores=20
params.mem=300


params.outdir="/rds/general/user/ah3918/projects/roche/live/ALEX/ANALYSIS/PIPELINE_RUN/nf_run/"

params.pseudobulk=false
params.genotype=false
params.gene_locs=false
params.matrixeqtl=false
params.summarize_outputs=false
params.count_data=false
params.merge=false
params.transeqtl=false

params.cov_file="/rds/general/user/ah3918/projects/roche/live/ALEX/METADATA/MASTER_COVARIATE_MATRIX.txt"
params.run_name="testing.nf"
params.full_sample_list="/rds/general/project/roche/live/ALEX/METADATA/full_samples.txt"
params.disease_sample_list="/rds/general/project/roche/live/ALEX/METADATA/disease_samples.txt"
params.control_sample_list="/rds/general/project/roche/live/ALEX/METADATA/control_samples.txt"
params.seurat_assay="decontXcounts"


process create_genotype{

    module "anaconda3/personal"
    conda "/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS"
    executor="pbspro"
    clusterOptions = "-lselect=1:ncpus=${params.ncores}:mem=${params.mem}gb -l walltime=1:00:00"
    publishDir "${params.outdir}/MatrixEQTL_IO", mode: "copy"

    input:
    path genofile

    output:
    path "genotype_012mat.csv", emit: genomat
    path "snp_chromlocations.csv", emit: snplocs 
    path "MAF_mat.csv", emit: mafmat


    """
    #!/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS/bin/Rscript

    library(SCGsuite)

    get_genotype_matrix(inputfile="$genofile")



    """


}

process pseudobulk_singlecell{

    module "anaconda3/personal"
    conda "/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS"
    executor="pbspro"
    clusterOptions = "-lselect=1:ncpus=40:mem=480gb -l walltime=5:00:00"
    publishDir "/rds/general/user/ah3918/projects/roche/live/ALEX/ANALYSIS/PIPELINE_RUN/nf_run/MatrixEQTL_IO", mode: "copy"


    output:
    path "*"



    """
    #!/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS/bin/Rscript

    library(SCGsuite)
                    
    seuratpath="/rds/general/user/ah3918/projects/roche/live/ALEX/ANALYSIS/SEURAT_OBJECTS/"

    input_seurat_list_files=c(paste0(seuratpath,"BRYOIS_92_AD_final_annotated_decontXfiltered.rds"),
    paste0(seuratpath,"BRYOIS_92_MS_final_annotated_decontXfiltered.rds"),
    paste0(seuratpath,"MRC_60_final_annotated_decontXfiltered.rds"),
    paste0(seuratpath,"ROCHE_MPD92_final_annotated_decontXfiltered.rds"),
    paste0(seuratpath,"MATTHEWS_final_annotated_decontXfiltered.rds"))



    input_seurat_list = TRUE
    celltype_colname="CellType"
    indiv_colname="Individual_ID"
    min.cells=10

    message("Using ${params.seurat_assay} assay to pseudobulk")


    message("First, making sure that all Seurat objects have the same genes..")
      so=readRDS(input_seurat_list_files[[1]])
      Seurat::DefaultAssay(so)="${params.seurat_assay}"
      commongenes=rownames(so)
      for(i in 2:length(input_seurat_list_files)){
        so=readRDS(input_seurat_list_files[[i]])
        Seurat::DefaultAssay(so)="${params.seurat_assay}"
        genes=rownames(so)
        commongenes=intersect(genes,commongenes)
      }

      
      n_seurat_objs=length(input_seurat_list_files)
      message(paste0(n_seurat_objs," seurat objects were provided. Reading in Seurat obj 1 .."))
      seuratobj<-readRDS(input_seurat_list_files[[1]])
      seuratobj=seuratobj[commongenes,]
      
      celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
      agg_count_list_full<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay="${params.seurat_assay}")

      for(i in 2:length(input_seurat_list_files)){
        message(paste0(" Reading in Seurat obj ",i,".."))
        seuratobj<-readRDS(input_seurat_list_files[[i]])
        seuratobj=seuratobj[commongenes,]
        celltypelist<-Seurat::SplitObject(seuratobj,split.by=celltype_colname)
        agg_count_list<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay="${params.seurat_assay}")
        agg_count_list<-agg_count_list[names(agg_count_list) %in% names(agg_count_list_full)]
        agg_count_list<-agg_count_list[names(agg_count_list_full)]
        agg_count_list_full<-Map(cbind,agg_count_list_full,agg_count_list)
      }
      nullvec=unlist(lapply(agg_count_list_full,is.null))
      dropped_celltypes=names(nullvec[which(nullvec==TRUE)])
      agg_count_list_full=Filter(Negate(is.null),agg_count_list_full)


      agg_count_list_full=normalize_pseudobulk(agg_count_list_full)


      if(length(dropped_celltypes)>0){
        message(paste0(dropped_celltypes," celltype was dropped due to no individuals passing min.cells criteria."))
      }


      for(i in 1:length(agg_count_list_full)){
          write.table(agg_count_list_full[[i]],paste0(names(agg_count_list[i]),"_pseudobulk.csv"))
        }
    

    ### 16th Jan 2024 - Add in a "fake" pseudobulk, sum counts for ALL cells

    print("Adding in a fake pseudobulk, summing counts for all cells")
    seuratobj<-readRDS(input_seurat_list_files[[1]])
    celltypelist<-list(seuratobj)
    names(celltypelist)<-"Bulk"
      agg_count_list_full<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay="${params.seurat_assay}")

      for(i in 2:length(input_seurat_list_files)){
        message(paste0(" Reading in Seurat obj ",i,".."))
        seuratobj<-readRDS(input_seurat_list_files[[i]])
        seuratobj=seuratobj[commongenes,]
        celltypelist<-list(seuratobj)
        names(celltypelist)<-"Bulk"
        agg_count_list<-pseudobulk_counts(celltypelist,min.cells=min.cells,indiv_col=indiv_colname,assay="${params.seurat_assay}")
        agg_count_list<-agg_count_list[names(agg_count_list) %in% names(agg_count_list_full)]
        agg_count_list<-agg_count_list[names(agg_count_list_full)]
        agg_count_list_full<-Map(cbind,agg_count_list_full,agg_count_list)
      }
      nullvec=unlist(lapply(agg_count_list_full,is.null))
      dropped_celltypes=names(nullvec[which(nullvec==TRUE)])
      agg_count_list_full=Filter(Negate(is.null),agg_count_list_full)


      agg_count_list_full=normalize_pseudobulk(agg_count_list_full)

      for(i in 1:length(agg_count_list_full)){
          write.table(agg_count_list_full[[i]],paste0(names(agg_count_list[i]),"_pseudobulk.csv"))
        }




    """


}

