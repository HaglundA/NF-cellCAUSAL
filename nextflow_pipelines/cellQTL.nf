
nextflow.enable.dsl=2


//general parameters
params.conda_env_path = "/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS"
params.outdir="/rds/general/user/ah3918/ephemeral/cellCAUSAL_TEST/"

// Pseudobulking parameters
params.min_cells = 10
params.indiv_column = "Individual_ID"
params.celltype_column = "CellType"
params.assay = "decontXcounts"

// Genotype parameters
params.matrix = "path/to/your/matrix"
params.vcf = "/rds/general/user/ah3918/projects/roche/live/ALEX/PROCESSED_DATA/PROCESSED_GENOTYPE/FINAL/final_geno_440samples_renamedsamples_sorted.vcf.gz "


process pseudobulk {
    executor 'pbspro'
    clusterOptions = '-lselect=1:ncpus=40:mem=480gb -l walltime=5:00:00'
    conda "$params.conda_env_path"
    publishDir "${params.outdir}/MatrixEQTL_IO", mode: 'copy'

    input:
    path matrix
    path rscript

    output:
    path "*_pseudobulk.csv", emit: pseudobulk_results
    path "*_gene_locations.csv", emit: gene_locations
    

    script:
    """
    Rscript ${baseDir}/cellQTL_scripts/expression/generate_pseudobulk.r \
    --input_seurat_files ${seuratdir}/BRYOIS_92_AD_final_annotated_decontXfiltered.rds \
        ${seuratdir}/BRYOIS_92_MS_final_annotated_decontXfiltered.rds \
        ${seuratdir}/2023-08-23_MRC_60_final_annotated_decontXfiltered.rds \
        ${seuratdir}/2023-08-23_ROCHE_MPD_92_final_annotated_decontXfiltered.rds \
        ${seuratdir}/2023-08-23_MATTHEWS_final_annotated_decontXfiltered.rds \
    --min_cells ${params.min_cells} \
    --indiv_column ${params.indiv_column} \
    --celltype_column ${params.celltype_column} \
    --assay ${params.assay} \
    --script_dir ${baseDir}/cellQTL_scripts/expression/ \

    """
}

process generateGenotype {
    module "anaconda3/personal"
    conda "$params.conda_env_path"

    executor 'pbspro'
    clusterOptions = '-lselect=1:ncpus=20:mem=240gb -l walltime=3:00:00'
    publishDir "${params.outdir}/MatrixEQTL_IO", mode: 'copy'

    input:
    path vcf

    output:
    path "genotype_012mat.csv", emit: genotype_results
    path "snp_chromlocations.csv", emit: snp_chromlocations

    script:
    """
    Rscript ${baseDir}/cellQTL_scripts/genotype/generate_genotype_matrix.r \
    --script_dir ${baseDir}/cellQTL_scripts/genotype/ \
    --vcf $vcf \
    --ncores 10
    """
}

workflow pseudobulk_wf {
    pseudobulk(params.matrix)
}

workflow genotype_wf {
    generateGenotype(params.vcf)
}

workflow {
    // pseudobulk_wf()
    genotype_wf()
}


process eQTL {
    executor 'pbspro'
    clusterOptions = '-lselect=1:ncpus=40:mem=480gb -l walltime=5:00:00'
    conda '/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS'

    input:
    path pseudobulk_output
    path genotype_output
    path rscript

    output:
    path "eQTL_output", emit: eqtl_output

    script:
    """
    Rscript $rscript $pseudobulk_output $genotype_output
    """
}

workflow eqtl_wf {
    eQTL(params.rscript3)
}