
nextflow.enable.dsl=2

params.rscript1 = "path/to/your/first/rscript"
params.rscript2 = "path/to/your/second/rscript"
params.rscript3 = "path/to/your/third/rscript"

params.matrix = "path/to/your/matrix"
params.vcf = "path/to/your/vcf"

params.outdir="/rds/general/user/ah3918/ephemeral/cellCAUSAL_TEST/"
process pseudobulk {
    executor 'pbspro'
    clusterOptions = '-lselect=1:ncpus=40:mem=480gb -l walltime=5:00:00'
    conda '/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS'

    input:
    path matrix
    path rscript

    output:
    path "*_pseudobulk.csv", emit: pseudobulk_results
    path "*_gene_locations.csv", emit: gene_locations

    script:
    """
    Rscript $rscript $matrix
    """
}

process generateGenotype {
    executor 'pbspro'
    clusterOptions = '-lselect=1:ncpus=20:mem=240gb -l walltime=3:00:00'
    conda '/rds/general/user/ah3918/home/anaconda3/envs/OSIRIS'
    publishDir "${params.outdir}/MatrixEQTL_IO", mode: 'copy'

    input:
    path vcf
    path rscript

    output:
    path "genotype_012mat.csv", emit: genotype_results
    path "snp_chromlocations.csv", emit: snp_chromlocations

    script:
    """
    Rscript ${baseDir}/cellQTL_scripts/genotype/generate_genotype_matrix.r \
    --script_dir ${baseDir}/cellQTL_scripts/genotype \
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