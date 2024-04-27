
params.rscript1 = "path/to/your/first/rscript"
params.rscript2 = "path/to/your/second/rscript"
params.rscript3 = "path/to/your/third/rscript"

params.matrix = "path/to/your/matrix"
params.vcf = "path/to/your/vcf"

process pseudobulk {
    input:
    path matrix from params.matrix
    path rscript from params.rscript1

    output:
    path "*_pseudobulk.csv" into pseudobulk_results
    path "*_gene_locations.csv" into gene_locations

    script:
    """
    Rscript $rscript $matrix
    """
}

process generateGenotype {
    input:
    path vcf from params.vcf
    path rscript from params.rscript2

    output:
    path "genotype_012mat.csv" into genotype_results
    path "snp_chromlocations.csv" into snp_chromlocations

    script:
    """
    Rscript $rscript $vcf
    """
}

process eQTL {
    input:
    path pseudobulk_output from pseudobulk_results
    path genotype_output from genotype_results
    path rscript from params.rscript3

    output:
    path "eQTL_output"

    script:
    """
    Rscript $rscript $pseudobulk_output $genotype_output
    """
}


process pseudobulk {
    input:
    path matrix from params.matrix
    path rscript from params.rscript1

    output:
    path "*_pseudobulk.csv" into pseudobulk_results
    path "*_gene_locations.csv" into gene_locations

    script:
    """
    Rscript $rscript $matrix
    """
}

process generateGenotype {
    input:
    path vcf from params.vcf
    path rscript from params.rscript2

    output:
    path "genotype_012mat.csv" into genotype_results
    path "snp_chromlocations.csv" into snp_chromlocations

    script:
    """
    Rscript $rscript $vcf
    """
}

process eQTL {
    input:
    tuple path(pseudobulk_output), path(gene_locations) from pseudobulk_and_gene_locations
    path genotype_output from genotype_results
    path snp_chromlocations from snp_chromlocations
    path rscript from params.rscript3

    output:
    path "eQTL_output"

    script:
    """
    Rscript $rscript $pseudobulk_output $genotype_output $snp_chromlocations $gene_locations
    """
}

workflow {
    pseudobulk(params.matrix, params.rscript1)
    generateGenotype(params.vcf, params.rscript2)
    pseudobulk_and_gene_locations = pseudobulk_results.zip(gene_locations)
    eQTL.each(pseudobulk_and_gene_locations, generateGenotype.out[0], generateGenotype.out[1], params.rscript3)
}