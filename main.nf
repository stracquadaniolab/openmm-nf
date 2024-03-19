// enabling nextflow DSL v2
nextflow.enable.dsl=2

process PdbfixerMutants {
    publishDir "${params.resultsDir}/pdbfixer/", pattern: "*_fixed.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/pdbfixer/", pattern: "*_reformat.csv", mode: 'copy'
    input:
        path incsv
        path inpdb
    output:
        path '*_fixed.pdb'
        path '*_reformat.csv', emit: csv_reformat
    shell:
    """
    mutant_maker.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}

process OpenmmMinimise {
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "*folded.pdb", mode: 'copy'
    publishDir "${params.resultsDir}/openmm-minimise/", pattern: "data.csv", mode: 'copy'
    input:
        path "*_fixed.pdb"
    output: 
        path '*_unfolded.pdb', emit: unfolded_pdbs
        path '*_folded.pdb', emit: folded_pdbs
        path 'data.csv', emit: data
    shell:
    """
    openmm-minimise.py --i $fixed_pdbs
    """
}

process OutputData {
    publishDir "${params.resultsDir}/output_data/", pattern: "data_ΔΔG.csv", mode: 'copy'
    publishDir "${params.resultsDir}/output_data/", pattern: "data_ΔΔG-spearman.csv", mode: 'copy'
    input:
        path benchcsv
        path testcsv
    output:
        path 'data_ΔΔG.csv', emit: data_ΔΔG
        path 'data_ΔΔG-spearman.csv', emit: spearman
    shell:
    """
    output_data.py --bench $benchcsv --test= $testcsv
    """
}


workflow {
    inpath_ch = channel.fromPath("${params.inputFile}")
    incsv_ch = channel.fromPath("${params.inputCsv}")
    PdbfixerMutants(incsv_ch, inpath_ch)
    OpenmmMinimise()
    OutputData(PdbfixerMutants.out.csv_reformat, OpenmmMinimise.out.data)
}

