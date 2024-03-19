// enabling nextflow DSL v2
nextflow.enable.dsl=2

process PdbfixerMutants {
    input:
        path incsv
        path inpdb
    output:
        path '*_fixed.pdb', emit: fixed_pdbs
        path '*_reformat.csv', emit: csv_reformat
    shell:
    """
    mutant_maker.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}

process OpenmmMinimise {
    input:
        path fixed_pdbs
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
    PdbfixerMutants(inpath_ch, incsv_ch)
    OpenmmMinimise(PdbfixerMutants.out.fixed_pdbs)
    OutputData(PdbfixerMutants.out.csv_reformat, OpenmmMinimise.out.data)
}

