// enabling nextflow DSL v2
nextflow.enable.dsl=2

process pdbfixer-mutants {
    input:
        path incsv
        path inpdb
    output:
        path '*_fixed.pdb', emit: fixed_pdbs
    shell:
    """
    mutant_maker.py --incsv $incsv --from-col ${params.csv.col} --in-pdb $inpdb ${params.mutant.maker.args}
    """
}

process openmm-minimise {
    input:
       path fixed_pdbs
    output: 
       path '*_unfolded.pdb', emit: unfolded_pdbs
       path '*_folded.pdb', emit: folded_pdbs
    shell:
    """
    openmm-minimise.py --i $fixed_pdbs
    """
}


workflow {
    inpath_ch = channel.fromPath("${params.inputFile}")
    incsv_ch = channel.fromPath("${params.inputCsv}")
    pdbfixer-mutants(inpath_ch, incsv_ch)
    openmm-minimise(pdbfixer-mutants.out.fixed_pdbs)
}

