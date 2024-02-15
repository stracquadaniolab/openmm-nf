// enabling nextflow DSL v2
nextflow.enable.dsl=2

process openmm {

    publishDir "${params.resultsDir}",  mode: 'copy'

    input:
        path inpdb
        val  temp
        val  pres
        val  solvmol
        val  nptrun
        val  nvtrun
        val  reprate
    output:
        path 'traj.pdb'
        path 'md_log.txt', emit: md_log
    
shell:
    """
        openmm-runner.py $inpdb $temp $pres $solvmol $nptrun $nvtrun $reprate  traj.pdb md_log.txt
    """

}

process hydro_plot {
    input:
      path md_log
      val  skipsteps
    
    shell:
      """
      plot.py $md_log $skipsteps
      """
}

workflow {
    inpath_ch = channel.fromPath("${params.inputFile}")
    temp = Channel.value(300)
    pres = Channel.value(1)
    solvmol = Channel.value(1000)
    nvtrun = Channel.value(1000)
    nptrun = Channel.value(1000)
    reprate = Channel.value(100)
    skipsteps = Channel.value(100)
    openmm(inpath_ch, temp, pres, solvmol, nvtrun, nptrun, reprate)
    hydro_plot(openmm.out.md_log, skipsteps)
}

