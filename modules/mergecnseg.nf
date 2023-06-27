/* 
 * merge segments for noisy profiles
 * 
 */

process mergecnseg{
    tag {"$seriesName"}
    
    input:
    val  previous_ready
    val  datadir 
    val  workdir
    val  seriesName
    val  genome

    output:
    val "done"

    script: 
    """
    mergecnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome
    """
}


