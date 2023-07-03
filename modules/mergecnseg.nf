/* 
 * merge segments for noisy profiles
 * 
 */

process mergecnseg{
    tag {"$seriesName"}
    
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  genome

    output:
    val "$seriesName"

    script: 
    """
    mergecnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome
    """
}


