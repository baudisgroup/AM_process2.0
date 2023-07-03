/* 
 * adjust median for total copy number
 * 
 */
process adjustcnmedian{
    tag {"$seriesName"}
    
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  adjust_probe
    val  genome

    
    output:
    val "$seriesName"


    script:       
    """
    adjustcnmedian.R -datadir $datadir -workdir $workdir -seriesName $seriesName -adjust_probe $adjust_probe -genome $genome 
    """
}