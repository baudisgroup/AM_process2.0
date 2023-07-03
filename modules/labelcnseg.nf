/* 
 * label and evaluate segments for total copy number 
 * 
 */
process labelcnseg{
    
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
    labelcnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome 
    """
}


