/* 
 * label and evaluate segments for total copy number 
 * 
 */
process labelcnseg{
    
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
    labelcnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome 
    """
}


