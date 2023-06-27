/* 
 * segmentation for total copy number 
 * 
 */
process cnseg{
    tag {"$seriesName"}
    
    input:
    val  previous_ready
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  undosd
    val  genome
    
    output:
      val "done"


    script:       
    """
    cnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -force $force -undosd $undosd -genome $genome 
    """
}



