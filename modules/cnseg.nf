/* 
 * segmentation for total copy number 
 * 
 */
process cnseg{
    tag {"$seriesName"}
    
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  undosd
    val  genome
    
    output:
    val "$seriesName"


    script:       
    """
    cnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -force $force -undosd $undosd -genome $genome 
    """
}



