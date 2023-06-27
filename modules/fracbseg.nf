/* 
 * segmentation for allele-specific copy number 
 * 
 */
process fracbseg{
    tag {"$seriesName"}
    
    input:
    val  previous_ready
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  genome
    
    output:
      val "done"


    script:       
    """
    fracbseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -force $force -genome $genome
    """
}



