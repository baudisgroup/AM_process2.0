/* 
 * Extract probe values for allele-specific copy number
 * 
 */
process ACNE{
    
    tag {"$seriesName"}
        
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  memory
    val  cleanup
    
    output:
      val "done"


    script:       
    """
    ACNE.R -datadir $datadir -workdir $workdir -seriesName $seriesName -force $force -memory $memory -cleanup $cleanup 
    """
}



