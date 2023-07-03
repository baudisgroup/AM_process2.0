/* 
 * Extract probe values for total copy number
 * 
 */
process CRMAv2{

    tag {"$seriesName"}

    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  memory
    val  cleanup
    
    output:
    val "$seriesName"


    script:       
    """
    CRMAv2.R -datadir $datadir -workdir $workdir -seriesName $seriesName -force $force -memory $memory -cleanup $cleanup 
    """
}



