/* 
 * reliability assessment for total copy number profiles 
 * 
 */
process CNARA{

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
    run_CNARA.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome 
    """
}