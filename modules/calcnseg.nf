
/* 
 * calibrate CNV calling based on prior knowledge 
 * 
 */

process calcnseg{
    tag {"$seriesName"}
    
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  genome
    val  docker 

    output:
    val "$seriesName"

    script: 
    """
    calcnseg.R -datadir $datadir -workdir $workdir -seriesName $seriesName -genome $genome -docker $docker
    """
}


