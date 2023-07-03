/* 
 * Apply Mecan4CNA to calibrate logR values 
 * 
 */

process mecan4cna{

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
    dir="$datadir/processed/$seriesName"
    filename="labelsegments,cn,$genome\.tsv"

    mapfile=\$(find "\$dir" -type f -name "\$filename")
    
    for file in \$mapfile; do
        python3 $workdir/bin/mecan4cna.py -i "\$file"        
    done

    """
}

