/* 
 * Apply Mecan4CNA to calibrate logR values 
 * 
 */

process mecan4cna{

    tag {"$seriesName"}
    
    input:
    val  previous_ready
    val  datadir 
    val  workdir
    val  seriesName

    output:
    val  "done"

    script: 
    """
    dir="$datadir/processed/$seriesName"
    filename="labelsegments,cn.tsv"

    mapfile=\$(find "\$dir" -type f -name "\$filename")
    
    for file in \$mapfile; do
        python3 $workdir/bin/mecan4cna.py -i "\$file"        
    done

    """
}

