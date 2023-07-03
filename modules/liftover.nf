/* 
 * Map probes from hg19 to hg38 
 * 
 */
process liftover{
    tag {"$seriesName"}
    
    input:
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  liftover_loc 
    val  probe_file
    
    output:
    val "$seriesName"

    script: 
    """
    mkdir -p $datadir/processed/logs/liftover-log    


    segmentLiftover.py -i $datadir/processed/$seriesName -o $datadir/processed/$seriesName -c hg19ToHg38 -pi $probe_file -po probes,cn,hg38.tsv -l $liftover_loc --log_path $datadir/processed/logs/liftover-log --force $force
    """
}


