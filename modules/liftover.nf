/* 
 * Map probes from hg19 to hg38 
 * 
 */
process liftover{
    tag {"$seriesName"}
    
    input:
    val  previous_ready
    val  datadir 
    val  workdir
    val  seriesName
    val  force
    val  liftover_loc 
    
    output:
    val "done"

    script: 
    """
    mkdir -p $datadir/processed/logs/liftover-log    

    segmentLiftover.py -i $datadir/processed/$seriesName -o $datadir/processed/$seriesName -c hg19ToHg38 -pi probes,cn.tsv -po probes,cn,hg38.tsv -l $liftover_loc --log_path $datadir/processed/logs/liftover-log --force $force

    segmentLiftover.py -i $datadir/processed/$seriesName -o $datadir/processed/$seriesName -c hg19ToHg38 -pi probes,fracb.tsv -po probes,fracb,hg38.tsv -l $liftover_loc --log_path $datadir/processed/logs/liftover-log --force $force
    """
}


