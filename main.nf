/* 
 * pipeline input parameters 
 */
params.datadir = '/Users/hangjia/Data/test-pipeline/work_dir6'
params.workdir = '/Users/hangjia/switchdrive/Private/work/AM-nextflow'
params.seriesName  
params.force = false 
params.memory = 50
params.cleanup = true
params.undosd = 1
params.genome= "hg38"
params.docker = false


log.info """\
         CNV CALLING PIPELINE    
         ===================================
         datadir       : $params.datadir
         workdir       : $params.workdir
         series        : $params.seriesName
         force         : $params.force 
         memory        : $params.memory
         cleanup       : $params.cleanup
         undosd_cnseg  : $params.undosd
         genome        : $params.genome
         docker        : $params.docker
         """
         .stripIndent()

 

 include {CRMAv2} from './modules/CRMAv2.nf'
 include {liftover} from './modules/liftover.nf'
 include {cnseg} from './modules/cnseg.nf'
 include {adjustcnmedian} from './modules/adjustcnmedian.nf'
 include {labelcnseg} from './modules/labelcnseg.nf'
 include {mergecnseg} from './modules/mergecnseg.nf'
 include {calcnseg} from './modules/calcnseg.nf'
 include {mecan4cna} from './modules/mecan4cna.nf'





workflow {

    seriesName = Channel.of(params.seriesName.split(',')) 

    if (params.docker){
        datadir = '/app/datadir'
        workdir = '/app/workdir' 
        liftover_loc = 'liftOver'
    } else{
        datadir = params.datadir
        workdir = params.workdir
        liftover_loc = "$workdir/bin/liftOver"
    }    
               
    CRMAv2(datadir, workdir, seriesName, params.force, params.memory, params.cleanup)
    if (params.genome == "hg38"){
        liftover(datadir, workdir, CRMAv2.out, params.force, liftover_loc,"probes,cn,hg19.tsv")
        cnseg(datadir,  workdir, liftover.out, params.force, params.undosd, params.genome)
    } else if (params.genome == "hg19"){
        cnseg(datadir, workdir, CRMAv2.out, params.force, params.undosd, params.genome)        
    } else{
        error "input genome is not supported"
    }

    adjustcnmedian(datadir, workdir, cnseg.out, "True", params.genome)
    labelcnseg(datadir, workdir, adjustcnmedian.out, params.genome)
    mergecnseg(datadir, workdir, labelcnseg.out, params.genome)
    calcnseg(datadir, workdir, mergecnseg.out, params.genome, params.docker)
    mecan4cna(datadir, workdir, calcnseg.out, params.genome)
    
    }





    
