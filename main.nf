/* 
 * pipeline input parameters 
 */
params.datadir 
params.workdir 
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

    if (params.docker){
        datadir = '/app/datadir'
        workdir = '/app/workdir' 
        liftover_loc = 'liftOver'
    } else{
        datadir = params.datadir
        workdir = params.workdir
        liftover_loc = "$workdir/bin/liftOver"
    }
    
    seriesName = Channel.from(params.seriesName.split(','))
    
    CRMAv2(datadir, workdir, seriesName, params.force, params.memory, params.cleanup)

    if (params.genome == "hg38"){
        liftover(CRMAv2.out.collect(), datadir, workdir, seriesName, params.force, liftover_loc)
        cnseg(liftover.out.collect(),datadir, workdir, seriesName, params.force, params.undosd, params.genome)
    } else if (params.genome == "hg19"){
        cnseg(CRMAv2.out.collect(),datadir, workdir, seriesName, params.force, params.undosd, params.genome)        
    } else{
        error "input genome is not supported"
    }

    adjustcnmedian(cnseg.out.collect(), datadir, workdir, seriesName, "True", params.genome)
    labelcnseg(adjustcnmedian.out.collect(), datadir, workdir, seriesName, params.genome)
    mergecnseg(labelcnseg.out.collect(), datadir, workdir, seriesName, params.genome)
    calcnseg(mergecnseg.out.collect(), datadir, workdir, seriesName, params.genome, params.docker)
    mecan4cna(calcnseg.out.collect(), datadir, workdir, seriesName)
}

    
