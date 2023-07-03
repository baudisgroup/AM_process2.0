# AM_process2.0

Updated pipelines for processing raw data from 9 SNP array platforms including Mapping 10K (GPL2641), Mapping 50K (Hind240 and Xba240) (GPL2004 and GPL2005), Mapping 250K (Nsp and Sty) (GPL3718 and GPL3720), Genomewide SNP (5.0 and 6.0) (GPL6894 and GPL6801, respectively), CytoScan (750K and HD) (GPL18637 and GPL16131) 

## Pre-requisites

###  Clone the repositoty

```
git clone https://github.com/baudisgroup/AM_process2.0.git
```

### Install nextflow 

Prerequisites for nextflow: Java 11 or later is required

```
cd AM_process2.0
curl -s https://get.nextflow.io | bash
```
This will create the `nextflow` main executable file in the current directory. 

Make scripts executable on your system 

```
chmod +x ./bin/*
```

### Install required packages

The pipeline can be run with Docker or local Mac OSX system. The installation is different for execution way.

#### Local execution

##### R libraries

Code in terminal

```
# for installation of mongolite

brew install openssl
```
Code in R

```
install.packages(c('BiocManager','remotes','aroma.affymetrix','ACNE','pastecs','genlasso','R.utils','matrixStats','tibble','plyr','mongolite'))

remotes::install_github('baudisgroup/LabelSeg')

BiocManager::install(c('aroma.core','DNAcopy','affxparser'))

install.packages('sfit',repos = 'https://henrikbengtsson.r-universe.dev')
```

##### Python libraries

Code in terminal

```
pip3 install click pandas numpy
```

#### Docker 

##### Install Docker 

Guide: [https://docs.docker.com/get-docker/](https://docs.docker.com/get-docker/)

##### Pull the following Docker image

```
docker pull hj12345/am_processv2.0
```

### Setup local MongoDB database

Guide: see internal group website `Computer Setup` and `Databases & MongoDB Cheats`

## Organization of data directory and working directory

Data directory is the directory where raw data and output data are stored. Working directory is the directory where the nextflow pipeline code resides. You must run this pipeline under the working directory. 

### data dir

The data directory structure is strict due to requirements of aroma project. See this [link](https://aroma-project.org/setup) for details.

The necessary directories can be found in `switchdrive/baudisgroup/arraymapPipeline/AM_process2.0/datadir`

```
data_dir/
├── PlatformInfo                                 # annotation files for each array platform (must keep)
├── ReferenceFile                                # processed data for external reference samples (must keep)
├── annotationData                               # annotation files for each array platform (must keep)
│   └── chipTypes
│       └── <chip type>
│           └──CDF and other annotation files
├── plmData                                      # processed data for external reference samples (must keep)
├── probeData                                    # intermediate directory when processing probe data
├── processed                                    # OUTPUT data directory
│   ├── logs                                       # log files
│   ├── data_quality_report                        # metrics of labeledsegment profiles and record of calibration steps 
│   │   └── cnseg
│   │       ├── <series1>.txt
│   │       └── <series2>.txt
│   ├── seriesFreq                                 # CNV frequency across samples in non-Progenetix series (used in calibration)
│   │   └── <series1>-freq.rds
│   ├── <series1>                                  # calling results of specific series
│   │   ├── <sample1>
│   │   │   ├── probes,cn.tsv
│   │   │   ├──	probes,cn,hg38.tsv
│   │   │   ├── segments,cn,provenance.tsv
│   │   │   ├── segments,cn.tsv
│   │   │   ├── labelsegments,cn.tsv
│   │   │   └── label_cnsegment.pdf
│   │ 	└── <sample2>
│   └── <series2>
│ 
└── rawData                                      # INPUT data directory                                        
    ├── <series1>
    │   ├── <chip type>
    │   │   ├── <sample1>.CEL
    │   │   └── <sample2>.CEL
    │   └── <chip type>   
    │       └── <sample3>.CEL
    └── <series2>
```

### work dir 

* `main.nf`: This is the main script. There are several parameters you can set directly by modifying the script or specify by command in terminal.

* `nextflow.config`: when running with docker, you need to modify `nextflow.config` file like this, making it possible to modify files outside of container. Other items in the config file can be left unchanged.

```
docker.runOptions='-v <your_datadir>:/app/datadir -v <your_workdir>:/app/workdir'
```

## Run the pipeline

**Parameters**

```
--datadir: location of data directory. Use absolute path and "/Users/usrname/..." instead of "~" in Mac OS.

--workdir: location of working directory. Use absolute path and "/Users/usrname/..." instead of "~" in Mac OS

--seriesName: series (GSE-XXX) to be analysed

--force: logical value to determine whether to continue to run and overwrite if there exists correponding output files. It is applied to some time-consuming calling steps including probe extraction, liftover, and segmentation. Default is false.

--memory: memory usage option in probe extraction. Generally for at most 8 threads, RAM 64GB -> 50, 96GB -> 80, 128GB -> 100. Default is 50.

--cleanup: logical value to determine whether to clean up intermediate files in the `plmData/` and `probeData/`. Default is true.

--undosd: a parameter for circular binary segmentation of relative copy numebr probe files (from `DNAcopy` R pcakge). Default is 1.

--genome: Version of reference genome . Because platform annotation files are based on hg19. The extracted probe files will be hg19. If specified as "hg38", this pipeline will convert the probes from hg19 to hg38 first and the following generated segment profiles will be hg38. If specified as "hg19", all of the output files are hg19. Default is "hg38".

--docker: logical value to determine whether to run with docker. Default is false.
```

### Local execution

```
cd <workdir>
./nextflow run main.nf --datadir <datadir> --workdir <workdir> --seriesName <series1>,<series2>
```

Example: `./nextflow run main.nf --datadir /Users/hangjia/datadir --datadir /Users/hangjia/workdir --seriesName GSE123,GSE456`

Default values for parameters can be set in `main.nf`. If in this case, you don't need to specify in terminal otherwise you want to set different values.

### Docker 

Before running, make sure you have modified the `nextflow.config` file and docker engine is running on your computer. 

```
cd <workdir>
./nextflow run main.nf --datadir <datadir> --workdir <workdir> --seriesName <series1>,<series2> --docker true -with-docker
```

Note: Runtime with Docker is much slower than local execution. Better to use if you just want to analyse individual series.

## Running tips

### Trace

* `-with-trace`: generate execution trace file (more useful for local execution) 
* `-with-report`: generate execution report (more useful for Docker execution)
* `-with-timeline`: generate timeline report

### Error strategy

`process.errorStrategy` option in `nextflow.config` file allows you to define how an error condition is managed by all processes. "terminate" (for debugging) and "ignore" (for batch processing) are usually used.




