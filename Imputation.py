import os
import subprocess
import pandas as pd

class Phasing:
    def __init__(self, inpath='./', indt='', outpath='./', outdt='', initCode='/mnt/tmp/Eagle_v2.4.1/eagle', endCode='--outPrefix='):
        self.inpath=inpath
        self.indt=indt
        self.outpath=outpath
        self.outdt=outdt
        self.initCode=initCode
        self.endCode=endCode
    
    def phasePlink(self, build=38, nThreads=4):

        if build==38:
            geneMap='--geneticMapFile=/mnt/tmp/Eagle_v2.4.1/tables/genetic_map_hg38_withX.txt.gz' 
        elif build==19:
            geneMap='--geneticMapFile=/mnt/tmp/Eagle_v2.4.1/tables/genetic_map_hg19_withX.txt.gz' 
        else:
            print('Build version is not supported.')
            return
        logfile=self.outdt+'.log'
        threadNum='--numThreads='+str(nThreads)
        
        for chr in range(1,23):
            print('Current Chr is: ',str(chr))
            curChr='--chrom='+str(chr)
            outname=self.outdt+'_'+str(chr)
            logfile=self.outdt+'.log'
            codebase=' '.join([self.initCode, '--bfile='+self.inpath+self.indt, geneMap, curChr, threadNum, self.endCode+self.outpath+outname, '2>&1 | tee -a', self.outpath+logfile])
            os.system(codebase)
        return codebase
    
    def openLog(self):
        f=open(self.outpath+self.outdt+'.log','r')
        contents=f.read()
        print(contents)
        f.close()
    
    def pahseBCF(self, build=38, ref=False):
        return
        
        
class Impute:
    def __init__(self, inpath='./', indt='', outpath='./', outdt=''):
        self.inpath=inpath
        self.indt=indt
        self.outpath=outpath
        self.outdt=outdt