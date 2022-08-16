import os
import subprocess
import pandas as pd
import numpy as np
#from qqman import qqman
import matplotlib.pyplot as plt
import adjustText
from adjustText import adjust_text
import seaborn as sns
from qmplot import qqplot


#os.system("sh test.sh")
#subprocess.call(['sh', './test.sh'])
#subprocess.call(['tools/plink', '--file', 'tools/toy', '--freq', '--out', 'tools/toy_analysis'])

class Genetics:
    def __init__(self, version=1.9, inpath='./', indt='', outpath='./', outdt='', starter='/mnt/tools/plink --file', finisher='--out', imputed=False):
        
        if version==2:
            plk="plink2"
        else:
            plk='plink'
        
        if inpath[len(inpath)-1]=='/':
            self.inpath=inpath
        else:
            self.inpath=inpath+'/'
        
        self.indt=indt
        self.outpath=outpath
        suffix=''
        
        if outpath[len(outpath)-1]=='/':
            self.outpath=outpath
        else:
            self.outpath=outpath+'/'
        
        self.outdt=outdt
        
        if os.path.exists(inpath+indt+'.bed'):
            self.starter='/mnt/tools/'+plk+' --bfile'
        elif os.path.exists(inpath+indt+'.ped'):
            self.starter='/mnt/tools/'+plk
        elif os.path.exists(inpath+indt+'.vcf'):
            self.starter='/mnt/tools/'+plk+' --vcf'
            suffix='.vcf'
        elif os.path.exists(inpath+indt+'.vcf.gz'):
            self.starter='/mnt/tools/'+plk+' --vcf'
            suffix='.vcf.gz'
        elif os.path.exists(inpath+indt+'.gen'):
            self.starter='/mnt/tools/'+plk+' --data'
        elif os.path.exists(inpath+indt+'.gen.gz'):
            self.starter='/mnt/tools/'+plk+' --gen'
        else:
            print('Plink file does not exist')
        self.finisher=finisher
        print(self.starter)
        print(self.finisher)
        #self.initCode=' '.join([self.starter, self.inpath+self.indt+suffix, '--keep-allele-order'])
        if imputed:
            self.initCode=' '.join([self.starter, self.inpath+self.indt+suffix+' dosage=HDS'])
        else:
            self.initCode=' '.join([self.starter, self.inpath+self.indt+suffix])
        self.endCode=' '.join([self.finisher, self.outpath+self.outdt])
    
    def Plink(self, extraClause='', snpFilter='--maf 0.01 --geno 0.05 --hwe 0.000001', applySnpFilter=0):
        
        if applySnpFilter==0:
            codebase=' '.join([self.initCode, extraClause, self.endCode])
        else:
            codebase=' '.join([self.initCode, snpFilter, extraClause, self.endCode])
        os.system(codebase)
        return codebase
    
    def calcFreq(self, extraClause=''):
        codebase=' '.join([self.initCode, extraClause, '--freq',self.endCode])
        os.system(codebase)
        return codebase
    
    def openLog(self):
        f=open(self.outpath+self.outdt+'.log','r')
        contents=f.read()
        print(contents)
        f.close()
        
    def loadResult(self, suffix='', header=None, commentSymbol='#'):
        if len(suffix)==0:
            suffix=suffix
        elif suffix[0]=='.':
            suffix=suffix
        elif suffix[0]!='.':
            suffix='.'+suffix
        print(self.outpath, self.outdt, suffix)
        inname=self.outpath+self.outdt+suffix
        dt=pd.read_csv(inname, sep=',|\t|\s+',engine='python',comment=commentSymbol,header=header)
        return dt
    
    def writeResult(self,outname,header=True):
        outname=self.outpath+outname
        self.dt.to_csv(outname,header=header)
    
    def manhattan(dt, path='./', mantitle='_Manhattan_Plot', label=True, sigPval=7.3, sugSigPval=5):
        running_pos = 0
        cumulative_pos = []

        for chrom, group_df in dt.groupby('#CHROM'):  
            cumulative_pos.append(group_df['POS'] + running_pos)
            running_pos += group_df['POS'].max()

        dt['cumulative_pos'] = pd.concat(cumulative_pos)
        
        dt['SNP number'] = dt.index
        
        plt.figure(figsize=(12, 10))
        #sns.set_context("paper", rc={"font.size":14,"axes.titlesize":12,"axes.labelsize":8})
        sns.set(font_scale = 1.5)
        g=sns.relplot(
            data = dt,
            x = 'cumulative_pos',
            y = 'LOG10_P',
            aspect = 2.5,
            hue = '#CHROM',
            palette = 'dark',
            linewidth=0,
            legend=None
        ).set(title=mantitle)
        
        g.ax.set_xlabel('Chromosome')
        g.ax.set_xticks(dt.groupby('#CHROM')['cumulative_pos'].median())
        g.ax.set_xticklabels(dt['#CHROM'].unique())
        g.ax.axhline(sigPval, linestyle='--', linewidth=1, color='r')
        g.ax.axhline(sugSigPval, linestyle='--', linewidth=1)
        if label:
            try:
                annotations = dt[dt['LOG10_P'] >sugSigPval].apply(lambda p : g.ax.annotate(p['ID'], (p['cumulative_pos'], p['LOG10_P'])), axis=1).to_list()
                adjust_text(annotations, arrowprops = {'arrowstyle' : '->', 'color' : 'blue'})
            except:
                print('No SNPs reached Genome wide significannce.')
        #fig = g.get_figure()
        g.fig.savefig(path+mantitle+'.png', dpi=250)
        
        return g

    
    def assoc(self, modeltype='linear', phefile='', covfile='', phenames=None, mac=20, covnames=None, scalePhe=False, scaleCov=False, extraClause='', log10Sig=7.3,log10SugSig=5, label=True, sampling=False):
        
        if len(phefile)==0:
            print("No phenotype file is specified.")
            return
        model='--glm "hide-covar" "omit-ref"'
        macFilt='--mac '+str(mac)
        incldPhe='--pheno '+phefile
        
        
        if not covnames:
            print('No covariates included.')
        else:
            if len(covfile)==0:
                print('No covar file specified., going to skip.')
                incldCov=''
                covnames=''
            else:
                if type(covnames) is str:
                    incldCov='--covar '+covfile
                    covnames='--covar-name '+covnames
                elif type(covnames) is list:
                    incldCov='--covar '+covfile
                    covnames='--covar-name '+', '.join(covnames)
        if not phenames:
            print("Phenotype name is missing")
            return
        elif type(phenames) is str:
            phenames=[phenames]
        
        scalePheFlag=0
        
        if scalePhe==False and scaleCov==True:
            scaling='--covar-variance-standardize'
        elif scalePhe==True and scaleCov==True:
            scaling='--variance-standardize'
        elif scalePhe==True and scaleCov==False:
            scalePheFlag=1
        else:
            scaling=''
        
        for x in phenames:
            phe2test='--pheno-name '+x
            if scalePheFlag==1:
                scaling='--covar-variance ' + x
            orioutdt=self.outdt
            oriendCode=self.endCode
            self.outdt=self.outdt+'_'+x
            self.endCode=self.endCode+'_'+x
            codebase=' '.join([self.initCode, model, macFilt, incldPhe, phe2test, incldCov, covnames, scaling, extraClause, self.endCode])
            print(codebase)
            os.system(codebase)
            assocRes=self.loadResult(suffix=x+'.glm.'+modeltype,header=0, commentSymbol=None)
            if sampling is not False:
                assocRes=assocRes.sample(10000)
            assocRes['LOG10_P']=-np.log10(assocRes['P'])
            sugSig=assocRes[assocRes['LOG10_P']>=log10SugSig]
            gwsig=assocRes[assocRes['LOG10_P']>=log10Sig]
            sugSig.to_csv(self.outpath+self.outdt+'_SugSigSNPs.tsv',index=False,sep='\t')
            gwsig.to_csv(self.outpath+self.outdt+'_GWSigSNPs.tsv',index=False,sep='\t')
            assocRes=assocRes.dropna(how="any", axis=0)
            #figure, axes = plt.subplots(nrows=1, ncols=2, figsize = (21,8))
            Genetics.manhattan(dt=assocRes,path=self.outpath, mantitle=self.outdt+'_Manhattan_plot', sigPval=log10Sig, sugSigPval=log10SugSig)
            qqplot(data=assocRes["P"], figname=self.outpath+self.outdt+"_QQ_plot.png",is_show=True,title=x,dpi=100)
            #qqman.qqplot(assocRes,title=x, col_p='LOG10_P')
            #figure.tight_layout()
            #plt.savefig(self.outpath+x+'_Manhattan_QQplot.png',format="png")
            self.outdt=orioutdt
            self.endCode=oriendCode
        return gwsig,sugSig,assocRes
        
            
        
class Shell:
    
    def __init__(self, inpath='./', indt='', outpath='./', outdt=''):
        
        if inpath[len(inpath)-1]=='/':
            self.inpath=inpath
        else:
            self.inpath=inpath+'/'
        self.indt=indt
        if outpath[len(outpath)-1]=='/':
            self.outpath=outpath
        else:
            self.outpath=outpath+'/'
        self.outdt=outdt
    
    def runShell(self, starter='cat', code='', order=1):
        """
        Run Shell script through python.
        
        Key parameters:
        starter: default is cat, you can provide command such as zcat, gunzip, zip.
        
        code: body of the code.
        
        Access to class parameters:
        self.inpath: 
        self.indt:
        self.outpath:
        self.outdt:
        
        """
        if order==1:
            codebase=' '.join([starter, self.inpath+self.indt, '|',code, '>', self.outpath+self.outdt])
        elif order==2:
            codebase=' '.join([starter, code, self.inpath+self.indt, '>', self.outpath+self.outdt])
        print(codebase)
        os.system(codebase)
        return codebase
    
    def Command(self, starter='', bodyCode='', endCode=''):
        """
        Run general command for any software through python
        Key paramters:
        starter: Beginning of the code, usuall is to call the software, for example, bcftools
        bodyCode: any additional parameters for the software after starter and input data
        endCode: output options
        
        Access to class parameters:
        self.inpath: 
        self.indt:
        self.outpath:
        self.outdt:
        
        """
        codebase=' '.join([starter, self.inpath+self.indt, bodyCode, self.outpath+self.outdt, endCode])
        os.system(codebase)
        return codebase