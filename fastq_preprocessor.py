#!/usr/bin/env python2
# Usage: (python) preprocessor.py /path/to/FASTQ_DIR
# Example: preprocessor.py /media/pw_synology3/PW_HiSeq_data/RNA-seq/Raw_data/testONLY/133R/BdPIFs-32747730/133E_23_DN-40235206
# Purpose: Download fastq files from the supplied path and 
#    combine them into R1.fastq and R2.fastq
# 
# Created:7  OCT 2016, Hui@SLCU map-RNA-seq.py
# Update: 11 Feb 2017. Hui@SLCU
# Update: 29 May 2018. Feng@SLCU preprocessor.py
# Update: 08 Apr 2020. Feng

import tempfile,subprocess
import os, sys, datetime, glob, re, collections
import multiprocessing as mp
import pandas as pd
import itertools
import re

class pyext(object):
    def df__iterdict(self, into=dict):
        '''
        ### See github issue https://github.com/pandas-dev/pandas/issues/25973#issuecomment-482994387
        '''
        it = self.iterrows()
        for index, series in it:
            d = series.to_dict(into=into)
            d['index'] = index
            yield d

class ptn(object):
    Rcomp=re.compile

    ridL = r'(^|.+/)(\d{1,4}[RCQ]{1,2}|SRR\d{7,8}'
    # ridL = r'(\d{1,4}[RC]'
    end = r'[_/\-]'
    runID = Rcomp('(?=%s)([^RCQ].*|$))'%(ridL,))
    runCond = Rcomp('%s%s.*)'%(ridL,end))
    ridPATH = Rcomp('%s%s.*)'%(ridL,end))
    runID_RNA = Rcomp('.*%s).*'%(ridL.replace('C',''),))

    srr = Rcomp('(?P<lead>SRR.*)_(?P<read>[012]).(?P<ext>.+)')
    baseSpace = Rcomp('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
    baseSpaceSimple = Rcomp('(?P<lead>.*)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
    sampleID = Rcomp('[_\/](S\d{1,3}|SRR\d{7,8})[_/\-\.]')

    BdAcc = Rcomp('(Bradi[\da-zA-Z]+)')
# shellexec = os.system
def paste0(ss,sep=None,na_rep=None,castF=unicode):
    '''Analogy to R paste0
    '''
    if sep is None:
        sep=u''
    
    L = max([len(e) for e in ss])
    it = itertools.izip(*[itertools.cycle(e) for e in ss])
    res = [castF(sep).join(castF(s) for s in next(it) ) for i in range(L)]
    res = pd.Series(res)
    return res
                 
def check_all_samples(d):
    '''
    [Deprecated] Use LeafFiles(DIR) instead, kept for a reference
    Recursively walking a directory tree
    TBM to return useful output
    '''
    for k in sorted(d.keys()):
        rpath = d[k]
        result = check_directory(rpath)
        if result != '':
            print('%s:%s [%s]' % (k, rpath, result))
            rpath = rpath.rstrip('/')
            if glob.glob(rpath) is []: # rpath does not exist or it is an incomplete path
                path_lst = glob.glob(rpath + '*/')
                assert len(path_lst) == 1,'Found multiple directories under %s'%rpath
                rpath = path_lst[0]
                d[k] = rpath.rstrip('/') # update rpath



def nTuple(lst,n,silent=1):
    """ntuple([0,3,4,10,2,3], 2) => [(0,3), (4,10), (2,3)]
    
    Group a list into consecutive n-tuples. Incomplete tuples are
    discarded e.g.
    
    >>> group(range(10), 3)
    [(0, 1, 2), (3, 4, 5), (6, 7, 8)]
    """
    if not silent:
        L = len(lst)
        if L % n != 0:
            print '[WARN] nTuple(): list length %d not of multiples of %d, discarding extra elements'%(L,n)
    return zip(*[lst[i::n] for i in range(n)])

def LinesNotEmpty(sub):
    sub = [ x for x in sub.splitlines() if x]
    return sub

def LeafFiles(DIR,check=False):
    ''' Drill down to leaf files of a directory tree if the path is unique.
    '''
    assert os.path.exists(DIR),'%s not exist'%DIR
    DIR = DIR.rstrip('/')
    if not os.path.isdir(DIR):
        return [DIR]
    else:
        cmd = 'ls -LR %s'%DIR
        res = subprocess.check_output(cmd,shell=1)
        res = re.split(r'([^\n]*):',res)[1:]
        it = nTuple(res,2,silent=0)
        DIR, ss = it[0];
        for dd,ss in it[1:]:
            NEWDIR, ALI = dd.rsplit('/',1)
            if check:
                assert NEWDIR == DIR, 'Next directory %s not contained in %s'%(dd,DIR)
            DIR = dd 
        res = [ '%s/%s'%(DIR,x) for x in LinesNotEmpty(ss)]
        return res                

retype = type(re.compile('hello, world'))
def revSub(ptn, dict):
    '''Reverse filling a regex matcher.
    Adapted from: https://stackoverflow.com/a/13268043/8083313
'''
    if isinstance(ptn, retype):
        ptn = ptn.pattern
    ptn = ptn.replace(r'\.','.')
    replacer_regex = re.compile(r'''
        \(\?P         # Match the opening
        \<(.+?)\>
        (.*?)
        \)     # Match the rest
        '''
        , re.VERBOSE)
    res = replacer_regex.sub( lambda m : dict[m.group(1)], ptn)
    return res

def write_log(fname, s):
    f = open(fname, 'a')
    f.write(s + '\n')
    f.close()

def gnuPara(cmd,debug=0,ncore = 6):
    '''
    [Deprecated] Bad and does not wait for tasks to finish
    '''
    tmp = tempfile.NamedTemporaryFile(delete=True) if not debug else open('temp.sh','w')
    with tmp as tmp:
        print cmd
        tmp.write(cmd)
        E = shellexec('parallel --gnu -j%d <%s &>>parallel.log'%(
            ncore,
            tmp.name
            )
        )
    return E

def mp_para(f,lst,ncore = 6):
    if ncore ==1:
        res = map(f,lst)
    else:
        p = mp.Pool(ncore)
        res = p.map_async(f,lst,)
        res = res.get(10000000)
        p.close()
        p.join()
    return res

datenow = lambda: datetime.datetime.now().strftime("%Y_%m_%d_%H:%M:%S")

#### Regex for downloaded .fastq(.gz) files
# PTN = re.compile('(?P<lead>.*)_S(?P<sample>\d{1,3})_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')
# PTN = re.compile('(?P<lead>.*)_L(?P<chunk>\d+)_R(?P<read>[012])_(?P<trail>\d{1,4})\.(?P<ext>.+)')


def shellexec(cmd,debug=0):
    print(cmd) 
    if not debug:
        return subprocess.call(cmd,shell=1)
#         return os.system(cmd)

def process_rna_sample(samplePATH, debug=0,checkMatch=0, timestamp=1, newDIR=None,
                       newDir= None,
#                        force=0,
                       NCORE=6,
                      moveRaw=1,rename=0,**kwargs):
    '''
    Pull together raw reads from an input folder
    Args:
        samplePATH: Folder of .fastq(.gz) fot. be processed
    Comment: Refactored based on Hui's map-RNA-seq.py process_rna_sample().    
    
    [ToDo]
        This function is way too messy to be understood
    '''
    del kwargs
    ### legacy
    if newDir is not None:
        newDIR = newDir
        
    #     return os.system('/bin/bash -c `%s`'%cmd)
    #     cmd = '/bin/bash -c `%s`'%cmd

    #     return subprocess.call(cmd,env=os.environ,cwd=os.getcwd(),
    #                           shell=True)
    

    samplePATH = samplePATH.rstrip('/')
    shellexec('echo $SHELL')
    
    
    ### Extract  RunID from samplePATH
    samplePATH = samplePATH.rstrip('/')
    sp = samplePATH.rsplit('/',2)
    DataAccPath = ridPath = '/'.join(sp[-2:])
    OLDPATH = sp[0]
    os.system('echo %s>OLDPATH' % OLDPATH)
    
    
    
    if checkMatch:
        #### OLD 
        pats = [ptn.runCond, 'SRR/SRR\d{7,8}', 'R\d{3,5}/S\d{1,5}']
        res = [ re.findall(pat, ridPath,) for pat in pats]
        assert max(map(len,res))==1,\
        '[ERROR] Cannot extract RunID from Accession:"{ridPath}"\
        \nRegex result:{res}\
        \nsamplePath:{samplePATH}'.format(**locals())
    
    print ('[ridPath]%s'%ridPath)


    
    def filterFastq(FILES):
        FILES = filter(re.compile('.*\.(fastq)(\.gz|)').match, FILES)
        return FILES

    #### Create a temporary directory 
    ODIR = os.getcwd()
    if newDIR is not None:
        temp_dir = newDIR
    else:
        DIR = [ridPath.replace('/','-')]
        if timestamp:
            DIR +=  [str(datenow())]
        DIR = '-'.join(DIR)
        temp_dir = DIR
    os.system('mkdir -p %s'%temp_dir)          
    os.chdir(temp_dir) #     shellexec('cd %s'%temp_dir)
    
    
    #### Download raw read .fastq from samplePATH
    if 1:
        FILES = glob.glob('%s/*' % samplePATH)
        FILES = sum(map(LeafFiles,FILES),[])
        FILES = filterFastq(FILES)
    #     ccmd = '%s/* -t %s'%(samplePATH,temp_dir) 
        ccmd = '%s -t ./'%(' '.join(FILES),) 
        cmd1 = 'cp -lr %s'%ccmd; 
        cmd2 = 'cp -r %s'%ccmd
        shellexec(cmd1) ==0 or shellexec(cmd2) 

        print ('[ODIR]%s'%ODIR )
        
        
        if 1:
            #### Parse .fastq filenames and assert quality checks
            if debug:
                FS = [x.rsplit('/')[-1] for x in  FILES]
                print( FS[:5] )
            else:
                FS = glob.glob('*')
                FS = filterFastq(FS)

            BUF = '\n'.join(FS)
            BUFHEAD = '\n'.join(FS[:5])
            ##### Process baseSpace files
            res = {
                x:len(re.findall( getattr(ptn,x),BUFHEAD)) for x in ['baseSpace',
                                                                     'baseSpaceSimple',
                                                                     'srr'] }

            if debug:
                print res.items()
            assert max(res.values())>0,'Cannot identify format of files:\n%s'%BUF
            patName, pat = [(x,getattr(ptn,x)) for x,y in res.items() if y > 0][0]
            PARSED = [dict(m.groupdict().items() + 
                           [('fname',m.group(0))])
                     for m in re.finditer(pat,BUF) ]
            meta = pd.DataFrame(PARSED)

            if patName =='baseSpace':
                meta = check_L004(meta)
                meta = meta.sort_values(['lead','read','chunk'])
            elif patName == 'baseSpaceSimple':
                meta = meta            
            elif patName=='srr':
                meta = meta
            meta  = meta__unzip(meta,debug=debug, NCORE=NCORE)
            meta = meta__concat(meta,debug=debug, NCORE=NCORE)
            if moveRaw:
                meta = meta__moveRaw(meta,debug=debug, NCORE=NCORE)
            if rename:
                meta = meta__rename(meta,debug=debug, NCORE=NCORE)


        if 1:
            ### Bookkeeping, not used
            meta.insert(0,'DataAccPath', DataAccPath)
            meta.insert(1,'DataAcc', DataAccPath.replace('/','-'))
            print ('[OLDDIR]',ridPath,os.system('echo %s | tee OLDDIR | tee DATAACC'%ridPath))
            meta__dumpFileArg(meta, debug=debug,NCORE=NCORE)
            print '[DONE!]:%s'%samplePATH
            meta.to_json('META.json',orient='records')
            meta.to_csv('META.csv')
            with open('FILE.json','w') as f:
                json.dump(next(pyext.df__iterdict(meta)), f)
                

        if debug:
            print( meta[['read','fname']] )
            assert 0
#             return meta
        else:
            pass
        
#         unzipAndConcat(meta)
#         exit(0)
#     except Exception as e:        
#         exc_type, exc_obj, exc_tb = sys.exc_info()
#         fname = os.path.split(exc_tb.tb_frame.f_code.co_filename)[1]
#         print(exc_type, fname, exc_tb.tb_lineno)
#         raise e
#     finally:
#     import os
    temp_dir = os.path.basename(os.getcwd())
    sys.stdout.write('%s\n'%temp_dir)
    if 1:
        os.chdir(ODIR)
    #### Stop here
    return temp_dir
def meta__dumpFileArg(meta,debug=None,NCORE=None):
    del debug,NCORE
    res = len(meta.read.unique())
    if len(meta)==2:
        paired = 1
        read1,read2 = meta.sort_values('read').fname
        res = '-1 {read1} -2 {read2}\n'.format(**locals())
    elif len(meta==1):
        paired = 0
        read1  = meta.fname[0]
        read2  = ''
        res = '-1 {read1}\n'.format(**locals())
#     res = '{read1} {read2}\n'.format(**locals())            
    with open('FILEARG','w') as f:
        f.write(res)
    return 

def check_L004(meta):
    g = meta.groupby(['lead','read'],as_index=0)
    ct = g.count()

    mout = meta.merge(ct[['lead','read','chunk']] ,on=['lead','read'],suffixes=['','_count'])
    idx = mout['chunk_count'] ==4
    if not idx.all():
        print '[WARN] following reads are discarded due to chunk_count != 4'
        print mout[~idx][['fname','chunk']] 
        mout = mout[idx]
    return mout
def meta__unzip(meta,debug=0,NCORE=None):
    assert NCORE is not None
    idx= [x.endswith('gz') for x in meta['fname']]
    if any(idx):
        #### unzip .gz where applicable
        mcurr = meta.iloc[idx]
        cmds = [cmd_ungzip(x) for x in mcurr['fname']]
        if debug:
            print '\n'.join(cmds[:1])
        else:
            mp_para(shellexec,cmds, ncore=NCORE)            
        #### Remove .gz in DataFrame accordingly
        meta.loc[idx,'ext'] = [ x.rstrip('.gz')  for x in mcurr['ext'] ]
        meta.loc[idx,'fname'] = [ x.rstrip('.gz')  for x in mcurr['fname'] ]
    return meta

def meta__moveRaw(meta,debug=0, NCORE=  None):
    assert NCORE is not None
    shellexec('mkdir -p raw/')
    cmds = ['mv %s -t raw/'%(' '.join(meta.fname))]
    meta['fname'] = meta.eval('@paste0([["raw/"],fname])')
    if debug:
        print '\n'.join(cmds[:1])
    else:
        mp_para(shellexec,cmds, ncore=NCORE)  
    return meta

def meta__rename(meta,debug=0, NCORE=None):
    assert NCORE is not None
#     shellexec('mkdir -p raw/')
    cmds = []
    for d in pyext.df__iterdict(meta):
        d['ofname'] = 'R{read}.{ext}'.format(**d)
        cmd = 'mv {fname} {ofname}'.format(**d)
        cmds += [cmd]
        meta.loc[d['index'],'fname'] = d['ofname']
#     cmds = ['mv %s -t raw/'%(' '.join(meta.fname))]
#     meta['fname'] = meta.eval('@paste0([["raw/"],fname])')
    if debug:
        print '\n'.join(cmds)
    else:
        mp_para(shellexec, cmds, ncore=NCORE)  
    return meta

def meta__concat(meta,debug= 0, NCORE=None):
    assert NCORE is not None
    ### Map metas to fnames after decompression 
    if 'chunk' not in meta.keys():
        return meta
    else:
        mapper = lambda x: revSub(ptn.baseSpace,x)
        meta['fname'] = meta.apply(mapper,axis=1)
        meta = meta.groupby(['lead','read']).apply(cmd_combineFastq).reset_index()

        print meta
    #     assert 0
    #     ofnames,cmds = zip([cmd_combineFastq(x[1]['fname']) for x in g])

        cmds = meta.cmd
        if debug:
            print '\n'.join(cmds[:1])
        else:
            mp_para( shellexec, cmds, ncore=NCORE)
    #     os.system('sleep 5;')
    return meta
def cmd_combineFastq(df,run=0):
    df = df.sort_values('fname')
#     fnames = sorted(list(fnames))
    fflat = ' '.join(df.fname)
#     d = ptn.baseSpace.match(fnames[0]).groupdict()
    lead,read,ext=df.head(1)[['lead','read','ext']].values.ravel()
    ofname = '{lead}_R{read}_raw.{ext}'.format(**locals())
    cmd = 'cat {fflat} >{ofname} ; sleep 0; rm {fflat} '.format(
        **locals())
#     df['cmd'] = cmd
#     df['ofname'] = ofname
    return pd.Series({'fname':ofname,'cmd':cmd})

# def cmd_combineFastq(fnames,run=0):
#     fnames = sorted(list(fnames))
#     d = ptn.baseSpace.match(fnames[0]).groupdict()
#     cmd = 'cat {IN} >{lead}_R{read}_raw.{ext} ; sleep 0; rm {IN} '.format(
#         IN=' '.join(fnames),
#                                                  **d)
#     return cmd

def cmd_ungzip(F,):
    cmd = 'gzip -d <{IN} >{OUT} ; sleep 0 ; rm {IN} '.format(IN=F,OUT=F.rstrip('.gz'))
    return cmd



import argparse
parser=argparse.ArgumentParser()
parser.add_argument('samplePATH',)
parser.add_argument('--timestamp',action='store_true')
parser.add_argument('--newDIR',default='.',)
parser.add_argument('--newDir',default=None,)
parser.add_argument('--moveRaw',default=0,type=bool)
parser.add_argument('--rename',default=0,type=bool)

parser.add_argument('--force',default = 0, type = int,) #### deprecated and not used
parser.add_argument('--checkMatch',default = 0, type = int,)
parser.add_argument('--debug',action='store_true')
parser.add_argument('--NCORE',default=6,type=int)
# 'raw-data/900R/S1'
# argparse.Aru

main = process_rna_sample
if __name__=='__main__':
    
#     NCORE = int(os.environ.get('NCORE',6))
#     print '[NCORE]=',NCORE
    # NCORE = 1
# assert len(sys.argv) >= 2,'''
#     Usage: (python) map-RNA-seq.py /path/to/folder/
#         The folder should contains raw reads in .fastq(.gz) format
# '''
    samplePATH = sys.argv[1]
    args = parser.parse_args()
#     NCORE = args.NCORE
    temp_dir = process_rna_sample(**vars(args))
    sys.exit(0)
