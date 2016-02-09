__author__ = 'stevebertolani'
import sys
import os
from os import path, environ
import tempfile
import argparse
import shutil
import glob
from hmmer.Parse import AlipidLine
import urllib2
from Bio import SeqRecord,SeqIO,AlignIO
from io import BytesIO
from hmmer import API
import gzip
import urllib
import json
import prody
import subprocess
# add in options argpars
import templates.Templates as templates
from geoutil.util import *
from rosetta import run
from Bio.SeqUtils.CheckSum import seguid
from random import randint
from time import sleep

class SmartRedirectHandler(urllib2.HTTPRedirectHandler):
    def http_error_302(self, req, fp, code, msg, headers):
        return headers



def filter_pid(mylist,cutoff=None):
    if cutoff == None: cutoff = 105;
    myfilteredlist = []
    for eachline in mylist:
        print eachline.pid
        if eachline.pid < cutoff:
            print "%s passed, pid to target is %f" %(eachline.name,eachline.pid)
            myfilteredlist.append(eachline)
    return myfilteredlist

def reorder_aln(myalnfile,name):
        '''
            Allows you to specify which aln (by name) should be at the top of the file
        '''
        with open(myalnfile,'rb') as fh, open('temp','wb') as outfh:
            myalignmentdict = SeqIO.to_dict( AlignIO.read(fh,'fasta') )
            target = myalignmentdict[name]
            outfh.write(">%s\n%s\n" %(target.id,target.seq))
            for i in myalignmentdict:
                if myalignmentdict[i].id != name:
                    outfh.write(">%s\n%s\n" %(myalignmentdict[i].id,myalignmentdict[i].seq))
        os.system("cp temp %s" %myalnfile)


def remove_dup_seqs(records):
    checksums = set()
    for record in records:
        checksum = seguid(record.seq)
        if checksum in checksums:
            print "Ignoring %s" % record.id
            continue
        checksums.add(checksum)
        yield record
        
def read_zipped_pdb_file(fn):
    lines = os.popen( 'zcat %s'%(fn),'r').readlines()
    return lines

def run_phmmer(seq,database):
    mydatabase = database
    res_params ={'output':'json','range':'1,500'}                   # don't forget to grab search params
    json_file_name = "%s.%s" %(database,res_params['output'])

    if path.isfile(mydatabase):
        if path.isfile(json_file_name):
            print 'It seems this is already completed'
            return json_file_name
            #
#    if not path.isfile(mydatabase+".gz"):
    print "Running a phmmer search on %s" %database
    method = "phmmer"

    search_params = {'seqdb':"%s" %database,'seq': ">%s\n%s" %(seq.id,seq.seq) }
    print "Using search parameters %s" %search_params

    opener = urllib2.build_opener(SmartRedirectHandler())
    urllib2.install_opener(opener);
    parameters = search_params
    enc_params = urllib.urlencode(parameters)
    request = urllib2.Request('http://www.ebi.ac.uk/Tools/hmmer/search/phmmer',enc_params)
    results_url = urllib2.urlopen(request).getheader('location')
    print "results url %s " %results_url
    enc_res_params = urllib.urlencode(res_params)
    modified_res_url = results_url + '?' + enc_res_params
    results_request = urllib2.Request(modified_res_url)
    data = urllib2.urlopen(results_request)
    download_seq_param = urllib.urlencode(  {'format':'fullfasta'}  )
    jobid = results_url.split('/')[-2]
    print "job id %s" %jobid
    download_url = 'http://www.ebi.ac.uk/Tools/hmmer/download/%s/score' %jobid + '?' + "%s" %download_seq_param
    download_request = urllib2.Request(download_url)
    print "downloading from %s" %download_url
    print "now trying to download the full sequences ..."
    datagz = urllib2.urlopen(download_request, timeout=500)

    with open('%s.gz' %database, 'w' ) as fh:
        fh.write(datagz.read())

    with open(json_file_name,'w') as fh:
        fh.write(data.read())
        print "Got the json file to write"

    if not path.isfile(database):
        output = gzip.open(database+".gz", 'rb')
        outfh = open('%s' %database,'wb')
        outfh.write( output.read() )
        outfh.close()
        
    return json_file_name

def read_json(jsonfile):
    json_data=open(jsonfile)
    data = json.load(json_data)
    print "Number of hits found: %s" %( len(data['results']['hits']) )
    mylist = []
    for i in range( 0, int(len( data['results']['hits'] ))  ):
        name = data['results']['hits'][i]['acc']
        mylist.append(name)

    json_data.close()
    return mylist

def align_templates_promals3d(infastafile):
    mycmd = '/home/bertolan/software/promals_package/bin/promals %s > %s.log' %(infastafile,infastafile) 
    print mycmd
    try:
        os.system(mycmd)
    except:
        print "error running promals3d aln"

def get_pid(name):
    mylist = []
    print 'checking for %s' %name
    cmd = ['esl-alipid','combo.msa.fasta.clean']
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in iter(p.stdout.readline,''):
        if line.startswith(name):
            myline = AlipidLine(line)
            print "%s %s" %(myline.name,str(myline.pid))
            mylist.append(myline)

    mytemplatesequences = filter_pid(mylist)
    return mytemplatesequences

def main():


    nofrags = False

    mymaindir = os.getcwd()

    mypidpath = mymaindir+'/pid'
    if not os.path.isdir(mypidpath):
        os.mkdir(mypidpath)

    mypdbstructurespath = mymaindir+'/pdbstructures'
    if not os.path.isdir(mypdbstructurespath):
        os.mkdir(mypdbstructurespath)

    mymodelspath = mymaindir+'/models'
    if not os.path.isdir(mymodelspath):
        os.mkdir(mymodelspath)

    myseqstructurespath = mymaindir+'/sequences'
    if not os.path.isdir(myseqstructurespath):
        os.mkdir(myseqstructurespath)

    homologuesdir = mymaindir+'/homologoussequences'
    if not os.path.isdir(homologuesdir):
        os.mkdir(homologuesdir)

    if not os.path.isdir(mymaindir+'/templates'):
        os.mkdir(mymaindir+'/templates')

    fasta_seq_in = sys.argv[1].rstrip('\n').rstrip(' ')
    my_target = sys.argv[2].rstrip('\n').rstrip(' ')
    
    with open(fasta_seq_in,'rb') as fh:
        myseq_dict = SeqIO.to_dict( SeqIO.parse(fh,'fasta') )

    model_seq = myseq_dict[my_target] #first and only record from the file
    name = model_seq.id
    seq = model_seq.seq

    mypidfh = open(mypidpath+'/%s.pid' %name,'wb')
    mypidfh.write("TargetName  HighestPIDFound \n")  

    print "Name is %s" %name
    print "Seq is %s" %seq

    with open(myseqstructurespath+'/%s.fasta' %name,'wb') as fhseq:
        fhseq.write('>%s\n%s\n' %(name,seq))
        
    os.chdir(homologuesdir)
        
    print "%s" %os.getcwd()
    if not os.path.isdir(name):
        print "making %s" %name
        os.mkdir(name)

    os.chdir(homologuesdir+'/'+name)

    with open('combo.fasta','wb') as fhfasta:
        fhfasta.write('>%s\n%s\n' %(name,seq))

    mytrys = 0
    for db in ['pdb']:
            #create lock file
        if not isfile(mymaindir+'/phmmer.lock'):
            nosuccess = True
            while nosuccess == True:
                try:
                    phmmerfh = open(mymaindir+'/phmmer.lock','wb')
                    phmmerfh.write('%s' %name)
                    phmmerfh.close()
                    print "starting lock file for phmmer %s " %db
                    jsonfile = run_phmmer(model_seq,db)                        
                    print "removing lock file for phmmer"
                    os.remove(mymaindir+'/phmmer.lock')
                    nosuccess = False

                except:
                    print 'ERROR while trying to run %s phmmer on %s' %(name,db)
                    print 'Going to try again'
                    mytrys +=1
                    sleep(randint(1,10))
                    if mytrys >10:
                        print "Exiting, cannot get phmmer to work, tried 10 times :("
                        sys.exit(-1)
                    
        mylisthits = read_json(jsonfile)
        mytemplateseq = []
            
        records = remove_dup_seqs( SeqIO.parse('%s' %db, 'fasta'))
        count = SeqIO.write(records,'%s.clean' %db,'fasta')
        print "Kept %i records" %count

        fh = open('%s.clean' %db,'rb')
        myseqdict = SeqIO.to_dict( SeqIO.parse(fh,'fasta') )
        for i in mylisthits:
            if i in myseqdict:
                mytemplateseq.append( myseqdict[i] )

        cmd = 'cat %s.clean >> combo.fasta' %db
        os.popen(cmd)
        # closes for each database loop

    print "About to align using promals 3d"

    if not isfile('combo.fasta.promals.aln'):
        align_templates_promals3d('combo.fasta')

        #the above line runs promals and maked combo.fasta.promals.aln
    mypromalsaln = 'combo.fasta.promals.aln' #this is in clustal format
    AlignIO.convert('combo.fasta.promals.aln','clustal','combo.msa.fasta','fasta')

#    os.chdir(mymaindir + '/homologoussequences/'+pdbcode)
    with open ('combo.msa.fasta','rb') as seqfh, open('combo.msa.fasta.clean','wb') as myfh:
        alnrec = remove_dup_seqs(  SeqIO.parse(seqfh,'fasta') )
        for record in alnrec:
            try:
                cleanname = record.id.split('|')[1]
                record.id = cleanname
            except:
                print "Error cleaning up names, or maybe they don't need it"
            myfh.write('>%s\n%s\n' %(record.id, record.seq))

    myfh.close()
    print name
    reorder_aln('combo.msa.fasta.clean',name)
    
    alignfh = open('combo.msa.fasta.clean','rb')
    myalignmentdict = SeqIO.to_dict( AlignIO.read(alignfh,'fasta') )
    

        #trim_length
    print "Sorting all pdbs and calculating the PID to target"

    mytemplateseqs = get_pid(name)
    mysortedpdb = sorted([x for x in mytemplateseqs if x.pdb == True],key=lambda x: x.pid,reverse=True)
#        mysortednr = sorted([x for x in mytemplateseqs if x.pdb == False],key=lambda x: x.pid,reverse=True)
    print "My Sorted PDBS"

    if len(mysortedpdb) == 0:
        print "Somethings wrong, I dont have any templates :( "
        sys.exit(-1)

    os.chdir(mymaindir+'/templates')

    if not os.path.isdir(mymaindir+'/templates'+'/'+name):
        os.mkdir(mymaindir+'/templates'+'/'+name)

    mytemplatedir = mymaindir+'/templates'+'/'+name
    print mytemplatedir
    os.chdir(mytemplatedir)

    mytemplates = list()
    max_num_templates_allowed = 10
    mycutoff_pid = 105 #Just has to be larger than 100 do to round off error

    num_templates_below_cutoff = int(len(mysortedpdb))
    print "Using cutoff = %s" %str(mycutoff_pid)
    mycount = 0
    mymax = 0

    stop = False
    while not stop:
        for i in mysortedpdb:
            if float(i.pid) < float(mycutoff_pid) and len(mytemplates) != max_num_templates_allowed:
                if mycount < 1:
                    mymax = float(i.pid)
                    print mymax
                    mypidfh.write("%s %.2f %s \n" %(name,i.pid,i.name))
                    mycount +=1

                tempi = templates.ModelingTemplate(maxpid=mymax, template_dir= mytemplatedir, seq = myalignmentdict[i.name].seq.ungap('.').ungap('-').upper(), name = i.name, alipidrecord = i, target = name, aln_to_target = myalignmentdict[i.name], target_aln_rec = myalignmentdict[name]  )
                mytemplates.append(tempi)
                print "Adding the following template::: %s %s for %s " %(i.name,i.pid,i.target)
                print len(mytemplates)

                if len(mytemplates) == num_templates_below_cutoff:
                    stop = True

            elif len(mytemplates) == max_num_templates_allowed:
                stop = True


    mypartialthread = run.RosettaRun(rosetta_bin='/home/bertolan/rosetta_new/Rosetta/main/source/bin',rosetta_db='/home/bertolan/rosetta_new/Rosetta/main/database',binary='partial_thread.default.linuxgccrelease')
    origseq = mymaindir+'/sequences/'+name+'.fasta'
    mypartialthread.add_flags('-in:file:fasta %s' %origseq)
    mypartialthread.add_flags('-in:file:alignment alignment.grishin')

    try:
        os.remove('alignment.grishin')
    except:
        pass
    mytemp = list()
    for each in mytemplates:
        each.download_pdb()
        each.download_seq()
        each.show_my_aln(make_file=True) #writes each template grishin alignmnet to file
        mytemp.append(each.pdbfile[:-3]+'pdb')

    s = ' '
    spdb = s.join(mytemp)
    mypartialthread.add_flags('-in:file:template_pdb %s' %spdb )
    mypartialthread.add_flags('-ignore_unrecognized_res T')
    mypartialthread.add_flags('-overwrite T')
    print " Running Partial threading of target sequence onto templates found"
    mypartialthread.run()

    os.chdir(mymaindir)
    os.chdir(mymodelspath)
    if not os.path.isdir(name):
        os.mkdir(name)
    os.chdir(name)

    if not os.path.isdir('model'):
        os.mkdir('model')
    os.chdir('model')

       ## now calc the weight for each model, setting the largest equal to 1

    print "Writing qsub engine submission for the hybridize run"
    with open('epiph_%s.sh' %name, 'wb') as fh:
        myhybridizerun = run.RosettaRun(rosetta_bin='/home/bertolan/rosetta_new/Rosetta/main/source/bin',rosetta_db='/home/bertolan/rosetta_new/Rosetta/main/database',binary='rosetta_scripts.default.linuxgccrelease')
        myhybridizerun.add_flags('@flags')
        myhybridizerun.add_flags('-in:file:fasta %s' %origseq)
        myhybridizerun.add_flags('-parser:protocol hybridize.xml')
        myhybridizerun.add_flags('-out:file:silent %s.out' %name)
        myhybridizerun.add_flags('-frag3 %s' %name[0:4]+'_.200.3mers' )
        myhybridizerun.add_flags('-frag9 %s' %name[0:4]+'_.200.9mers')

        fh.write( myhybridizerun.test_run() )


    write_hybridize_flags()
    write_stage1_wts('stage1.wts')
    write_stage2_wts('stage2.wts')
    write_stage3_wts('stage3.wts')
    write_xml(mytemplates,mytemplatedir)

    

    print "About to run James Thompson Constraints on templates"
    mycmd = "%s/cm_scripts/bin/predict_distances_multi.pl %s/alignment.grishin %s/%s  -outfile alignment.grishin.dist_csts -min_seqsep 5 -max_dist 10 -aln_format grishin -max_e_value 10000 -pdb_dir %s " %(environ['ROBETTA'],(mymaindir+'/templates/'+name), (mymaindir+'/sequences'), ( name +'.fasta'),mytemplatedir)
    mytmp = (mymaindir+'/templates/'+ name)
    mycmd2 = 'grep " CA " %s/alignment.grishin.dist_csts.bb_sc > alignment.grishin.dist_csts.bb_sc.CA' %mytmp
    mycmd3 = 'grep " CA " %s/alignment.grishin.dist_csts > alignment.grishin.dist_csts.CA' %mytmp
    print mycmd

    if not isfile('%s/alignment.grishin.dist_csts' %mytemp):
        try:
            os.system(mycmd)
            os.system(mycmd2)
            os.system(mycmd3)
        except OSError as e:
            print "error running predict_distances_multi.pl (James Thompson Evolutionary Constraints) \n %s" %e

    myfragname = name[:4]+'_'
    print "About to make Fragments for the target sequence"
    if (nofrags == False):
        if not isfile( (myfragname+'.make_fragments.success') ):
# Once I get the new one to work, use this
            mycmdfrag = '/share/backup/bertolan/frag_picker_do_not_touch/Rosetta/tools/fragment_tools/make_fragments.pl -verbose -id %s %s ' %(myfragname,(mymaindir+'/sequences/'+name+'.fasta'))
#        mycmdfrag = '/share/tmp-data-1/siegellab/frag_picker_do_not_touch/Rosetta/tools/fragment_tools/make_fragments.pl -nohoms -verbose -id %s %s ' %(myfragname,(mymaindir+'/sequences/'+name+'.fasta'))
            print mycmdfrag
#should add a try 5 times before dying here
            try:
                os.system(mycmdfrag)#? why this and then the next one?
#            subprocess.Popen(mycmdfrag,shell=True)
            except OSError as e:
                print "uh oh! error generating fragments; good luck, maybe you should just have a beer \n %s" %e
    else:
        print "Not making fragments per your request, set nofrags=False to change this"
            #each.update_seq_with_pdb_seq() ## Doesn't appear to be different for my 1 test case, can add in later, just update and then realign
        ## realign new sequences that match PDBS, then write out all the files needed/ convert format to grishin, write xml, rosetta flags, etc


if __name__ == '__main__':
    main()
