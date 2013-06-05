"""
    Function:  addMolweight

    adds the column Molweight to paralogFull.tab and writes out to 
    paralogMolweight.tab
    --------------------

    Felix Kruger
    momo.sander@googlemail.com
"""

def addMolweight(path, vsChembl, user, pword, host, port):
    import queryDevice
    import os
    import time    
                                                      
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()

    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tmolweight\tlogP\n'%lines[0].rstrip('\n')) 

    els = lines[0].split('\t')
    for i, el in enumerate(els):
        if el == 'molregno':
            idx = i
            break    
    for line in lines[1:]:
        time.sleep(0.03)
        elements = line.split('\t')
        molregno = elements[idx]
        try:
            mw = queryDevice.queryDevice("SELECT mw_freebase FROM compound_properties WHERE molregno = %s "% molregno, vsChembl, user, pword, host, port)[0][0]
        except IndexError:
            mw = None
        try:
            logP = queryDevice.queryDevice("SELECT alogp FROM compound_properties WHERE molregno = %s" % molregno, vsChembl, user, pword, host, port)[0][0]
        except IndexError:
            mw = None
        out.write("%s\t%s\t%s\n"%(line.rstrip('\n'), mw, logP  ))                                                                      
    out.close()
    os.system('mv  %s %s'% ('_'.join([path,"sed"]), path))  


def addTargetClass(level, key, path, vsChembl, user, pword, host, port):
    import queryDevice
    import os
    import time

    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()

    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\ttargetClass_%s\n'%(lines[0].rstrip('\n'),level) )
    header = lines[0].split('\t')
    for i, col in enumerate(header):
        if col == key:
            idx = i
            break
    for line in lines[1:]:
        time.sleep(0.03)     
        elements = line.split('\t')
        uniprot = elements[idx]
        targetClass = queryDevice.queryDevice("SELECT %s FROM target_class tc JOIN target_dictionary td ON td.tid = tc.tid  WHERE protein_accession = '%s' "%(level, uniprot),vsChembl, user, pword, host, port)
        try:
            targetClass = targetClass[0][0]
        except IndexError:
            targetClass = None
        out.write("%s\t%s\n"%(line.rstrip('\n'), targetClass ))
    out.close()
    os.system('mv %s %s'% ('_'.join([path,"sed"]), path))



def addSeq100(path):
    import os

    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()

    out = open('_'.join([path,"sed"]) ,'w')
    out.write('%s\tseqId\n'%lines[0].rstrip('\n'))
    for line in lines[1:]:   
        out.write("%s\t%s\n"%(line.rstrip('\n'), 100 ))    
    out.close()
    os.system('mv  %s %s'% ('_'.join([path,"sed"]), path))






