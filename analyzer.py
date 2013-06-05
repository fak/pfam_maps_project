"""
Function:  master

Goes through all necessary steps.
    --------------------
    Author:
    Felix Kruger
    fkrueger@ebi.ac.uk
"""

def readfile(path, keykey, valkey):
    infile = open(path, 'r')
    lines = infile.readlines()
    infile.close()
    lkp = {}
    els = lines[0].rstrip().split('\t')
    print els
    for i, el in enumerate(els):
        if el == keykey:
            idx = i
        if el == valkey:
            vdx = i
    for line in lines[1:]:
        elements = line.rstrip().split('\t')
        lkp[elements[idx]] = elements[vdx]
    return  lkp



def master(release, version):


    ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
    import os
    import yaml
    # Read config file.
    #paramFile = open('mpf.yaml')
    paramFile = open('mpf_local.yaml')
    params = yaml.safe_load(paramFile)
    user = params['user']
    pword = params['pword']
    host = params['host']
    port = params['port']
    th = params['threshold']



    # Load the list of domains.
    domains = readfile('data/valid_pfam_v_%(version)s.tab' % locals(), 'pfam_a', 'pfam_a')
    dom_string = "','".join(domains.keys())

    # Load the list of resolved conflicts.
    conflicts = readfile('data/resolved_conflicts_v_%(version)s.tab' % locals(), 'activity_id', 'pfam_id')


    # Load activities from chembl and keep those that are active -
    # for this need standard_type, standard_value, stadard_units, standard_relation
    acts = queryDevice.queryDevice(
         """SELECT DISTINCT standard_value, standard_type, standard_units, act.standard_relation, act.activity_id
                      FROM activities act
                      JOIN assays ass
                          ON ass.assay_id = act.assay_id
                      JOIN target_dictionary td
                          ON ass.tid = td.tid
                     WHERE ass.assay_type IN('B')
                     AND td.target_type IN('PROTEIN COMPLEX', 'SINGLE PROTEIN')
                     AND act.standard_relation ='='
                     AND ass.relationship_type = 'D'
                     AND act.standard_type IN(
                       'Ki', 'Kd', 'IC50', 'EC50','AC50'
                       'log Ki', 'log Kd', 'log IC50', 'Log EC50', 'Log AC50'
                       'pKi', 'pKd', 'pIC50', 'pEC50', 'pAC50')"""
                            ,release, user, pword, host, port)

    import numpy as np
    for data in acts:
        try:
            standard_value = float(data[0])
        except ValueError:
            continue
        standard_type = data[1]
        standard_units = data[2]
        act_id = data[0]
        # p-scaling.
        if standard_type in ['Ki','Kd','IC50','EC50', 'AC50'] and standard_units == 'nM':
            standard_value = -(np.log10(standard_value)-9)
            standard_type = 'p' + standard_type
        # p-scaling.
        if standard_type in ['log Ki', 'log Kd', 'log IC50', 'log EC50', 'logAC50'] and standard_units is None:
            standard_value = - standard_value
            standard_type = 'p' + standard_type.split(' ')[1]
        # Mixing types.
        if standard_type in ['pKi', 'pKd']:
            standard_value = standard_value - ki_adjust
        if standard_value <= threshold:
            #write out value, take molregno.


    # construct lookup dictionary - need a query that joins
    # act_id, tid, component_id, domain_id and domain_name.
    acts = queryDevice.queryDevice(
         """SELECT DISTINCT act.activity_id, ass.tid, tc.component_id, cd.compd_id, dm.domain_name
                      FROM activities act
                      JOIN assays ass
                          ON ass.assay_id = act.assay_id
                      JOIN target_dictionary td
                          ON ass.tid = td.tid
                      JOIN target_components tc
                          ON ass.tid = tc.tid
                      JOIN component_domains cd
                          ON tc.component_id = cd.component_id
                      JOIN domains dm
                          ON dm.domain_id = cd.domain_id
                     WHERE ass.assay_type IN('B','F')
                     AND td.target_type IN('PROTEIN COMPLEX', 'SINGLE PROTEIN')
                     AND act.standard_relation ='='
                     AND ass.relationship_type = 'D'
                     AND act.standard_type IN(
                       'Ki', 'Kd', 'IC50', 'EC50', 'AC50',
                       'log Ki', 'log Kd', 'log IC50', 'Log EC50', 'Log AC50'
                       'pKi', 'pKd', 'pIC50', 'pEC50', 'pAC50'
                        )
                     AND dm.domain_name IN('%(dom_string)s')"""
                           %locals() ,release, user, pword, host, port)

    # Map interactions to activity ids.
    lkp = {}
    for act in acts:
        (act_id, tid, component_id, compd_id, domain_name) = act
        try:
            lkp[act_id][compd_id]=domain_name
        except KeyError:
            lkp[act_id] ={}
            lkp[act_id][compd_id]=domain_name

    # Flag conflicts.
    flag_lkp  = {}
    for act_id in lkp.keys():
        if len(lkp[act_id].keys()) == 1:
            flag_lkp[act_id]=(0,0)
        elif len(lkp[act_id].keys()) > 1:
            if len(set(lkp[act_id].values())) > 1:
                flag_lkp[act_id] = (2,0)
            elif len(set(lkp[act_id].values())) == 1:
                flag_lkp[act_id] = (1,0)

    # Apply manually resolved conflicts.
    infile = 'data/resolved_conflicts_v_%(version)s.tab' % locals()
    infile = open(infile, 'r')
    lines = infile.readlines()
    infile.close()
    for line in lines[1:]:
        (act_id, compd_id, domain_name, comment) = line.rstrip().split('\t')
        try:
            lkp[act_id][compd_id][domain_name]
            flag_lkp[act_id] = (0,1)
        except KeyError:
            print 'lost entry:', line


    # write a table containing activity_id, domain_id, tid, conflict_flag, type_flag
    out = open('data/pfam_maps_v_%(version)s.tab' %locals(), 'w')
    out.write('activity_id\tcompd_id\tdomain_name\tconflict_flag\tmanual_flag\n')
    for act_id in lkp.keys():
        compd_ids = lkp[act_id]
        (conflict_flag, manual_flag) = flag_lkp[act_id]
        for compd_id in compd_ids.keys():
            domain_name = lkp[act_id][compd_id]
            out.write('%(act_id)s\t%(compd_id)s\t%(domain_name)s\t%(conflict_flag)s\t%(manual_flag)s\n'%locals())
    out.close()

    os.system("mysql -u%(user)s -p%(pword)s -h%(host)s -P%(port)s -e 'CREATE TABLE pfam_maps(activity_id  INT, compd_id INT , domain_name VARCHAR(100) , conflict_flag INT, manual_flag INT)' %(release)s"% locals())

    os.system("cp data/pfam_maps_v_%(version)s.tab data/pfam_maps.txt" % locals())

    os.system("mysqlimport -u%(user)s -p%(pword)s -h%(host)s -P%(port)s --lines-terminated-by='\n' --local %(release)s data/pfam_maps.txt"% locals())



    # Analyze the data.
    #analysis.analysis(th, release, user, pword, host, port)

if __name__ == '__main__':
    import sys

    if len(sys.argv) < 5:  # the program name and the two arguments

        sys.exit("Must specify release, user, pword, host, port")


    release = sys.argv[1]
    user = sys.argv[2]
    pword = sys.argv[3]
    host = sys.argv[4]
    port = int(sys.argv[5])

    master(release, user, pword, host, port)
