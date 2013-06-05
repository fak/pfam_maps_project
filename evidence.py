"""
Function:  evidence

Retrieves evidence for small molecule domain interactions from the database and feeds them forward.
    --------------------
    Author:
    Felix Kruger
    fkrueger@ebi.ac.uk
"""
import numpy as np
import os
impor yaml


def master(release, version):


    ## Get all human protein coding genes from ensembl and a count of all uniqe domains.
    # Read config file.
    #paramFile = open('mpf.yaml')
    paramFile = open('mpf_local.yaml')
    params = yaml.safe_load(paramFile)
    user = params['user']
    pword = params['pword']
    host = params['host']
    port = params['port']
    th = params['threshold']



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
