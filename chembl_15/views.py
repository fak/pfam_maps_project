from django.template import Context, loader
from chembl_15.models import PfamMaps
from django.http import HttpResponse
from django.db import connection
import queryDevice
import numpy as np
import yaml
import simplejson as json
from operator import itemgetter


def standardize_acts(acts):
    paramFile = open('local.yaml')
    params = yaml.safe_load(paramFile)
    ki_adjust = params['ki_adjust']
    std_acts = []
    lkp = {}
    for data in acts:
        try:
            standard_value = float(data[0])
        except TypeError:
            continue
        standard_units = data[1]
        standard_type = data[2]
        act_id = data[3]
        molregno = data[4]
        accession = data[5]

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
        # Filtering inactives.
        if standard_value >= 3:
            std_acts.append((molregno, standard_value, accession, act_id))
            try:
                lkp[molregno] += 1
            except KeyError:
                lkp[molregno] = 1
    return (std_acts, lkp)

def add_meta(top_acts):
    act_ids = [str(x[3]) for x in top_acts]
    act_str = "','".join(act_ids)
    data = queryDevice.custom_sql("""SELECT act.activity_id, td.pref_name, dcs.title
                                            FROM activities act
                                            JOIN assays ass
                                                ON act.assay_id = ass.assay_id
                                            JOIN target_dictionary td
                                                ON ass.tid = td.tid
                                            JOIN docs dcs
                                                ON act.doc_id = dcs.doc_id
                                            WHERE activity_id IN('%s')""" % act_str, [])
    lkp = {}
    for meta in data:
        act = meta[0]
        lkp[act] = (meta[1], meta[2])
    top_acts = [x + lkp[x[3]]  for x in top_acts]
    return top_acts

def filter_acts(std_acts, lkp):
    top_mols = [key for key,value in sorted(lkp.items(), key=itemgetter(1), reverse = True)][:10]
    top_acts = [x for x in std_acts if x[0] in top_mols]
    top_acts = add_meta(top_acts)
    top_mols = json.dumps(top_mols)
    top_acts = json.dumps(top_acts)
    return(top_mols, top_acts)


def index(request):
    data = queryDevice.custom_sql('SELECT DISTINCT domain_name FROM pfam_maps', [])
    names = sorted([x[0] for x in data])
    t = loader.get_template('chembl_15/index.html')
    c = Context({
        'names': names,
    })
    return HttpResponse(t.render(c))


def evidence(request, pfam_name):
    acts = queryDevice.custom_sql(
"""SELECT DISTINCT act.standard_value, act.standard_units, act.standard_type, act.  activity_id, act.molregno, single_domains.accession
    FROM pfam_maps pm
        JOIN activities act
            ON act.activity_id = pm.activity_id
        JOIN assays ass
            ON act.assay_id = ass.assay_id
        JOIN (SELECT  tid, cs.accession
                FROM component_domains cd
                    JOIN component_sequences cs
                        ON cd.component_id = cs.component_id
                    JOIN target_components tc
                        ON tc.component_id = cs.component_id
                    GROUP BY tid
                    HAVING COUNT(compd_id) =1)
        AS single_domains
        ON single_domains.tid = ass.tid
        WHERE domain_name = %s AND standard_relation= '=' AND assay_type = 'B' AND relationship_type = 'D' LIMIT 1500""" , [pfam_name])

    (std_acts, lkp) = standardize_acts(acts)
    (top_mols, top_acts) = filter_acts(std_acts, lkp)
    t = loader.get_template('chembl_15/evidence.html')
    c = Context({
        'top_mols' : top_mols,
        'top_acts' : top_acts,
        'pfam_name': pfam_name,
        })
    return HttpResponse(t.render(c))
