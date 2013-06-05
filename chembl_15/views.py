from django.template import Context, loader
from chembl_15.models import PfamMaps
from django.http import HttpResponse
from django.db import connection
import queryDevice


def process_sql(acts):
    import numpy as np
    import simplejson as json
    import yaml
    paramFile = open('local.yaml')
    params = yaml.safe_load(paramFile)
    ki_adjust = params['ki_adjust']
    lkp  = {}
    for data in acts:
        try:
            standard_value = float(data[0])
        except ValueError:
            continue
        standard_type = data[1]
        standard_units = data[2]
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
        try:
            lkp[molregno].append((standard_value, accession, act_id))
        except KeyError:
            lkp[molregno] = [(standard_value, accession, act_id)]
    lkp = json.dumps(lkp)
    return(lkp)


def index(request):
    data = queryDevice.custom_sql('SELECT DISTINCT domain_name FROM pfam_maps', [])
    names = sorted([x[0] for x in data])
    t = loader.get_template('chembl_15/index.html')
    c = Context({
        'names': names,
    })
    return HttpResponse(t.render(c))


def evidence(request, pfam_name):
    acts = queryDevice.custom_sql('SELECT DISTINCT act.standard_value, act.standard_type, act.standard_units, act.activity_id, act.molregno, cs.accession FROM pfam_maps pm JOIN activities act ON act.activity_id = pm.activity_id JOIN component_domains cd ON cd.compd_id = pm.compd_id JOIN component_sequences cs ON cd.component_id = cs.component_id WHERE domain_name = %s',[pfam_name])

    query_out = process_sql(acts)
    t = loader.get_template('chembl_15/evidence.html')
    c = Context({
        'acts' : query_out,
        'query_out' : query_out,
    })
    return HttpResponse(t.render(c))
