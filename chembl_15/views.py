from django.template import Context, loader, RequestContext
from chembl_15.models import PfamMaps
from django.core.urlresolvers import reverse
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, render_to_response
from django.db import connection
import queryDevice
import numpy as np
import yaml
import simplejson as json
from operator import itemgetter
import requests

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

def process_arch(data):
    lkp = {}
    for clash in data:
        act_id = clash[0]
        dom = str(clash[1])
        try:
            lkp[act_id] = lkp[act_id] +'  vs. ' + dom
        except KeyError:
            lkp[act_id] = dom
    inv_lkp = dictinvert(lkp)
    return inv_lkp

def dictinvert(d):
        inv = {}
        for k, v in d.iteritems():
            try:
                inv[v].append(k)
            except KeyError:
                inv[v] = [k]
        return inv

def process_acts(arch_acts, data):
    arch_mols = {}
    for clash in data:
        act_id = clash[0]
        if act_id in arch_acts:
             molregno = clash[2]
             try:
                 arch_mols[molregno][act_id]=0
             except KeyError:
                 arch_mols[molregno] = {act_id:0}
    return arch_mols


##############################
### Definition of views ######
##############################

def index(request):
    t = loader.get_template('chembl_15/index.html')
    c = Context({
        })
    return HttpResponse(t.render(c))


def evidence_portal(request):
    data = queryDevice.custom_sql('SELECT DISTINCT domain_name FROM pfam_maps', [])
    names = sorted([x[0] for x in data])
    t = loader.get_template('chembl_15/evidence_portal.html')
    c = Context({
        'names': names,
    })
    return HttpResponse(t.render(c))


def evidence(request, pfam_name):
    acts = queryDevice.custom_sql("""
    SELECT DISTINCT act.standard_value, act.standard_units, act.standard_type, act.  activity_id, act.molregno, single_domains.accession
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
        WHERE domain_name = %s AND standard_relation= '=' AND assay_type = 'B' AND relationship_type = 'D' LIMIT 1500
        """ , [pfam_name])

    (std_acts, lkp) = standardize_acts(acts)
    (top_mols, top_acts) = filter_acts(std_acts, lkp)
    t = loader.get_template('chembl_15/evidence.html')
    n_acts = len(std_acts)
    c = Context({
        'top_mols'  : top_mols,
        'top_acts'  : top_acts,
        'pfam_name' : pfam_name,
        'n_acts'    : n_acts
        })
    return HttpResponse(t.render(c))


def conflicts_portal(request):
    act_count = queryDevice.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 0
    """, [])[0][0]
    dub_count = queryDevice.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 1
    """, [])[0][0]
    clash_count = queryDevice.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 2
            """, [])[0][0]
    clash_arch = queryDevice.custom_sql("""
    SELECT DISTINCT activity_id, domain_name FROM pfam_maps WHERE conflict_flag=2""", [])
    clash_arch = process_arch(clash_arch)
    t = loader.get_template('chembl_15/conflicts_portal.html')
    c = Context({
        'act_count'   : act_count,
        'dub_count'   : dub_count,
        'clash_count' : clash_count,
        'clash_arch'  : clash_arch.keys(),
        })
    return HttpResponse(t.render(c))

def vote(request, conflict_id, act):
    try:
        domain_name = request.POST['choice']
    except KeyError:
        return render_to_response('index.html', {
            'error_message': "You didn't select a choice.",
        }, context_instance=RequestContext(request))
    data = queryDevice.custom_sql("""
    SELECT DISTINCT compd_id
        FROM pfam_maps
        WHERE activity_id = %s AND domain_name = %s
        """ ,[act, domain_name])[0]
    for compd_id in data:
        new_entry = PfamMaps(activity_id= act, compd_id = compd_id, domain_name = domain_name, conflict_flag = 2, manual_flag = 1 )
        new_entry.save()
    return HttpResponseRedirect(reverse('conflicts', args=(conflict_id,)))

from django.core.paginator import Paginator, InvalidPage, EmptyPage
def conflicts(request, conflict_id):

    data = queryDevice.custom_sql("""
    SELECT DISTINCT pm.activity_id, pm.domain_name, act.molregno
        FROM pfam_maps pm
        JOIN activities act
          ON act.activity_id = pm.activity_id
        WHERE pm.conflict_flag=2
        """, [])
    clash_arch = process_arch(data)
    arch_acts = clash_arch[conflict_id]
    arch_mols = process_acts(arch_acts, data)
    paginator = Paginator(arch_mols.keys(), 25)
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1
    try:
        arch_idx = paginator.page(page)
    except (EmptyPage, InvalidPage):
        arch_idx = paginator.page(paginator.num_pages)
    print arch_idx.object_list
    print arch_idx.number
    domL = conflict_id.split(' vs. ') # check for defaultdict syntax, this should make this code cleaner.
    domstr = "','".join(domL)
    counts = queryDevice.custom_sql("""
    SELECT activity_id, domain_name
        FROM pfam_maps
        WHERE manual_flag = 1 AND domain_name IN('%s')
        """%domstr, [])
    doms = {}
    for act in arch_acts:
        doms[act]={}
        for dom in domL:
            doms[act][dom]=0
    for ent in counts:
        dom = ent[1]
        act = ent[0]
        try:
            doms[act][dom] +=1
        except KeyError:
            pass
    arch_mols_page = dict((k, arch_mols[k]) for k in arch_idx.object_list)
    c = {'arch'     : conflict_id,
         'arch_mols': arch_mols_page,
         'doms'     : doms,
         'arch_idx' : arch_idx
        }
    return render_to_response('chembl_15/conflict.html',c, context_instance=RequestContext(request))



def get_arch(data):
    trc = []
    #[(domain_name, start, end), ...]
    for ent in data:

        domain_name = ent[6]
        start = ent[4]
        end = ent[5]
        trc.append((domain_name, start, end))
    return trc

def details(request, act):
    data = queryDevice.custom_sql("""
    SELECT DISTINCT dcs.chembl_id, dcs.pubmed_id, act.molregno, cs.accession, td.pref_name, dcs.title, ass.chembl_id, ass.description
        FROM activities act
        JOIN docs dcs
          ON act.doc_id = dcs.doc_id
        JOIN assays ass
          ON ass.assay_id = act.assay_id
        JOIN target_dictionary td
          ON td.tid = ass.tid
        JOIN target_components tc
          ON tc.tid = td.tid
        JOIN component_sequences cs
          ON tc.component_id = cs.component_id
        WHERE activity_id = %s
        """, [act])

    doc_id =  data[0][0]
    pubmed_id =  data[0][1]
    molregno = data[0][2]
    accession = data[0][3]
    pref_name = data[0][4]
    title = data[0][5]
    ass_id = data[0][6]
    desc = data[0][7]
    dom_arch = {}
    for ent in data:
        uniprot = ent[3]
        try:
            dom_arch[uniprot]
        except KeyError:
            #doms = '[{"length":"950","regions":[{"colour":"#2dcf00","endStyle":"jagged","end":"361","startStyle":"jagged","text":"Peptidase_S8","href":"/family/PF00082","type":"pfama","start":"159"},]}]'
            r = requests.get('http://pfam.sanger.ac.uk/protein/%s/graphic' % uniprot)
            doms = r.content
            dom_arch[uniprot] = doms

    if not pubmed_id:
        pubmed_id = 'NA'
    if not title:
        title = 'NA'
    if not pref_name:
        pref_name = 'NA'
    if not desc:
        desc = 'NA'
    c = {'doc_id'       : doc_id,
         'pubmed_id'    : pubmed_id,
         'molregno'     : molregno,
         'title'        : title,
         'pref_name'    : pref_name,
         'accession'    : accession,
         'ass_id'       : ass_id,
         'desc'         : desc,
         'act'          : act,
         'dom_arch'     : dom_arch
        }
    return render_to_response('chembl_15/details.html',c,                          context_instance=RequestContext(request))
