from django.template import Context, loader, RequestContext
from chembl_15.models import PfamMaps
from django.core.paginator import Paginator, InvalidPage, EmptyPage
from django.core.urlresolvers import reverse
from django.http import HttpResponse, HttpResponseRedirect
from django.shortcuts import get_object_or_404, render, render_to_response
from django.db import connection
import helper



def index(request):
    t = loader.get_template('chembl_15/index.html')
    c = Context({
        })
    return HttpResponse(t.render(c))


def evidence_portal(request):
    data = helper.custom_sql('SELECT DISTINCT domain_name FROM pfam_maps', [])
    names = sorted([x[0] for x in data])
    t = loader.get_template('chembl_15/evidence_portal_ebi.html')
    c = Context({
        'names': names,
    })
    return HttpResponse(t.render(c))


def evidence(request, pfam_name):
    acts = helper.custom_sql("""
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
    (std_acts, lkp) = helper.standardize_acts(acts)
    (top_mols, top_acts) = helper.filter_acts(std_acts, lkp)
    t = loader.get_template('chembl_15/evidence_ebi.html')
    n_acts = len(std_acts)
    c = Context({
        'top_mols'  : top_mols,
        'top_acts'  : top_acts,
        'pfam_name' : pfam_name,
        'n_acts'    : n_acts
        })
    return HttpResponse(t.render(c))


def conflicts_portal(request):
    act_count = helper.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 0
    """, [])[0][0]
    dub_count = helper.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 1
    """, [])[0][0]
    clash_count = helper.custom_sql("""
    SELECT COUNT(DISTINCT activity_id) FROM pfam_maps WHERE conflict_flag = 2
            """, [])[0][0]
    clash_arch = helper.custom_sql("""
    SELECT DISTINCT activity_id, domain_name FROM pfam_maps WHERE conflict_flag=2""", [])
    clash_arch = helper.process_arch(clash_arch)
    t = loader.get_template('chembl_15/conflict_portal_ebi.html')
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
    data = helper.custom_sql("""
    SELECT DISTINCT compd_id
        FROM pfam_maps
        WHERE activity_id = %s AND domain_name = %s
        """ ,[act, domain_name])[0]
    for compd_id in data:
        new_entry = PfamMaps(activity_id= act, compd_id = compd_id, domain_name = domain_name, conflict_flag = 2, manual_flag = 1 )
        new_entry.save()
    return HttpResponseRedirect(reverse('conflicts', args=(conflict_id,)))


def vote_on_assay(request, conflict_id, assay_id):
    try:
        domain_name = request.POST['choice']
    except KeyError:
        return render_to_response('conflict_ebi.html', {
            'error_message': "You didn't select a choice.",
        }, context_instance=RequestContext(request))
    data = helper.custom_sql("""
    SELECT DISTINCT act.activity_id, compd_id
        FROM pfam_maps pm
        JOIN activities act
          ON act.activity_id = pm.activity_id
        JOIN assays ass
          ON ass.assay_id = act.assay_id
        WHERE ass.chembl_id = %s AND domain_name = %s
        """ ,[assay_id, domain_name])
    for ent in data:
        act = ent[0]
        compd_id = ent[1]
        new_entry = PfamMaps(activity_id= act, compd_id = compd_id, domain_name = domain_name, conflict_flag = 2, manual_flag = 1 )
        new_entry.save()
    return HttpResponseRedirect(reverse('conflicts', args=(conflict_id,)))


def conflicts(request, conflict_id):
    data = helper.custom_sql("""
    SELECT DISTINCT pm.domain_name, ass.chembl_id
        FROM pfam_maps pm
        JOIN activities act
          ON act.activity_id = pm.activity_id
        JOIN assays ass
          ON act.assay_id = ass.assay_id
        WHERE pm.conflict_flag=2
        """, [])
    clash_arch = helper.arch_assays(data)
    arch_assays = clash_arch[conflict_id]
    assay_hier = {}
    for ass_id in arch_assays:
        assay_hier[ass_id] = {}
    paginator = Paginator(assay_hier.keys(), 1)
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1
    try:
        arch_idx = paginator.page(page)
    except (EmptyPage, InvalidPage):
        arch_idx = paginator.page(paginator.num_pages)
    assay_hier_page = dict((k, assay_hier[k]) for k in arch_idx.object_list)
    assay_hier_page = helper.get_assay_meta(assay_hier_page)
    assay_hier_page = helper.get_pfam_arch(assay_hier_page)
    dom_l = conflict_id.split(' vs. ')
    c = {'arch'         : conflict_id,
         'assay_hier'   : assay_hier_page,
         'doms'         : dom_l,
         'arch_idx'     : arch_idx
        }
    return render_to_response('chembl_15/conflict_ebi.html',c, context_instance=RequestContext(request))



# in progress
def resolved(request, conflict_id):
    data = helper.custom_sql("""
    SELECT DISTINCT pm.domain_name, ass.chembl_id
        FROM pfam_maps pm
        JOIN activities act
          ON act.activity_id = pm.activity_id
        JOIN assays ass
          ON act.assay_id = ass.assay_id
        WHERE pm.manual_flag=1
        """, [])
    clash_arch = helper.arch_assays(data)
    arch_assays = clash_arch[conflict_id]
    assay_hier = {}
    for ass_id in arch_assays:
        assay_hier[ass_id] = {}
    paginator = Paginator(assay_hier.keys(), 1)
    try:
        page = int(request.GET.get('page', '1'))
    except ValueError:
        page = 1
    try:
        arch_idx = paginator.page(page)
    except (EmptyPage, InvalidPage):
        arch_idx = paginator.page(paginator.num_pages)
    assay_hier_page = dict((k, assay_hier[k]) for k in arch_idx.object_list)
    assay_hier_page = helper.get_assay_meta(assay_hier_page)
    assay_hier_page = helper.get_pfam_arch(assay_hier_page)
    dom_l = conflict_id.split(' vs. ')
    c = {'arch'         : conflict_id,
         'assay_hier'   : assay_hier_page,
         'doms'         : dom_l,
         'arch_idx'     : arch_idx
        }
    return render_to_response('chembl_15/conflict_ebi.html',c, context_instance=RequestContext(request))



def details(request, act):
    data = helper.custom_sql("""
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
    pfam_arch = helper.get_pfam_arch(data)
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
         'pfam_arch'     : pfam_arch
        }
    return render_to_response('chembl_15/details_ebi.html',c,                          context_instance=RequestContext(request))
