from flask import Flask, render_template, json, request, url_for
import requests, sys, time, threading, re, math, random, cgi
from jinja2 import Template
from gp_clear_code import get_transcripts, get_exons, get_primers, get_seq_len, get_reverse_complement_seq, get_gene_id, get_gc_content_f, primers_design_time_counter, pb_server_status_checker, get_clear_seq_for_textarea

application = Flask(__name__)


@application.route('/reverse_complement', methods=['GET', 'POST'])
def reverse_complement():
    return render_template('rev_comp.html')


@application.route('/', methods=['GET', 'POST'])
@application.route('/GRCh37', methods=['GET', 'POST'])
def GRCh37():
    gp_request_id = random.randint(1, 2048)
    return render_template('index.html', gp_request_id=gp_request_id)


@application.route('/GRCh38', methods=['GET', 'POST'])
def GRCh38():
    gp_request_id = random.randint(1, 2048)
    return render_template('GRCh38.html', gp_request_id=gp_request_id)


@application.route('/about')
def about():
    return render_template('about.html')


@application.route('/get_clear_seq', methods=['GET', 'POST'])
def get_clear_seq():
    seq_for_rc = request.form['seq_for_rc'];
    seq = get_clear_seq_for_textarea(seq_for_rc)
    print(seq)
    return seq


@application.route('/get_len_seq', methods=['GET', 'POST'])
def get_len_seq():
    seq_for_rc = request.form['seq_for_rc'];
    return render_template('get_len_seq.html', len_seq=get_seq_len(seq_for_rc))


@application.route('/get_reverse_complement', methods=['GET', 'POST'])
def get_reverse_complement():
    seq_for_rc = request.form['seq_for_rc'];
    rev_comp_seq = get_reverse_complement_seq(seq_for_rc)
    return rev_comp_seq


@application.route('/get_len_seq_2', methods=['GET', 'POST'])
def get_len_seq_2():
    seq_for_rc_2 = request.form['seq_for_rc_2'];
    return render_template('get_len_seq.html', len_seq=get_seq_len(seq_for_rc_2))


@application.route('/get_reverse_complement_2', methods=['GET', 'POST'])
def get_reverse_complement_2():
    seq_for_rc_2 = request.form['seq_for_rc_2'];
    return render_template('reverse_complement.html', reverse_complement=get_reverse_complement_seq(seq_for_rc_2))


@application.route('/get_gc_content_page', methods=['GET', 'POST'])
def get_gc_content_fa():
    seq_for_gc_cont = request.form['seq_for_rc'];
    return render_template('gc_cont.html', gc_content=get_gc_content_f(seq_for_gc_cont))


@application.route('/get_gc_content_page_2', methods=['GET', 'POST'])
def get_gc_content_fa_2():
    seq_for_gc_cont = request.form['seq_for_rc_2'];
    return render_template('gc_cont.html', gc_content=get_gc_content_f(seq_for_gc_cont))


@application.route('/get_pb_server_status', methods=['GET', 'POST'])
def get_pb_server_status():
    print('PB_Status')
    return render_template('pb_server_status.html', pb_server_status=pb_server_status_checker())


@application.route('/get_transcripts', methods=['GET', 'POST'])
def get_transcripts_app():
    #global file_name_pb
    gp_request_id = request.form['gp_request_id'];
    file_name_pb = 'pb_file_' + str(gp_request_id)
    print(file_name_pb)
    with open('params_dir/{}.json'.format(file_name_pb), 'w') as params:
        params.write('')
    SPECIES = request.form['SPECIES'];
    gene = (request.form['gene']).replace(' ', '');
    if gene == '':
        return render_template('wrong_gene.html', gene=gene.upper())
    GCS = request.form['GCS'];
    print('\nGenomic coordinate system - {}'.format(GCS))
    gene_id = get_gene_id(SPECIES, gene, GCS)
    if gene_id != 'Error. Wrong gene.':
        gene_id = gene_id['gene_id']
    if gene_id == 'Error. Wrong gene.':
        return render_template('wrong_gene.html', gene=gene.upper())
    data = get_transcripts(SPECIES, gene, GCS)
    if data == 'Error. Wrong gene.':
        return render_template('wrong_gene.html', gene=gene.upper())
    else:
        return render_template('transcripts.html', data=data, gene_id=gene_id, GCS_data=GCS, SPECIES=SPECIES)


@application.route('/get_exons', methods=['GET', 'POST'])
def get_exons_app():
    gp_request_id = request.form['gp_request_id'];
    file_name_pb = 'pb_file_' + str(gp_request_id)
    transcript = request.form['transcript'];
    print(transcript)
    transcript = re.sub(r' - \w*.\w*', '', transcript)
    print(transcript)
    GCS = request.form['GCS'];
    if transcript == 'Select transcript':
        return render_template('pls_select_transcript.html')
    else:
        return render_template('exons.html', data=get_exons(file_name_pb, transcript, GCS)['exons_count'])


@application.route('/results', methods=['GET', 'POST'])
def results():
    form_params_dict = {}

    gp_request_id = request.form['gp_request_id'];
    print(gp_request_id)

    SPECIES = request.form['SPECIES'];
    form_params_dict['SPECIES']=SPECIES

    gene = (request.form['gene']).replace(' ', '');
    form_params_dict['gene']=gene

    GCS = request.form['GCS'];
    form_params_dict['GCS']=GCS

    transcript = request.form['transcript'];
    form_params_dict['transcript']=transcript

    taken_exons = request.form.getlist('exon');
    form_params_dict['taken_exons']=taken_exons

    SEARCH_SPECIFIC_PRIMER = request.form.getlist('SEARCH_SPECIFIC_PRIMER');
    form_params_dict['SEARCH_SPECIFIC_PRIMER']=SEARCH_SPECIFIC_PRIMER

    auto_distances = request.form.getlist('auto_distances');
    form_params_dict['auto_distances']=auto_distances

    CROSS_SEARCH = request.form.getlist('CROSS_SEARCH');
    form_params_dict['CROSS_SEARCH']=CROSS_SEARCH

    NO_SNP = request.form.getlist('NO_SNP');
    form_params_dict['NO_SNP']=NO_SNP

    SHOW_PB_LINK = request.form.getlist('SHOW_PB_LINK');
    form_params_dict['SHOW_PB_LINK']=SHOW_PB_LINK

    PRIMER_MIN_TM_PB = request.form['PRIMER_MIN_TM'];
    form_params_dict['PRIMER_MIN_TM_PB']=PRIMER_MIN_TM_PB

    PRIMER_OPT_TM_PB = request.form['PRIMER_OPT_TM'];
    form_params_dict['PRIMER_OPT_TM_PB']=PRIMER_OPT_TM_PB

    PRIMER_MAX_TM_PB = request.form['PRIMER_MAX_TM'];
    form_params_dict['PRIMER_MAX_TM_PB']=PRIMER_MAX_TM_PB

    PRIMER_MAX_DIFF_TM_PB = request.form['PRIMER_MAX_DIFF_TM'];
    form_params_dict['PRIMER_MAX_DIFF_TM_PB']=PRIMER_MAX_DIFF_TM_PB

    PRIMER_MIN_SIZE_PB = request.form['PRIMER_MIN_SIZE'];
    form_params_dict['PRIMER_MIN_SIZE_PB']=PRIMER_MIN_SIZE_PB

    PRIMER_OPT_SIZE_PB = request.form['PRIMER_OPT_SIZE'];
    form_params_dict['PRIMER_OPT_SIZE_PB']=PRIMER_OPT_SIZE_PB

    PRIMER_MAX_SIZE_PB = request.form['PRIMER_MAX_SIZE'];
    form_params_dict['PRIMER_MAX_SIZE_PB']=PRIMER_MAX_SIZE_PB

    CROSS_EXONS_MAX_SIZE_PB = request.form['CROSS_EXONS_MAX_SIZE'];
    form_params_dict['CROSS_EXONS_MAX_SIZE_PB']=CROSS_EXONS_MAX_SIZE_PB

    POLYX_PB = request.form['POLYX'];
    form_params_dict['POLYX_PB']=POLYX_PB

    PRIMER_PRODUCT_MIN_PB = request.form['PRIMER_PRODUCT_MIN'];
    form_params_dict['PRIMER_PRODUCT_MIN_PB']=PRIMER_PRODUCT_MIN_PB

    PRIMER_PRODUCT_MAX_PB = request.form['PRIMER_PRODUCT_MAX'];
    form_params_dict['PRIMER_PRODUCT_MAX_PB']=PRIMER_PRODUCT_MAX_PB

    PRIMER_MIN_GC_PB = request.form['PRIMER_MIN_GC'];
    form_params_dict['PRIMER_MIN_GC_PB']=PRIMER_MIN_GC_PB

    PRIMER_MAX_GC_PB = request.form['PRIMER_MAX_GC'];
    form_params_dict['PRIMER_MAX_GC_PB']=PRIMER_MAX_GC_PB

    FIVE_SAVE_EXON_DISTANCE_PB = request.form['FIVE_SAVE_EXON_DISTANCE'];
    form_params_dict['FIVE_SAVE_EXON_DISTANCE_PB']=FIVE_SAVE_EXON_DISTANCE_PB

    THREE_SAVE_EXON_DISTANCE_PB = request.form['THREE_SAVE_EXON_DISTANCE'];
    form_params_dict['THREE_SAVE_EXON_DISTANCE_PB']=THREE_SAVE_EXON_DISTANCE_PB

    F_SEARCH_DISTANCE_PB = request.form['F_SEARCH_DISTANCE'];
    form_params_dict['F_SEARCH_DISTANCE_PB']=F_SEARCH_DISTANCE_PB

    R_SEARCH_DISTANCE_PB = request.form['R_SEARCH_DISTANCE'];
    form_params_dict['R_SEARCH_DISTANCE_PB']=R_SEARCH_DISTANCE_PB

    MAX_MAF = request.form['MAX_MAF'];
    form_params_dict['MAX_MAF']=MAX_MAF

    pb_server_status=pb_server_status_checker()

    #global form_params_file_name_pb
    form_params_file_name_pb = 'form_params_pb_file_' + str(gp_request_id)
    print(form_params_file_name_pb)
    with open('form_params_dir/{}.json'.format(form_params_file_name_pb), 'w') as form_params:
        json.dump(form_params_dict, form_params)

    statuses_dict = {}
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
        json.dump(statuses_dict, statuses_file)

    return render_template('loading.html', gp_request_id=gp_request_id, primers_design_time_counter=primers_design_time_counter(taken_exons, SEARCH_SPECIFIC_PRIMER, NO_SNP, pb_server_status), pb_server_status=pb_server_status, gene=gene.upper())


@application.route('/ajax_results', methods=['GET', 'POST'])
def ajax_results():
    data = request.get_json()
    gp_request_id = data['gp_request_id']
    file_name_pb = 'pb_file_' + str(gp_request_id)
    form_params_file_name_pb = 'form_params_pb_file_' + str(gp_request_id)
    print(form_params_file_name_pb)
    with open('form_params_dir/{}.json'.format(form_params_file_name_pb)) as form_params:
        form_params_data = json.load(form_params)
    SPECIES = form_params_data['SPECIES']
    gene = form_params_data['gene']
    GCS = form_params_data['GCS']
    gene_id = get_gene_id(SPECIES, gene, GCS)['gene_id']
    try:
        transcript = form_params_data['transcript']
        transcript_id = re.sub(r' - \w*', '', transcript)
        print(transcript_id)
    except:
        return render_template('no_params.html')
    if not transcript or transcript == 'Select transcript':
        return render_template('no_params.html')
    taken_exons = form_params_data['taken_exons']
    if not taken_exons or taken_exons == []:
        return render_template('no_params.html')
    taken_exons_count = 'Exon\'s count:  {}'.format(len(taken_exons));
    SEARCH_SPECIFIC_PRIMER = form_params_data['SEARCH_SPECIFIC_PRIMER']
    CROSS_SEARCH = form_params_data['CROSS_SEARCH']
    NO_SNP = form_params_data['NO_SNP']
    SHOW_PB_LINK = form_params_data['SHOW_PB_LINK']
    PRIMER_MIN_TM_PB = form_params_data['PRIMER_MIN_TM_PB']
    PRIMER_OPT_TM_PB = form_params_data['PRIMER_OPT_TM_PB']
    PRIMER_MAX_TM_PB = form_params_data['PRIMER_MAX_TM_PB']
    PRIMER_MAX_DIFF_TM_PB = form_params_data['PRIMER_MAX_DIFF_TM_PB']
    PRIMER_MIN_SIZE_PB = form_params_data['PRIMER_MIN_SIZE_PB']
    PRIMER_OPT_SIZE_PB = form_params_data['PRIMER_OPT_SIZE_PB']
    PRIMER_MAX_SIZE_PB = form_params_data['PRIMER_MAX_SIZE_PB']
    CROSS_EXONS_MAX_SIZE_PB = form_params_data['CROSS_EXONS_MAX_SIZE_PB']
    POLYX_PB = form_params_data['POLYX_PB']
    PRIMER_PRODUCT_MIN_PB = form_params_data['PRIMER_PRODUCT_MIN_PB']
    PRIMER_PRODUCT_MAX_PB = form_params_data['PRIMER_PRODUCT_MAX_PB']
    PRIMER_MIN_GC_PB = form_params_data['PRIMER_MIN_GC_PB']
    PRIMER_MAX_GC_PB = form_params_data['PRIMER_MAX_GC_PB']
    FIVE_SAVE_EXON_DISTANCE_PB = form_params_data['FIVE_SAVE_EXON_DISTANCE_PB']
    THREE_SAVE_EXON_DISTANCE_PB = form_params_data['THREE_SAVE_EXON_DISTANCE_PB']
    F_SEARCH_DISTANCE_PB = form_params_data['F_SEARCH_DISTANCE_PB']
    R_SEARCH_DISTANCE_PB = form_params_data['R_SEARCH_DISTANCE_PB']
    MAX_MAF = form_params_data['MAX_MAF']
    auto_distances = form_params_data['auto_distances']
    params_data = {};
    with open('params_dir/{}.json'.format(file_name_pb)) as params:
        params_data = json.load(params)
    return render_template('results.html', gp_request_id=gp_request_id, SPECIES=SPECIES, transcript_id=transcript_id, gene_id=gene_id, GCS_data=GCS, data=get_primers(gp_request_id, gene, SPECIES, params_data['chromosome'], params_data['strand'], taken_exons, params_data['exons_id'], params_data['dict_exons'], SEARCH_SPECIFIC_PRIMER, CROSS_SEARCH, NO_SNP, SHOW_PB_LINK, GCS, PRIMER_MIN_TM_PB, PRIMER_OPT_TM_PB, PRIMER_MAX_TM_PB, PRIMER_MAX_DIFF_TM_PB, PRIMER_MIN_SIZE_PB, PRIMER_OPT_SIZE_PB, PRIMER_MAX_SIZE_PB, POLYX_PB, CROSS_EXONS_MAX_SIZE_PB, PRIMER_PRODUCT_MIN_PB, PRIMER_PRODUCT_MAX_PB, PRIMER_MIN_GC_PB, PRIMER_MAX_GC_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, MAX_MAF, auto_distances), gene=get_gene_id(SPECIES, gene, GCS)['gene_display_name'].upper(), transcript=transcript, taken_exons_count=taken_exons_count)


@application.route('/results_table', methods=['GET', 'POST'])
def results_table():
    gp_request_id = request.form['gp_request_id'];
    result_dict_json_name = 'result_dict_{}'.format(gp_request_id)
    with open('result_dicts_dir/{}.json'.format(result_dict_json_name)) as result_dicts_json:
        result_dict = json.load(result_dicts_json)
    taken_primers = request.form.getlist('taken_primers');
    concentration = '0,04'
    taken_primers_dict = {}
    for taken_primer_pair in taken_primers:
        taken_primer_pair = taken_primer_pair.replace('primer_checkbox_', '')
        taken_primer_pair_splitted = taken_primer_pair.split('.')
        gene = taken_primer_pair_splitted[0]
        exon = taken_primer_pair_splitted[1]
        primer_pair = taken_primer_pair_splitted[2]
        exon_number = (exon.replace('Exon_', '')).replace('Exons_', '')
        gene_exon_number = '{} {}'.format(gene, exon_number)
        f_primer_name = '{}_{}F'.format(gene, exon_number)
        r_primer_name = '{}_{}R'.format(gene, exon_number)
        taken_primers_dict[taken_primer_pair]={}
        taken_primers_dict[taken_primer_pair]['gene_exon_number']=gene_exon_number
        taken_primers_dict[taken_primer_pair]['f_primer_name']=f_primer_name
        taken_primers_dict[taken_primer_pair]['seq_F']=result_dict[exon][primer_pair]['seq_F']
        taken_primers_dict[taken_primer_pair]['r_primer_name']=r_primer_name
        taken_primers_dict[taken_primer_pair]['seq_R']=result_dict[exon][primer_pair]['seq_R']
        taken_primers_dict[taken_primer_pair]['concentration']=concentration
    print(taken_primers_dict)
    if taken_primers != []:
        return render_template('results_table.html', taken_primers_dict=taken_primers_dict)
    elif taken_primers == []:
        return '<h1 style="margin: 25% 0 0; text-align: center;">Primer pairs were not taken</h1>'


@application.route('/statuses', methods=['GET', 'POST'])
def statuses():
    time.sleep(2)
    data = request.get_json()
    gp_request_id = data['gp_request_id']
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name)) as statuses_file:
        statuses_dict = json.load(statuses_file)

    completed_statuses_list = []
    for status in statuses_dict:
        if statuses_dict[status] == 'completed':
            completed_statuses_list.append(statuses_dict[status])

    current_statuses_list = []
    for status in statuses_dict:
        current_statuses_list.append(statuses_dict[status])

    if completed_statuses_list != current_statuses_list:
        return render_template('statuses.html', statuses_dict=statuses_dict)
    elif completed_statuses_list == current_statuses_list:
        #Данный код необходим для отображения статусов при перезагрузке страницы, чтобы в дикте не был записан статус completed для всех эзонов
        statuses_file_name = 'status_of_{}'.format(gp_request_id)
        with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
            json.dump('{}', statuses_file)
        return 'Completion'


if __name__ == '__main__':
    application.run(debug=True)
