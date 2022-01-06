from flask import Flask, render_template, json, request, url_for, jsonify
import requests, sys, time, threading, re, math, random, cgi, geocoder
from jinja2 import Template
from gp_clear_code import get_transcripts, get_exons, get_primers, get_seq_len, get_reverse_complement_seq, get_gene_id, get_gc_content_f, pb_server_status_checker, get_clear_seq_for_textarea

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
    transcript = re.sub(r'\.\w*', '', transcript)
    transcript = re.sub(r' - \w*.\w*', '', transcript)
    GCS = request.form['GCS'];
    if transcript == 'Select transcript':
        return render_template('pls_select_transcript.html')
    else:
        return render_template('exons.html', data=get_exons(file_name_pb, transcript, GCS)['exons_count'])


@application.route('/results', methods=['GET', 'POST'])
def results():
    gp_request_id = request.form['gp_request_id'];
    print(gp_request_id)

    statuses_dict = {}
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
        json.dump('{}', statuses_file)

    form_params_dict = {}

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

    include_UTRs = request.form.getlist('include_UTRs');
    form_params_dict['include_UTRs']=include_UTRs

    split_exons = request.form.getlist('split_exons');
    form_params_dict['split_exons']=split_exons

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

    #pb_server_status=pb_server_status_checker()

    #global form_params_file_name_pb
    form_params_file_name_pb = 'form_params_pb_file_' + str(gp_request_id)
    print(form_params_file_name_pb)
    with open('form_params_dir/{}.json'.format(form_params_file_name_pb), 'w') as form_params:
        json.dump(form_params_dict, form_params)

    statuses_dict = {}
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
        json.dump(statuses_dict, statuses_file)

    return render_template('loading.html', gp_request_id=gp_request_id, gene=gene.upper())


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
    gene_data = get_gene_id(SPECIES, gene, GCS)
    gene_name = gene_data['gene_display_name'].upper()
    gene_id = gene_data['gene_id']
    try:
        transcript = form_params_data['transcript']
        transcript_id = re.sub(r' - \w*.\w*', '', transcript)
        short_transcript_id = re.sub(r'\.\w*', '', transcript)
        short_transcript_id = re.sub(r' - \w*.\w*', '', short_transcript_id)
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
    include_UTRs = form_params_data['include_UTRs']
    split_exons = form_params_data['split_exons']
    params_data = {};
    with open('params_dir/{}.json'.format(file_name_pb)) as params:
        params_data = json.load(params)
    return render_template('results.html', gp_request_id=gp_request_id, SPECIES=SPECIES, transcript_id=transcript_id, gene_id=gene_id, GCS_data=GCS, data=get_primers(gp_request_id, gene, SPECIES, params_data['chromosome'], params_data['strand'], taken_exons, params_data['exons_id'], params_data['dict_exons'], SEARCH_SPECIFIC_PRIMER, CROSS_SEARCH, NO_SNP, SHOW_PB_LINK, GCS, PRIMER_MIN_TM_PB, PRIMER_OPT_TM_PB, PRIMER_MAX_TM_PB, PRIMER_MAX_DIFF_TM_PB, PRIMER_MIN_SIZE_PB, PRIMER_OPT_SIZE_PB, PRIMER_MAX_SIZE_PB, POLYX_PB, CROSS_EXONS_MAX_SIZE_PB, PRIMER_PRODUCT_MIN_PB, PRIMER_PRODUCT_MAX_PB, PRIMER_MIN_GC_PB, PRIMER_MAX_GC_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, MAX_MAF, auto_distances, short_transcript_id, include_UTRs, split_exons), gene=gene_name, transcript=transcript, taken_exons_count=taken_exons_count)


@application.route('/results_table', methods=['GET', 'POST'])
def results_table():
    gp_request_id = request.form['gp_request_id'];
    result_dict_json_name = 'result_dict_{}'.format(gp_request_id)
    with open('result_dicts_dir/{}.json'.format(result_dict_json_name)) as result_dicts_json:
        result_dict = json.load(result_dicts_json)
    auto_naming = request.form.getlist('auto_naming');
    show_scale = request.form.getlist('show_scale');
    scale_value = request.form['scale_value'];
    tm_and_product_length = request.form.getlist('tm_and_product_length');
    name_column = request.form.getlist('name_column');
    print(auto_naming, show_scale, scale_value, tm_and_product_length, name_column)
    taken_primers = request.form.getlist('taken_primers');
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
        taken_primers_dict[taken_primer_pair]['tm_F']=result_dict[exon][primer_pair]['tm_F']
        taken_primers_dict[taken_primer_pair]['r_primer_name']=r_primer_name
        taken_primers_dict[taken_primer_pair]['seq_R']=result_dict[exon][primer_pair]['seq_R']
        taken_primers_dict[taken_primer_pair]['tm_R']=result_dict[exon][primer_pair]['tm_R']
        taken_primers_dict[taken_primer_pair]['len_amp']=result_dict[exon][primer_pair]['len_amp']
        taken_primers_dict[taken_primer_pair]['scale_value']=scale_value
    print(taken_primers_dict)
    if taken_primers != []:
        return render_template('results_table.html', taken_primers_dict=taken_primers_dict, auto_naming=auto_naming, show_scale=show_scale, tm_and_product_length=tm_and_product_length, name_column=name_column)
    elif taken_primers == []:
        return '<h1 style="margin: 25% 0 0; text-align: center;">Primer pairs were not taken</h1>'


@application.route('/statuses', methods=['GET', 'POST'])
def statuses():
    #all_statuses = ['started', 'request_to_pb', 'additional primers design', 'request_from_pb_received', 'primers_were_found', 'primers_were_not_found', 'polymorphisms_checking', 'primers_do_not_contain_polymorphisms', 'primers_contain_polymorphisms', 'completed']
    time.sleep(2)
    data = request.get_json()
    gp_request_id = data['gp_request_id']

    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    try:
        with open('statuses/{}.json'.format(statuses_file_name)) as statuses_file:
            statuses_dict = json.load(statuses_file)
    except:
        print('Error in getting statuses_dict 1')
        get_dict = False
        while get_dict == False:
            try:
                with open('statuses/{}.json'.format(statuses_file_name)) as statuses_file:
                    statuses_dict = json.load(statuses_file)
                get_dict = True
            except:
                print('Error in getting statuses_dict 2')


    data_file_name = 'data_of_{}'.format(gp_request_id)
    try:
        with open('data/{}.json'.format(data_file_name)) as data_file:
            data_dict = json.load(data_file)
    except:
        print('Some error occured while opening data.json file 1 (application.py)')
        success = False
        while success == False:
            try:
                with open('data/{}.json'.format(data_file_name)) as data_file:
                    data_dict = json.load(data_file)
                success = True
            except:
                print('Some error occured while opening data.json file 2 (application.py)')
                time.sleep(1)
    taken_exons = data_dict['taken_exons']
    taken_exons_count = len(taken_exons)

    design_completed = False
    completed_exons = []
    try:
        for status in statuses_dict:
            if statuses_dict[status] == 'completed':
                completed_exons.append(statuses_dict[status])
    except:
        print('Some error in statuses 1')
        success = False
        while success == False:
            try:
                with open('statuses/{}.json'.format(statuses_file_name)) as statuses_file:
                    statuses_dict = json.load(statuses_file)
                completed_exons = []
                for status in statuses_dict:
                    if statuses_dict[status] == 'completed':
                        completed_exons.append(statuses_dict[status])
                success = True
                print('success')
            except:
                print('Some error in statuses 2')
    completed_exons_count = len(completed_exons)

    if taken_exons_count == completed_exons_count:
        design_completed = True
        print('Design comleted')
    else:
        print('Design is not comleted yet')

    form_params_file_name_pb = 'form_params_pb_file_' + str(gp_request_id)
    with open('form_params_dir/{}.json'.format(form_params_file_name_pb)) as form_params:
        form_params_data = json.load(form_params)

    def progress_bar(gp_request_id, statuses_dict, taken_exons_count, form_params_data):
        NO_SNP = form_params_data['NO_SNP']
        primers_design_completed_statuses_list = []
        primers_design_completed_statuses = ['primers_were_found', 'primers_were_not_found', 'polymorphisms_checking', 'primers_do_not_contain_polymorphisms', 'primers_contain_polymorphisms', 'completed']
        polymorphisms_checking_completed_statuses_list = []
        polymorphisms_checking_completed_statuses = ['primers_do_not_contain_polymorphisms', 'primers_contain_polymorphisms', 'completed']
        exon_started_statuses_list = []
        exon_started_statuses = ['request_to_pb', 'request_from_pb_received', 'primers_were_found', 'primers_were_not_found', 'polymorphisms_checking', 'primers_do_not_contain_polymorphisms', 'primers_contain_polymorphisms', 'completed']
        additional_primers_design = 'additional primers design'
        for status in statuses_dict:
            if statuses_dict[status] in primers_design_completed_statuses:
                primers_design_completed_statuses_list.append(statuses_dict[status])
            if statuses_dict[status] in polymorphisms_checking_completed_statuses:
                polymorphisms_checking_completed_statuses_list.append(statuses_dict[status])
            if statuses_dict[status] in exon_started_statuses or additional_primers_design in statuses_dict[status]:
                exon_started_statuses_list.append(statuses_dict[status])
        primers_design_completed_statuses_list_len = len(primers_design_completed_statuses_list)
        polymorphisms_checking_completed_statuses_list_len = len(polymorphisms_checking_completed_statuses_list)
        exon_started_statuses_list_len = len(exon_started_statuses_list)

        if NO_SNP == []:
            progress = int((exon_started_statuses_list_len * 9) / taken_exons_count) + int((primers_design_completed_statuses_list_len * 91) / taken_exons_count)
        elif NO_SNP == ['checked']:
            progress = int((exon_started_statuses_list_len * 9) / taken_exons_count) + int((primers_design_completed_statuses_list_len * 61) / taken_exons_count) + int((polymorphisms_checking_completed_statuses_list_len * 30) / taken_exons_count)
        print(progress)

        return progress

    try:
        progress = progress_bar(gp_request_id, statuses_dict, taken_exons_count, form_params_data)
    except:
        progress = 0
        print('Some error in progress function')

    def pb_server_status_for_users(form_params_data, data_dict):
        SEARCH_SPECIFIC_PRIMER = form_params_data['SEARCH_SPECIFIC_PRIMER']
        NO_SNP = form_params_data['NO_SNP']
        pb_server_status_for_each_exon = data_dict['pb_server_status']

        pb_server_status = ''
        for exon_status in pb_server_status_for_each_exon:
            if pb_server_status_for_each_exon[exon_status] == 'overloaded':
                pb_server_status = 'overloaded'
                break
            elif pb_server_status_for_each_exon[exon_status] == 'ok':
                pb_server_status = 'ok'

        return pb_server_status

    try:
        pb_server_status = pb_server_status_for_users(form_params_data, data_dict)
        #expected_waiting_time = server_status_and_time['expected_waiting_time']
        print(pb_server_status)
        #print(expected_waiting_time)
    except:
        pb_server_status = ''
        #expected_waiting_time = 'wait'
        print('PB status - no status yet')

    def expected_waiting_time_for_users(form_params_data, data_dict, taken_exons_count):
        SEARCH_SPECIFIC_PRIMER = form_params_data['SEARCH_SPECIFIC_PRIMER']
        NO_SNP = form_params_data['NO_SNP']

        if taken_exons_count == 1:
            time = 2
        else:
            time = 1
        additional_time = math.ceil(taken_exons_count / 10)

        expected_waiting_time = 0
        if SEARCH_SPECIFIC_PRIMER == ['checked'] and NO_SNP == ['checked']:
            expected_waiting_time = (time * taken_exons_count) + additional_time
        elif SEARCH_SPECIFIC_PRIMER == ['checked'] and NO_SNP == []:
            expected_waiting_time = math.ceil((time * 0.7 * taken_exons_count)) + additional_time
        elif SEARCH_SPECIFIC_PRIMER == [] and NO_SNP == ['checked']:
            expected_waiting_time = math.ceil((time * 0.5 * taken_exons_count)) + additional_time
        elif SEARCH_SPECIFIC_PRIMER == [] and NO_SNP == []:
            expected_waiting_time = math.ceil(0.025 * taken_exons_count)

        return expected_waiting_time

    try:
        expected_waiting_time = expected_waiting_time_for_users(form_params_data, data_dict, taken_exons_count)
        print(expected_waiting_time)
    except:
        expected_waiting_time = 'wait'
        print('Expected time - no time yet')

    if design_completed == False:
        return render_template('statuses.html', statuses_dict=statuses_dict, progress=progress, pb_server_status=pb_server_status, expected_waiting_time=expected_waiting_time)
    elif design_completed == True:
        #Данный код необходим для отображения статусов при перезагрузке страницы, чтобы в дикте не был записан статус completed для всех эзонов
        statuses_file_name = 'status_of_{}'.format(gp_request_id)
        with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
            json.dump('{}', statuses_file)
        return 'Completion'


@application.route('/tests', methods=['GET', 'POST'])
def tests():
    client_ip = request.environ.get('HTTP_X_REAL_IP', request.remote_addr)
    g = geocoder.ip('{}'.format(client_ip))
    g = g.json
    tests_dict = {}
    tests_dict['client_ip']=client_ip
    tests_dict['g']=g

    return tests_dict


if __name__ == '__main__':
    application.run(debug=True)
