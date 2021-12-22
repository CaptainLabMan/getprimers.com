from flask import Flask, render_template, json, request
import requests, sys, time, threading, re, math
from gp_clear_code import get_transcripts, get_exons, get_primers, get_seq_len, get_reverse_complement_seq, get_gene_id, get_gc_content_f, primers_design_time_counter, pb_server_status_checker

app = Flask(__name__)


@app.route('/reverse_complement', methods=['GET', 'POST'])
def reverse_complement():
    return render_template('rev_comp.html')


@app.route('/')
@app.route('/GRCh37')
def GRCh37():
    return render_template('index.html')


@app.route('/GRCh38')
def GRCh38():
    return render_template('GRCh38.html')


@app.route('/about')
def about():
    return render_template('about.html')


@app.route('/get_len_seq', methods=['GET', 'POST'])
def get_len_seq():
    seq_for_rc = request.form['seq_for_rc'];
    return render_template('get_len_seq.html', len_seq=get_seq_len(seq_for_rc))


@app.route('/get_reverse_complement', methods=['GET', 'POST'])
def get_reverse_complement():
    seq_for_rc = request.form['seq_for_rc'];
    return render_template('reverse_complement.html', reverse_complement=get_reverse_complement_seq(seq_for_rc))


@app.route('/get_len_seq_2', methods=['GET', 'POST'])
def get_len_seq_2():
    seq_for_rc_2 = request.form['seq_for_rc_2'];
    return render_template('get_len_seq.html', len_seq=get_seq_len(seq_for_rc_2))


@app.route('/get_reverse_complement_2', methods=['GET', 'POST'])
def get_reverse_complement_2():
    seq_for_rc_2 = request.form['seq_for_rc_2'];
    return render_template('reverse_complement.html', reverse_complement=get_reverse_complement_seq(seq_for_rc_2))


@app.route('/get_gc_content_page', methods=['GET', 'POST'])
def get_gc_content_fa():
    seq_for_gc_cont = request.form['seq_for_rc'];
    return render_template('gc_cont.html', gc_content=get_gc_content_f(seq_for_gc_cont))

@app.route('/get_gc_content_page_2', methods=['GET', 'POST'])
def get_gc_content_fa_2():
    seq_for_gc_cont = request.form['seq_for_rc_2'];
    return render_template('gc_cont.html', gc_content=get_gc_content_f(seq_for_gc_cont))


@app.route('/get_transcripts', methods=['GET', 'POST'])
def get_transcripts_app():
    with open('params.json', 'w') as params:
        params.write('')
    SPECIES = request.form['SPECIES'];
    gene = request.form['gene'];
    GCS = request.form['GCS'];
    print(GCS)
    gene_id = get_gene_id(SPECIES, gene, GCS)
    print(gene_id)
    data = get_transcripts(SPECIES, gene, GCS)
    if data == 'Error. Wrong gene.':
        return render_template('wrong_gene.html', gene=gene.upper())
    else:
        return render_template('transcripts.html', data=data, gene_id=gene_id, GCS_data=GCS, SPECIES=SPECIES)


@app.route('/get_exons', methods=['GET', 'POST'])
def get_exons_app():
    transcript = request.form['transcript'];
    transcript = re.sub(r' - \w*', '', transcript)
    GCS = request.form['GCS'];
    if transcript == 'Select transcript':
        return render_template('pls_select_transcript.html')
    else:
        return render_template('exons.html', data=get_exons(transcript, GCS))


@app.route('/results', methods=['GET', 'POST'])
def loading():
    form_params_dict = {}

    SPECIES = request.form['SPECIES'];
    form_params_dict['SPECIES']=SPECIES

    gene = request.form['gene'];
    form_params_dict['gene']=gene

    GCS = request.form['GCS'];
    form_params_dict['GCS']=GCS

    transcript = request.form['transcript'];
    form_params_dict['transcript']=transcript

    taken_exons = request.form.getlist('exon');
    form_params_dict['taken_exons']=taken_exons

    SEARCH_SPECIFIC_PRIMER = request.form.getlist('SEARCH_SPECIFIC_PRIMER');
    form_params_dict['SEARCH_SPECIFIC_PRIMER']=SEARCH_SPECIFIC_PRIMER

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

    with open('form_params.json', 'w') as form_params:
        json.dump(form_params_dict, form_params)
    return render_template('loading.html', primers_design_time_counter=primers_design_time_counter(taken_exons, SEARCH_SPECIFIC_PRIMER, NO_SNP), pb_server_status=pb_server_status_checker())


@app.route('/ajax_results', methods=['GET', 'POST'])
def results():
    with open('form_params.json') as form_params:
        form_params_data = json.load(form_params)
    SPECIES = form_params_data['SPECIES']
    gene = form_params_data['gene']
    GCS = form_params_data['GCS']
    gene_id = get_gene_id(SPECIES, gene, GCS)
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
    params_data = {};
    with open('params.json') as params:
        params_data = json.load(params)
    return render_template('results.html', SPECIES=SPECIES, transcript_id=transcript_id, gene_id=gene_id, GCS_data=GCS, data=get_primers(SPECIES, params_data['chromosome'], params_data['strand'], taken_exons, params_data['exons_id'], params_data['dict_exons'], SEARCH_SPECIFIC_PRIMER, CROSS_SEARCH, NO_SNP, SHOW_PB_LINK, GCS, PRIMER_MIN_TM_PB, PRIMER_OPT_TM_PB, PRIMER_MAX_TM_PB, PRIMER_MAX_DIFF_TM_PB, PRIMER_MIN_SIZE_PB, PRIMER_OPT_SIZE_PB, PRIMER_MAX_SIZE_PB, POLYX_PB, CROSS_EXONS_MAX_SIZE_PB, PRIMER_PRODUCT_MIN_PB, PRIMER_PRODUCT_MAX_PB, PRIMER_MIN_GC_PB, PRIMER_MAX_GC_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, MAX_MAF), gene=gene.upper(), transcript=transcript, taken_exons_count=taken_exons_count)

if __name__ == '__main__':
    app.run(debug=True)


{% extends 'clear_base.html' %}


{% block main_block %}
<div id="dna_loading" class='text-center' style='margin-top: 8%;'>
  <img src="https://svgshare.com/i/Vk4.svg" style='width: 20%; margin-bottom: 4%;'>
  <h2>Please wait</h2>
  <div class="mt-2 p-2" style="margin-left: 30%; margin-right: 30%;">
    Expected waiting time - {{ primers_design_time_counter }} min.
  </div>
  {% if pb_server_status == 'overloaded' %}
  <div class="mt-2" id='server_status' style="color: red;">
    The server is processing a lot of requests. Waiting time can be significantly longer than expected.
  </div>
  {% endif %}
</div>
<div id="results_here">
</div>
<script>
            function get_results() {
                $.ajax({
                    type: "POST",
                    url: "/ajax_results",
                    type: 'POST',
                    success: function(response) {
                        $('#results_here').html(response)
                        console.log(response);
                    },
                    error: function(error) {
                        console.log(error);
                    }
                });
            window.onload = get_results;
</script>
{% endblock %}


  {% if pb_server_status == 'overloaded' %}
  <div class="mt-2" id='server_status' style="color: red;">
    The server is processing a lot of requests. Waiting time can be significantly longer than expected.
  </div>
  {% endif %}
