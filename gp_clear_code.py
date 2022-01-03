
def state_statuses(statuses_dict, exon, gp_request_id, state):
    import json
    statuses_dict[exon]=state
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
        json.dump(statuses_dict, statuses_file)

def cyc_post_req_for_ens(server, link, headers, data):
    import requests, urllib3, time
    req_status = False
    try:
        request = requests.post(server+link, data=data, headers=headers, timeout=15.0)
    except Exception:
        print('An exception 1 occured in Ensembl post request')
        while req_status != True:
            time.sleep(2.5)
            try:
                request = requests.post(server+link, data=data, headers=headers, timeout=15.0)
                req_status = True
            except Exception:
                print('An exception 2 occured in Ensembl post request')
    return request


def cyc_post_req(link, data):
    import requests, urllib3, time
    req_status = False
    try:
        request = requests.post(link, data=data, timeout=30.0)
    except Exception:
        print('An exception 1 occured in Primer-BLAST post request')
        while req_status != True:
            time.sleep(2.5)
            try:
                request = requests.post(link, data=data, timeout=30.0)
                print('Попытка запроса 2')
                req_status = True
            except Exception:
                print('An exception 2 occured in Primer-BLAST post request')
    return request


def cyc_get_req_for_pb(link):
    import requests, urllib3, time
    req_status = False
    try:
        request = requests.get(link, timeout=15.0)
    except Exception:
        print('An exception 1 occured in Primer-BLAST get request')
        counter = 1
        while req_status != True:
            time.sleep(2.5)
            counter += 1
            try:
                request = requests.get(link, timeout=15.0)
                print('Попытка запроса {}'.format(counter))
                req_status = True
            except Exception:
                print('An exception 2 occured in Primer-BLAST get request')
    return request


def cyc_get_req_for_ens(server, link, headers):
    import requests, urllib3, time
    req_status = False
    try:
        request = requests.get(server+link, headers=headers, timeout=15.0)
    except Exception:
        print('An exception 1 occured in Ensembl request')
        while req_status != True:
            time.sleep(2.5)
            try:
                request = requests.get(server+link, headers=headers, timeout=15.0)
                req_status = True
            except Exception:
                print('An exception 2 occured in Ensembl request')
    return request


def get_gene_id(SPECIES, gene, GCS):
    import requests, sys, json

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"

    gene_path = '/xrefs/symbol/{}/{}?'.format(SPECIES, gene)
    headers={"Content-Type" : "application/json"}
    #get_gene_id = requests.get(server+gene_path, headers={ "Content-Type" : "application/json"})
    #get_gene_id = cycreq(requests.get(server+gene_path, headers={ "Content-Type" : "application/json"}, timeout=None))
    get_gene_id = cyc_get_req_for_ens(server, gene_path, headers)


    #Тут условие, что если введен неверный ген, то get_gene_id_decoded == [], то будет вызван близайший return и на сайт будет подаваться исключение.
    get_gene_id_decoded = get_gene_id.json()
    if get_gene_id_decoded == []:
        return 'Error. Wrong gene.'
    gene_id = ''
    for item in get_gene_id_decoded:
        if item['type'] == 'gene':
            gene_check_info_path = "/lookup/id/{}?".format(item['id'])
            headers={"Content-Type" : "application/json"}
            #gene_check_info_response = requests.get(server+gene_check_info_path, headers={ "Content-Type" : "application/json"})
            #gene_check_info_response = cycreq(requests.get(server+gene_check_info_path, headers={ "Content-Type" : "application/json"}, timeout=None))
            gene_check_info_response = cyc_get_req_for_ens(server, gene_check_info_path, headers)
            check_gene_name = gene_check_info_response.json()
            if check_gene_name['display_name'].upper() == gene.upper():
                gene_id = item['id']
                break
    #print('--------------------------\n')

    if gene_id == '':
        return 'Error. Wrong gene.'
    else:
        gene_dict = {}
        gene_dict['gene_display_name']=check_gene_name['display_name']
        gene_dict['gene_id']=gene_id
        return gene_dict



def get_transcripts(SPECIES, gene, GCS):
    import requests, sys, json, re

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"

    gene_path = '/xrefs/symbol/{}/{}?'.format(SPECIES, gene)
    headers={"Content-Type" : "application/json"}
    #get_gene_id = requests.get(server+gene_path, headers={ "Content-Type" : "application/json"})
    #get_gene_id = cycreq(requests.get(server+gene_path, headers={ "Content-Type" : "application/json"}, timeout=None))
    get_gene_id = cyc_get_req_for_ens(server, gene_path, headers)


    #Тут условие, что если введен неверный ген, то get_gene_id_decoded == [], то будет вызван близайший return и на сайт будет подаваться исключение.
    get_gene_id_decoded = get_gene_id.json()
    gene_id = ''
    print('\n-------------------------')
    print('\nEntered gene name - {}\n'.format(gene.upper()))
    for item in get_gene_id_decoded:
        if item['type'] == 'gene':
            gene_check_info_path = "/lookup/id/{}?".format(item['id'])
            headers={"Content-Type" : "application/json"}
            #gene_check_info_response = requests.get(server+gene_check_info_path, headers={ "Content-Type" : "application/json"})
            #gene_check_info_response = cycreq(requests.get(server+gene_check_info_path, headers={ "Content-Type" : "application/json"}, timeout=None))
            gene_check_info_response = cyc_get_req_for_ens(server, gene_check_info_path, headers)
            check_gene_name = gene_check_info_response.json()
            print('Found gene name - {}'.format(check_gene_name['display_name'].upper()))
            print('{} gene ID - {}\n'.format(check_gene_name['display_name'].upper(), item['id']))
            #print(check_gene_name)
            if check_gene_name['display_name'].upper() == gene.upper() or gene.upper() in check_gene_name['display_name'].upper():
                gene_id = item['id']
                break
    print('--------------------------\n')


    # Эта чатсь может доставать инфу по транскриптам в зависимости от features
    transcripts_path = "/overlap/id/{}?feature=transcript".format(gene_id)
    headers={"Content-Type" : "application/json"}
    #get_transcripts_id = requests.get(server+transcripts_path, headers={ "Content-Type" : "application/json"})
    #get_transcripts_id = cycreq(requests.get(server+transcripts_path, headers={ "Content-Type" : "application/json"}, timeout=None))
    get_transcripts_id = cyc_get_req_for_ens(server, transcripts_path, headers)

    get_transcripts_id_decoded = get_transcripts_id.json()

    transcripts_list = []
    for element in get_transcripts_id_decoded:
        if element['Parent'] == gene_id:
            transcript_and_version = '{}.{}'.format(element['transcript_id'], element['version'])
            transcripts_list.append(transcript_and_version)


    transcripts_and_refseq_list = []
    # Часть по вытаскиванию RefSeq и Часть по выводу на экран инфы по транскриптан с рефсек.
    #Для увеличения скорости работы, можно сделать несколько потоков (только вот скорее всего сервера энсембла не позволят делать так много запросов практически одномоментно).
    for transcript_id_element in transcripts_list:
        transcript_id_element_for_req = re.sub(r'\.\w*', '', transcript_id_element)
        transcript_id_element_for_req = re.sub(r' - \w*.\w*', '', transcript_id_element_for_req)
        ref_seq_path = "/xrefs/id/{}?".format(transcript_id_element_for_req)
        headers={"Content-Type" : "application/json"}
        #get_ref_seq = requests.get(server+ref_seq_path, headers={ "Content-Type" : "application/json"})
        #get_ref_seq = cycreq(requests.get(server+ref_seq_path, headers={ "Content-Type" : "application/json"}, timeout=None))
        get_ref_seq = cyc_get_req_for_ens(server, ref_seq_path, headers)

        get_ref_seq_decoded = get_ref_seq.json()
        for dict_element in get_ref_seq_decoded:
            if dict_element['dbname'] == 'RefSeq_mRNA':
                #primary_id = dict_element['primary_id']
                primary_id = dict_element['display_id']
                transcripts_and_refseq_list.append('{} - {}'.format(transcript_id_element, primary_id))
                if transcript_id_element in transcripts_list:
                    transcripts_list.remove(transcript_id_element)  # при наличии у транскрипта refseq он будет удалять из списка transcripts чтобы этот транскрипт не дублировался

    #Здесь объединяю списки транскриптов без РефСек и с РефСек.
    for singe_transcript_and_refseq in transcripts_and_refseq_list:
        transcripts_list.append(singe_transcript_and_refseq)
    #Переворачиваю лист транскриптов чтобы транскрипты с NM были вверху, а не внизу.
    transcripts_list.reverse()

    return transcripts_list



def get_exons(file_name_pb, transcript_id, GCS):
    import requests, sys, json

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"

    # Эта чатсь может доставать инфу по экзонам в зависимости от features
    exons_path = "/overlap/id/{}?feature=exon".format(transcript_id)
    headers={"Content-Type" : "application/json"}
    #get_exons = requests.get(server+exons_path, headers={ "Content-Type" : "application/json"})
    #get_exons = cycreq(requests.get(server+exons_path, headers={ "Content-Type" : "application/json"}, timeout=None))
    get_exons = cyc_get_req_for_ens(server, exons_path, headers)

    get_exons_decoded = get_exons.json()

    strand = ''
    chromosome = ''

    exons_id = []
    for element in get_exons_decoded:
        if element['Parent'] == transcript_id:
            exons_id.append(element['exon_id'])
            strand = element['strand']               # Определение направления цепи
            chromosome = element['seq_region_name']  # Определение [хромосомы]

    #Рассчитывает количество экзонов в транскрипте для передаче html шаблонизатору в exons.html
    exons_count = []
    for nubmer in range(1, len(exons_id) + 1):
        exons_count.append(nubmer)

    dict_exons = dict()
    for element in get_exons_decoded:
        if element['Parent'] == transcript_id:
            dict_exons[element['exon_id']] = {'start': element['start'], 'end': element['end']}

    #Формирование дикта для сериализации данных в json
    params_json = {'chromosome': chromosome, 'strand': str(strand), 'exons_id': exons_id, 'dict_exons': dict_exons, 'exons_count': exons_count}
    #print(params_json)
    #Запись номера хромосомы и цепи в файл params.json
    #print(file_name_pb)
    with open('params_dir/{}.json'.format(file_name_pb), 'w') as params:
        json.dump(params_json, params)

    return params_json



def get_primers(gp_request_id, gene_name, SPECIES, chromosome, strand, taken_exons, exons_id, dict_exons, SEARCH_SPECIFIC_PRIMER, CROSS_SEARCH, NO_SNP, SHOW_PB_LINK, GCS, PRIMER_MIN_TM_PB, PRIMER_OPT_TM_PB, PRIMER_MAX_TM_PB, PRIMER_MAX_DIFF_TM_PB, PRIMER_MIN_SIZE_PB, PRIMER_OPT_SIZE_PB, PRIMER_MAX_SIZE_PB, POLYX_PB, CROSS_EXONS_MAX_SIZE_PB, PRIMER_PRODUCT_MIN_PB, PRIMER_PRODUCT_MAX_PB, PRIMER_MIN_GC_PB, PRIMER_MAX_GC_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, MAX_MAF, auto_distances, short_transcript_id, include_UTRs, split_exons):
    import requests, sys, time, threading, re, json, math, copy
    from bs4 import BeautifulSoup
    from gp_clear_code import get_reverse_complement_seq, get_exons

    statuses_dict = {}
    statuses_file_name = 'status_of_{}'.format(gp_request_id)
    with open('statuses/{}.json'.format(statuses_file_name), 'w') as statuses_file:
        json.dump('{}', statuses_file)

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"


    gene_name = gene_name.upper()

    # Формирование списка id выбранных экзонов
    taken_exons_id = []
    for element in taken_exons:
        taken_exons_id.append(exons_id[int(element) - 1])

    #Здесь все выбранные экзоны собираются в словарь с ключом в виде НОМЕРА ЭКЗОНА и значением в виде координат начала и конца экзона.
    taken_exons_list = []
    taken_exons_dict = {}
    for element in taken_exons:
        taken_exons_list.append('Exon_{}'.format(int(element)))
        taken_exons_dict['Exon_{}'.format(int(element))]=dict_exons[exons_id[int(element) - 1]]
    taken_exons_id = taken_exons_list
    dict_exons = taken_exons_dict
    print(taken_exons_id)
    print(dict_exons)

    #Получение данных о расположении UTRs транскрипта и изменение позиций первого и последнего экзонов если задан чекбокс
    if include_UTRs != ['checked']:
        link = "/lookup/id/{}?expand=1;utr=1".format(short_transcript_id)
        headers = {"Content-Type" : "application/json"}
        UTRs_request = cyc_get_req_for_ens(server, link, headers)
        UTRs_request = UTRs_request.json()
        UTRs_request = UTRs_request['UTR']
        UTRs_dict = {}
        for UTR in UTRs_request:
            if UTR['object_type'] == 'five_prime_UTR':
                UTRs_dict['five_prime_UTR'] = UTR
            elif UTR['object_type'] == 'three_prime_UTR':
                UTRs_dict['three_prime_UTR'] = UTR

        for exon in dict_exons:
            if strand == '-1':
                if dict_exons[exon]['end'] == UTRs_dict['five_prime_UTR']['end']:
                    dict_exons[exon]['end'] = UTRs_dict['five_prime_UTR']['start'] - 1
                if dict_exons[exon]['start'] == UTRs_dict['three_prime_UTR']['start']:
                    dict_exons[exon]['start'] = UTRs_dict['three_prime_UTR']['end'] + 1
            elif strand == '1':
                if dict_exons[exon]['start'] == UTRs_dict['five_prime_UTR']['start']:
                    dict_exons[exon]['start'] = UTRs_dict['five_prime_UTR']['end'] + 1
                if dict_exons[exon]['end'] == UTRs_dict['three_prime_UTR']['end']:
                    dict_exons[exon]['end'] = UTRs_dict['three_prime_UTR']['start'] - 1

    #Здесь происходит формирование списка и словаря для cross_search, то ест здесь объединяются экзоны с максимальной заданной длиной - cross_search_max_distance
    if CROSS_SEARCH == ['checked']:
        cross_search_distance = 0
        cross_search_max_distance = int(CROSS_EXONS_MAX_SIZE_PB)
        cross_search_taken_exons_dict = {}
        cross_search_taken_exons_list = []
        for position, element in enumerate(taken_exons_list):
            selection_list = []
            selection_dict = {}
            exon_delta = 1
            while exon_delta < (len(taken_exons_list) - position):
                if strand == '-1':
                    dif = taken_exons_dict[element]['end'] - taken_exons_dict[taken_exons_list[position + exon_delta]]['start']
                elif strand == '1':
                    dif = taken_exons_dict[element]['start'] - taken_exons_dict[taken_exons_list[position + exon_delta]]['end']
                    dif = int(math.sqrt(dif ** 2)) #Это необходимо так как в случае +1 цепи хромосомы, в ходе вычитания получаются отрицательные значения.
                if int(dif) <= int(cross_search_max_distance) and int(dif) > 0:
                    selection_list.append(dif)
                    selection_dict[dif]=taken_exons_list[position + exon_delta]
                exon_delta += 1
            if selection_list != []:
                selection_list_max = max(selection_list)
                if strand == '-1':
                    cross_search_taken_exons_dict['Exons_{}-{}'.format(re.sub(r'Exon_', '', element), re.sub(r'Exon_', '', selection_dict[selection_list_max]))]={'start': taken_exons_dict[element]['end'] - selection_list_max, 'end': taken_exons_dict[element]['end']}
                elif strand == '1':
                    cross_search_taken_exons_dict['Exons_{}-{}'.format(re.sub(r'Exon_', '', element), re.sub(r'Exon_', '', selection_dict[selection_list_max]))]={'start': taken_exons_dict[element]['start'], 'end': taken_exons_dict[element]['start'] + selection_list_max}
                to_delete = []
                for element in range(position + 1, int(taken_exons_list.index(selection_dict[selection_list_max]) + 1)):
                    to_delete.append(taken_exons_list[element])
                for element in to_delete:
                    if element in taken_exons_list:
                        taken_exons_list.remove(element)
            elif selection_list == []:
                cross_search_taken_exons_dict[element]=taken_exons_dict[element]
        for element in cross_search_taken_exons_dict:
            cross_search_taken_exons_list.append(element)
        taken_exons_id = cross_search_taken_exons_list
        dict_exons = cross_search_taken_exons_dict
    print(dict_exons)

    #Разбивание больших экзонов на части
    if split_exons == ['checked']:
        print(dict_exons)
        new_taken_exons_list = []
        dict_exons_deepcopy = copy.deepcopy(dict_exons)
        max_exon_size = int(CROSS_EXONS_MAX_SIZE_PB)
        print(dict_exons)
        for exon in dict_exons:
            exon_len = dict_exons[exon]['end'] - dict_exons[exon]['start'] + 1
            if exon_len > max_exon_size:
                divider = math.ceil(exon_len / max_exon_size)
                optimal_exons_size = math.ceil(exon_len / divider)
                distance_counter = dict_exons[exon]['start']
                exons_part_counter = 1
                stop_flag = False
                while stop_flag == False:
                    if strand == '1':
                        if (distance_counter + optimal_exons_size) < dict_exons[exon]['end']:
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]={}
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]['start']=distance_counter
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]['end']=distance_counter + optimal_exons_size
                        else:
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]={}
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]['start']=distance_counter
                            dict_exons_deepcopy['{}-{}'.format(exon, exons_part_counter)]['end']=dict_exons[exon]['end']
                            stop_flag = True
                    elif strand == '-1':
                        if (distance_counter + optimal_exons_size) < dict_exons[exon]['end']:
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]={}
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]['start']=distance_counter
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]['end']=distance_counter + optimal_exons_size
                        else:
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]={}
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]['start']=distance_counter
                            dict_exons_deepcopy['{}-{}'.format(exon, divider + 1 - exons_part_counter)]['end']=dict_exons[exon]['end']
                            stop_flag = True
                    distance_counter += optimal_exons_size + 1
                    exons_part_counter += 1
                del dict_exons_deepcopy[exon]
                print(exon_len)
                print(divider)
                print(optimal_exons_size)
        for exon in dict_exons_deepcopy:
            new_taken_exons_list.append(exon)
        taken_exons_id = new_taken_exons_list
        dict_exons = dict_exons_deepcopy
        print(taken_exons_id)
        print(dict_exons)


    # Часть для получения последовательностей для ПраймерБласта + Часть по рассчетам для вставки в Праймер бласт
    print('\nPrimer design is going.\nPlease wait...\n')
    print('Count of exons (threads) - {}\n'.format(len(taken_exons)))


    result_dict = {}

    info_dict = {}

    def get_get_primers(element, result_dict, statuses_dict, gp_request_id, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, info_dict, auto_distances):

        state_statuses(statuses_dict, element, gp_request_id, 'started')
        #Здесь начало авто-дистант поиска
        if auto_distances == ['checked'] and SEARCH_SPECIFIC_PRIMER == ['checked']:
            F_SEARCH_DISTANCE_PB = 125
            R_SEARCH_DISTANCE_PB = 125

        def primers_search(element, result_dict, statuses_dict, gp_request_id, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, info_dict, request_status, auto_distances):
            primers_search_dict = {}
            #print('\nlook_here_2')
            #print(auto_distances)

            info_dict[element]={}
            info_dict[element]['search_distances']={}

            #Преобразование входящих переменных
            F_SEARCH_DISTANCE_PB = int(F_SEARCH_DISTANCE_PB)
            #print('look here')
            #print(F_SEARCH_DISTANCE_PB)
            R_SEARCH_DISTANCE_PB = int(R_SEARCH_DISTANCE_PB)
            FIVE_SAVE_EXON_DISTANCE_PB = int(FIVE_SAVE_EXON_DISTANCE_PB)
            THREE_SAVE_EXON_DISTANCE_PB = int(THREE_SAVE_EXON_DISTANCE_PB)

            info_dict[element]['search_distances']['F_SEARCH_DISTANCE_PB']=F_SEARCH_DISTANCE_PB
            info_dict[element]['search_distances']['R_SEARCH_DISTANCE_PB']=R_SEARCH_DISTANCE_PB
            info_dict[element]['search_distances']['FIVE_SAVE_EXON_DISTANCE_PB']=FIVE_SAVE_EXON_DISTANCE_PB
            info_dict[element]['search_distances']['THREE_SAVE_EXON_DISTANCE_PB']=THREE_SAVE_EXON_DISTANCE_PB

            exon_number_for_thread = element
            print('Thread for {} started'.format(exon_number_for_thread))
            #Тут по сути собирается запрос для получения последовательности для поиска праймеров (последовательонсти вместе с отступом для поиска).
            if strand == '1':
                begin_of_search_seq = int(dict_exons[element]['start']) - F_SEARCH_DISTANCE_PB
                end_of_search_seq = int(dict_exons[element]['end']) + R_SEARCH_DISTANCE_PB
            elif strand == '-1':
                begin_of_search_seq = int(dict_exons[element]['start']) - R_SEARCH_DISTANCE_PB
                end_of_search_seq = int(dict_exons[element]['end']) + F_SEARCH_DISTANCE_PB
            primers_search_dict['begin_of_search_seq']=begin_of_search_seq
            primers_search_dict['end_of_search_seq']=end_of_search_seq
            exon_and_distance_path = "/sequence/region/{}/{}:{}..{}:{}?expand_3prime=0;expand_5prime=0".format(SPECIES, chromosome, begin_of_search_seq, end_of_search_seq, strand)
            headers={"Content-Type" : "text/plain"}
            #get_exon_and_distance = requests.get(server+exon_and_distance_path, headers={ "Content-Type" : "text/plain"})
            #get_exon_and_distance = cycreq(requests.get(server+exon_and_distance_path, headers={ "Content-Type" : "text/plain"}, timeout=None))
            get_exon_and_distance = cyc_get_req_for_ens(server, exon_and_distance_path, headers)

            #Здесь часть похожая на часть сверху, по получению последовательности, но это последовательность строго самого кэзона (+UTR если он есть).
            #Нужна эта часть для того чтобы выделить большими буквами экзон, а маленькими интроны в последовательности.
            exon_only_path = "/sequence/region/{}/{}:{}..{}:{}?expand_3prime=0;expand_5prime=0".format(SPECIES, chromosome, int(dict_exons[element]['start']), int(dict_exons[element]['end']), strand)
            headers={"Content-Type" : "text/plain"}
            #get_exon_only = requests.get(server+exon_only_path, headers={ "Content-Type" : "text/plain"})
            #get_exon_only = cycreq(requests.get(server+exon_only_path, headers={ "Content-Type" : "text/plain"}, timeout=None))
            get_exon_only = cyc_get_req_for_ens(server, exon_only_path, headers)
            get_exon_only = get_exon_only.text
            primers_search_dict['get_exon_only']=get_exon_only

            #Строки ниже необходимы просто для вывода инфы о дистанции для поиска праймеров в консоль.
            distance = '1 - {}\n{} - {}'.format(int(F_SEARCH_DISTANCE_PB) - int(FIVE_SAVE_EXON_DISTANCE_PB), len(get_exon_and_distance.text) - (int(R_SEARCH_DISTANCE_PB) - int(THREE_SAVE_EXON_DISTANCE_PB)), len(get_exon_and_distance.text))
            print_dis = distance
            primers_search_dict['print_dis']=print_dis
            number_of_exon = element   # Номер экзона
            primers_search_dict['number_of_exon']=number_of_exon
            print_exon_number = '\n{}'.format(number_of_exon)     #Перенесено в конец с часть по принту
            primers_search_dict['print_exon_number']=print_exon_number

            #Переменная ниже - это последователность на которой будут искаться праймеры вместе с отступами для поиска праймеров.
            #Ниже две новые переменные (exon_sequence_for_users и seq_for_pb), которые несут в себе одно и то же, просто нужны для разных целей. В принципе можно и одну сдлеать.
            exon_sequence_for_users = (get_exon_and_distance.text).lower()
            exon_sequence_for_users = exon_sequence_for_users.replace(get_exon_only.lower(), get_exon_only.upper())
            primers_search_dict['exon_sequence_for_users']=exon_sequence_for_users
            exon_sequence = {'sequence': exon_sequence_for_users}  #Создание словаря, в котром будет храниться последовательность ЭКЗОНА.
            print(exon_sequence)

            seq_for_pb = '\'{}\''.format(get_exon_and_distance.text)
            if 'N' in seq_for_pb or 'n' in seq_for_pb:
                print('Sequence contains \'N\'')
                return 'Sequence contains \'N\''

            #Часть по рассчету позиций поиска праймеров для вставки в форму PB.
            PRIMER5_END_PB = int(F_SEARCH_DISTANCE_PB) - int(FIVE_SAVE_EXON_DISTANCE_PB)
            PRIMER3_START_PB = len(get_exon_and_distance.text) - (int(R_SEARCH_DISTANCE_PB) - int(THREE_SAVE_EXON_DISTANCE_PB))
            PRIMER3_END_PB = len(get_exon_and_distance.text)

            #Преобразование длины ампликона в минимально возможный если он меньше
            CROSS_SEARCH_LEN_ORDER = ''
            if CROSS_SEARCH == ['checked'] and CROSS_SEARCH_LEN_ORDER == ['checked']:
                min_possible_amplicon_size = int(len(get_exon_only)) + FIVE_SAVE_EXON_DISTANCE_PB + THREE_SAVE_EXON_DISTANCE_PB + (int(PRIMER_MAX_SIZE_PB) * 2)
                if int(PRIMER_PRODUCT_MAX_PB) < min_possible_amplicon_size:
                    PRIMER_PRODUCT_MAX_PB_remake = min_possible_amplicon_size
                else:
                    PRIMER_PRODUCT_MAX_PB_remake = PRIMER_PRODUCT_MAX_PB
            else:
                PRIMER_PRODUCT_MAX_PB_remake = PRIMER_PRODUCT_MAX_PB

            #print(PRIMER_PRODUCT_MAX_PB)
            #print(PRIMER_PRODUCT_MAX_PB_remake)

            # Добавление пустого дикта экзона чтобы в случае разрыва сети он не пропадал из результ дикта, потому что следующее добавление идет уже после получения ответа от ПБ
            result_dict['{}'.format(element)]={}

            #Часть по получению тупо страницы HTML с праймерами от PrimerBlast'a
            #Сначала я загоняю последовательность со всеми параметрами в ПБ и он возвращает мне ссылку
            #Датее эту ссылку я вставляю во второй запрос (r2) и он уже дает мне HTML страницу
            #которую я буду парсить на наличие праймеров (естественно со слипом)
            #refseq_mrna
            #PRIMERDB/genome_selected_species
            url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
            data = {
                'INPUT_SEQUENCE': seq_for_pb,
                'PRIMER5_START': '1',
                'PRIMER5_END': PRIMER5_END_PB,
                'PRIMER3_START': PRIMER3_START_PB,
                'PRIMER3_END': PRIMER3_END_PB,
                'PRIMER_SPECIFICITY_DATABASE': 'PRIMERDB/genome_selected_species',
                'SEARCHMODE': '1',
                'ORGANISM': SPECIES.replace('_', ' '),
                'TOTAL_PRIMER_SPECIFICITY_MISMATCH': '2',
                'PRIMER_3END_SPECIFICITY_MISMATCH': '1',
                'MISMATCH_REGION_LENGTH': '5',
                'TOTAL_MISMATCH_IGNORE': '6',
                'MAX_CANDIDATE_PRIMER': '3000',
                'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
                'HITSIZE': '100000',
                'EVALUE': '30000',
                'WORD_SIZE': '7',
                'LOW_COMPLEXITY_FILTER': 'checked',
                'PRIMER_MIN_SIZE': PRIMER_MIN_SIZE_PB,
                'PRIMER_OPT_SIZE': PRIMER_OPT_SIZE_PB,
                'PRIMER_MAX_SIZE': PRIMER_MAX_SIZE_PB,
                'POLYX': POLYX_PB,
                'PRIMER_MIN_TM': PRIMER_MIN_TM_PB,
                'PRIMER_OPT_TM': PRIMER_OPT_TM_PB,
                'PRIMER_MAX_TM': PRIMER_MAX_TM_PB,
                'PRIMER_MAX_DIFF_TM': PRIMER_MAX_DIFF_TM_PB,
                'PRIMER_PRODUCT_MIN': PRIMER_PRODUCT_MIN_PB,
                'PRIMER_PRODUCT_MAX': PRIMER_PRODUCT_MAX_PB_remake,
                'PRIMER_MIN_GC': PRIMER_MIN_GC_PB,
                'PRIMER_MAX_GC': PRIMER_MAX_GC_PB,
                'SEARCH_SPECIFIC_PRIMER': SEARCH_SPECIFIC_PRIMER
            }

            #r = requests.post(url, data=data)
            r = cyc_post_req(url, data)
            headers_pb = r.headers

            if request_status == 'first_request':
                state = 'request_to_pb'
            elif 'additional primers design' in request_status:
                state = request_status
            state_statuses(statuses_dict, element, gp_request_id, state)

            #Здесь условная конструкция, которая если SEARCH_SPECIFIC_PRIMER = [], говорит о том, что выбран неспецифический подбор праймеров
            #а значит нет смысла брать ссылку из первого запроса на второй запрос, который специфически подбирает праймеры. В этом случае второй запрос, который далее
            #необходим для выполнения алгоритма приравнивается к значению первого запроса, просто потому что второй нужен для дальнейшей работы алгоритма, как я уже
            #и написал.
            #В случае же если SEARCH_SPECIFIC_PRIMER = ['checked'], то из хедера первого запроса будет получена ссылка на страницу специфического подбора праймеров, которая
            #далее будет использована во втором запросе для получения хтмл страницы со специфиечскими праймерами.
            pb_link_checker = '' # Нужно для добавления ссылки если она есть
            r_timeout = 0
            max_waiting_time_for_pb = 1200

            def pb_timeout_checker(pb_response):
                #  Your request is waiting to be processed...our system has temporarily reached full capacity and the wait time can be much longer than usual.
                server_is_overloaded = 'our system has temporarily reached full capacity and the wait time can be much longer than usual.'
                max_waiting_time_for_pb = 600
                if server_is_overloaded in pb_response:
                    max_waiting_time_for_pb = 1100
                print('Max waiting time: {}'.format(max_waiting_time_for_pb))

                return max_waiting_time_for_pb


            if SEARCH_SPECIFIC_PRIMER == []:
                q = 'Making primers specific to your PCR template.'
                if q in r.text:
                    pb_link_checker = 'checked' # Нужно для добавления ссылки если она есть
                    pb_specific_link = '{}'.format(headers_pb['NCBI-RCGI-RetryURL'])
                    primers_search_dict['pb_specific_link'] = pb_specific_link
                    print('\nPrimer_BLAST link for {}:'.format(exon_number_for_thread))
                    #print(pb_specific_link)
                    #r2 = requests.get(pb_specific_link)
                    r2 = cyc_get_req_for_pb(pb_specific_link)
                    r_timeout = 0
                    count_of_requests = 1
                    while q in r2.text and r_timeout < max_waiting_time_for_pb:
                        max_waiting_time_for_pb = pb_timeout_checker(r.text)
                        time.sleep(15)
                        r_timeout += 15
                        #r2 = requests.get(pb_specific_link)
                        r2 = cyc_get_req_for_pb(pb_specific_link)
                        count_of_requests += 1
                        print('\nRequest {} for {}'.format(count_of_requests, exon_number_for_thread))
                        print('Time has pass - {} sec. Max waiting time - {} sec'.format(r_timeout, max_waiting_time_for_pb))
                        if r_timeout >= 240:
                            print(pb_specific_link)
                else:
                    r2 = r
            elif SEARCH_SPECIFIC_PRIMER == ['checked']:
                def pb_response_status_f(pb_response):
                    pb_response_status = ''
                    highly_similar_seq = 'Your PCR template is highly similar to the following sequence(s) from the search database. To increase the chance' #Your PCR template is highly similar to the following sequence(s) from the search database. To increase the chance of finding specific primers, please review the list below and select all sequences (within the given sequence ranges) that are intended or allowed targets.
                    primers_design_in_progress = 'Making primers specific to your PCR template.'
                    primers_is_ready = 'Detailed primer reports'

                    if primers_is_ready in pb_response.text:
                        pb_response_status = 'primers_is_ready'
                    elif primers_design_in_progress in pb_response.text:
                        pb_response_status = 'primers_design_in_progress'
                    elif 'NCBI-RCGI-RetryURL' in pb_response.headers:
                        pb_response_status = 'pb_redirect'
                    elif highly_similar_seq in pb_response.text:
                        pb_response_status = 'highly_similar_seq'
                    else:
                        pb_response_status = 'another_response'
                    print(pb_response_status)
                    return pb_response_status

                pb_response_status = pb_response_status_f(r)
                waiting_time_is_over = ''
                wrong_params = ''
                while pb_response_status != 'primers_is_ready' and r_timeout < max_waiting_time_for_pb:
                    max_waiting_time_for_pb = pb_timeout_checker(r.text)
                    if pb_response_status == 'primers_design_in_progress':
                        time.sleep(5)
                        pb_link_checker = 'checked'
                        pb_specific_link = '{}'.format(r.headers['NCBI-RCGI-RetryURL'])
                        #primers_search_dict['pb_specific_link'] = pb_specific_link
                        print(pb_specific_link)
                        #r = requests.get(pb_specific_link)
                        #print('Start primers_design_in_progress_for_{}'.format(element))
                        r = cyc_get_req_for_pb(pb_specific_link)
                        pb_response_status = pb_response_status_f(r)
                        #print('End primers_design_in_progress_for_{}'.format(element))
                    elif pb_response_status == 'pb_redirect':
                        pb_link_checker = 'checked' # Нужно для добавления ссылки если она есть
                        pb_specific_link = '{}'.format(r.headers['NCBI-RCGI-RetryURL'])
                        #primers_search_dict['pb_specific_link'] = pb_specific_link
                        print(pb_specific_link)
                        #r = requests.get(pb_specific_link)
                        #print('Start pb_redirect')
                        r = cyc_get_req_for_pb(pb_specific_link)
                        pb_response_status = pb_response_status_f(r)
                        #print('End pb_redirect')
                    elif pb_response_status == 'highly_similar_seq':
                        similar_genes = BeautifulSoup(r.text, 'html.parser')
                        similar_genes_items = similar_genes.find_all('table')
                        #print(similar_genes_items)

                        all_gene_values = re.findall(r'ref\|\w*_\w*\.\w*\|\?\w*\?\w*', str(similar_genes_items))
                        #print('\nall_gene_values')
                        #print(all_gene_values)
                        all_id_inputs = re.findall(r'<tr><td.*</td></tr>', str(similar_genes_items))
                        #print('\nall_id_inputs')
                        #print(all_id_inputs)
                        all_id_inputs = (all_id_inputs[0].replace('<tr><td class="c0">', '')).replace('</td></tr>', '')
                        #print('\nall_id_inputs')
                        #print(all_id_inputs)
                        all_id_inputs = all_id_inputs.split('<input ')
                        #print('\nall_id_inputs')
                        #print(all_id_inputs)
                        all_id_inputs_splitted = []
                        for id_input in all_id_inputs:
                            if id_input != '':
                                all_id_inputs_splitted.append(id_input)
                        #print('\nall_id_inputs_splitted')
                        #print(all_id_inputs_splitted)

                        ###
                        similar_genes_dict = {}
                        ###

                        match = 0
                        for index, gene_in_list in enumerate(all_id_inputs_splitted):
                            gene = re.search(r'target\d?="new_entrez\d?">\w*-*\w*</a>', str(gene_in_list))
                            gene = gene.group(0)
                            print(gene)
                            gene = re.sub(r'target\d?="new_entrez\d?">', '', gene)
                            gene = re.sub(r'</a>', '', gene)
                            #print(gene)
                            similarity = re.search(r'\d*%', str(gene_in_list))
                            similarity = similarity.group(0)
                            similarity = int(similarity.replace('%', ''))
                            #print(similarity)
                            similar_gene_checkbox_id = re.search(r'id="\d*', str(gene_in_list))
                            similar_gene_checkbox_id = similar_gene_checkbox_id.group(0)
                            similar_gene_checkbox_id = re.sub(r'id="', '', similar_gene_checkbox_id)
                            #print(similar_gene_checkbox_id)
                            similar_gene_value = re.search(r'ref\|\w*_\w*\.\w*\|\?\w*\?\w*', str(gene_in_list))
                            similar_gene_value = similar_gene_value.group(0)
                            similar_gene_value = re.sub(r'"/>', '', similar_gene_value)
                            #print(similar_gene_value)
                            match = match + index
                            index_of_match = 'match_{}'.format(match)
                            #print(index_of_match)
                            similar_genes_dict[index_of_match]={}
                            similar_genes_dict[index_of_match]['similar_gene_checkbox_id']=similar_gene_checkbox_id
                            similar_genes_dict[index_of_match]['gene']=gene
                            similar_genes_dict[index_of_match]['similarity']=similarity
                            similar_genes_dict[index_of_match]['value']=similar_gene_value

                        primers_search_dict['similar_genes_dict'] = similar_genes_dict
                        #print('\nsimilar_genes_dict')
                        #print(similar_genes_dict)

                        all_inputs = BeautifulSoup(r.text, 'html.parser')
                        all_inputs_item = all_inputs.find_all('input')
                        all_inputs_dict = {}
                        first_user_rid = BeautifulSoup(r.text, 'html.parser')
                        first_user_rid_item = first_user_rid.find_all('input')
                        #print(first_user_rid_item)
                        for item in first_user_rid_item:
                            try:
                                all_inputs_dict[item['name']]=item['value']
                            except:
                                print('No_name')

                        similarities_of_target_gene = []
                        for gene_checkbox_id in similar_genes_dict:
                            similarities_of_target_gene.append(similar_genes_dict[gene_checkbox_id]['similarity'])
                        highest_similarity_of_target_gene = max(similarities_of_target_gene)
                        #print('\nsimilarities_of_target_gene')
                        #print(similarities_of_target_gene)
                        #print('\nhighest_similarity_of_target_gene')
                        #print(highest_similarity_of_target_gene)
                        #print('\n')

                        choosen_gene_checkbox_id = ''
                        for gene_checkbox_id in similar_genes_dict:
                            if similar_genes_dict[gene_checkbox_id]['gene'].upper() == gene_name and similar_genes_dict[gene_checkbox_id]['similarity'] == highest_similarity_of_target_gene:
                                choosen_gene_checkbox_id = gene_checkbox_id
                                break
                            elif similar_genes_dict[gene_checkbox_id]['gene'].upper() == gene_name:
                                choosen_gene_checkbox_id = gene_checkbox_id
                                break
                            elif similar_genes_dict[gene_checkbox_id]['similarity'] == highest_similarity_of_target_gene:
                                choosen_gene_checkbox_id = gene_checkbox_id
                                break
                            elif similar_genes_dict[gene_checkbox_id]['gene'].upper() != gene_name:
                                choosen_gene_checkbox_id = gene_checkbox_id
                                break
                        all_inputs_dict['USER_SEQLOC']=similar_genes_dict[choosen_gene_checkbox_id]['value']
                        #print(similar_genes_dict[choosen_gene_checkbox_id]['value'])
                        primers_search_dict['gene_checkbox_id'] = gene_checkbox_id

                        url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
                        #r = requests.post(url, data=all_inputs_dict)
                        #r = cycreq(requests.post(url, data=all_inputs_dict, timeout=None))
                        #print('Start highly_similar_seq')
                        r = cyc_post_req(url, all_inputs_dict)
                        #print(r.text)
                        #print(r.headers)
                        pb_response_status = pb_response_status_f(r)
                        #print('End highly_similar_seq')
                    else:
                        #print('ELSE')
                        break
                    time.sleep(1)
                    r_timeout += 1
                    #print(r_timeout)
                else:
                    if r_timeout >= max_waiting_time_for_pb:
                        waiting_time_is_over = 'waiting_time_is_over'
                        print(waiting_time_is_over)
                r2 = r
                #print('CHECK_R_PRIMERS_PAGE_2')
                primers_search_dict['waiting_time_is_over'] = waiting_time_is_over
                primers_search_dict['wrong_params'] = wrong_params

            primers_search_dict['pb_link_checker'] = pb_link_checker
            try:
                primers_search_dict['pb_specific_link'] = pb_specific_link
            except:
                print('No_pb_specific_link')
            #Тут нужно грамотно разместить слип чтобы программа могла делать проверки
            #готовности праймеров и в случае не готовности делать слип и потом снова
            #повторять запрос
            #После того как q не будет в р2.текст, программа будет идти дальше
            html_primers = r2.text
            primers_search_dict['html_primers']=html_primers

            #_______________________________________________________________________________
            return primers_search_dict


        #Первый запуск функции поиска праймеров
        request_status = 'first_request'
        primers_search_dict = primers_search(element, result_dict, statuses_dict, gp_request_id, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, info_dict, request_status, auto_distances)
        primers_are_specific = 'Primer pairs are specific to input template as no other targets were found in selected database'
        non_specific_primers_message = 'Primers may <b>not</b> be specific to the input PCR template'
        non_specific_primers_status = False
        step = 100
        html_primers = primers_search_dict['html_primers']
        #Далее происходит повторный запуск функции поиска праймеров в сулчае если не найдены специичные с первой попытки
        if auto_distances == ['checked'] and primers_are_specific not in html_primers and SEARCH_SPECIFIC_PRIMER == ['checked']:
            get_exon_only = primers_search_dict['get_exon_only']
            max_distance_for_search = (int(PRIMER_PRODUCT_MAX_PB) - len(get_exon_only)) / 2
            max_requests = math.floor(max_distance_for_search / step)
            print('max_distance_for_search')
            print(max_distance_for_search)
            #max_distance_for_search = 300
            test_html_primers = ''
            request_number = 1
            while primers_are_specific not in test_html_primers and F_SEARCH_DISTANCE_PB <= max_distance_for_search:
                request_status = 'additional primers design (try {} / {})'.format(request_number, max_requests)
                request_number += 1
                test_primers_search_dict = primers_search(element, result_dict, statuses_dict, gp_request_id, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, info_dict, request_status, auto_distances)
                test_html_primers = test_primers_search_dict['html_primers'] #Для проверки нахождения нужной строки в ответе, чтобы не сохранялся не тестовый html_primers
                if primers_are_specific in test_html_primers:
                    primers_search_dict = test_primers_search_dict
                    first_F_SEARCH_DISTANCE_PB = F_SEARCH_DISTANCE_PB
                    first_R_SEARCH_DISTANCE_PB = R_SEARCH_DISTANCE_PB
                    break
                elif non_specific_primers_message in test_html_primers and non_specific_primers_status == False:
                    non_specific_primers_status = True #Для того чтобы если есть только неспецифичные праймеры он сохранял самый короткий найденный ампликон
                    primers_search_dict = test_primers_search_dict
                    first_F_SEARCH_DISTANCE_PB = F_SEARCH_DISTANCE_PB
                    first_R_SEARCH_DISTANCE_PB = R_SEARCH_DISTANCE_PB
                F_SEARCH_DISTANCE_PB += step
                R_SEARCH_DISTANCE_PB += step

        html_primers = primers_search_dict['html_primers']
        try:
            info_dict[element]['search_distances']['F_SEARCH_DISTANCE_PB'] = first_F_SEARCH_DISTANCE_PB
            info_dict[element]['search_distances']['R_SEARCH_DISTANCE_PB'] = first_R_SEARCH_DISTANCE_PB
        except:
            pass
        print_dis = primers_search_dict['print_dis']
        print_exon_number = primers_search_dict['print_exon_number']
        exon_sequence_for_users = primers_search_dict['exon_sequence_for_users']
        number_of_exon = primers_search_dict['number_of_exon']
        begin_of_search_seq = primers_search_dict['begin_of_search_seq']
        end_of_search_seq = primers_search_dict['end_of_search_seq']
        pb_link_checker = primers_search_dict['pb_link_checker']
        if 'similar_genes_dict' in primers_search_dict:
            similar_genes_dict = primers_search_dict['similar_genes_dict']
        if 'pb_specific_link' in primers_search_dict:
            pb_specific_link = primers_search_dict['pb_specific_link']
        if 'gene_checkbox_id' in primers_search_dict:
            gene_checkbox_id = primers_search_dict['gene_checkbox_id']
        if 'waiting_time_is_over' in primers_search_dict:
            waiting_time_is_over = primers_search_dict['waiting_time_is_over']
        if 'wrong_params' in primers_search_dict:
            wrong_params = primers_search_dict['wrong_params']

        state_statuses(statuses_dict, element, gp_request_id, 'request_from_pb_received')
        primers_for_status_were_not_found = '<p class="info">No primers were found'
        primers_for_status_were_found = '<th>Template strand</th><th>Length</th><th>Start</th><th>Stop</th><th>Tm</th><th>GC%</th>'
        # Full string: <tr><th></th><th>Sequence (5'->3')</th><th>Template strand</th><th>Length</th><th>Start</th><th>Stop</th><th>Tm</th><th>GC%</th><th>Self complementarity</th><th>Self 3' complementarity</th></tr>
        if primers_for_status_were_not_found in html_primers:
            state = 'primers_were_not_found'
        elif primers_for_status_were_found in html_primers:
            state = 'primers_were_found'
        else:
            state = 'primers_were_not_found'
        state_statuses(statuses_dict, element, gp_request_id, state)

        #Часть по парсингу и получению только нужных строк
        soup = BeautifulSoup(html_primers, 'html.parser')
        items = soup.find_all('div', class_='prPairInfo')

        #Тут создается словарь в который будут включены все элементы всех праймеров
        all_primers_dicts = {}

        #Часть по принту того что намного выше (тут возможны затупы, менял на абум)
        if F_SEARCH_DISTANCE_PB != 0 and R_SEARCH_DISTANCE_PB != 0 and FIVE_SAVE_EXON_DISTANCE_PB != 0 and THREE_SAVE_EXON_DISTANCE_PB != 0:
            print_print_dis = print_dis
        seq_string = '{}\n{}\n{}\n{}'.format(print_exon_number, print_print_dis, exon_sequence_for_users, '')

        all_primer_strings = []
        all_primer_strings.append(seq_string)
        print_everything = '\n' #Здесь я создаю строку с переносом строки, ниже я буду делать так:
            # print_everything = print_everything.join(all_primer_strings) чтобы при применении многопоточности
            # у меня принт делался с помощью одной переменной, а не цикла фор, потому что тогда разные
            # потоки могут пересекаться при принте и смешивать его в кашу

        #Тут я полученный контейнер 'div', class_='prPairInfo' (в контейнере сразу все праймеры и их параметры)
        # разбиваю по индивидуальным спискам сначала через all, затем каждый элемент в all вызываю циклом
        # и с помощью get_text избавляюсь от тегов, оставляя только текст, все это уходит в новый список all_td_text
        for i in items:
            number_of_primerpair = i.find('h2').get_text()
            all_td = []
            all_td_text = []
            all_td.append(i.find_all('td'))
            for i in all_td[0]:
                all_td_text.append(i.get_text())

        #Тут уже просто каждый элемент спика all_td_text я добавляю в словарь с соответствующим ключом как для реверса,
        #так и для форварда и их параметров и общей длины ампликона
            primers_dict = {}
            primers_dict = {
                'direct_F': 'Forward primer',
                'seq_F': all_td_text[0],
                'template_strand_F': all_td_text[1],
                'len_F': all_td_text[2],
                'start_F': all_td_text[3],
                'stop_F': all_td_text[4],
                'tm_F': all_td_text[5],
                'gc_F': all_td_text[6],
                'sc_F': all_td_text[7],
                's3c_F': all_td_text[8],
                'direct_R': 'Reverse primer',
                'seq_R': all_td_text[9],
                'template_strand_R': all_td_text[10],
                'len_R': all_td_text[11],
                'start_R': all_td_text[12],
                'stop_R': all_td_text[13],
                'tm_R': all_td_text[14],
                'gc_R': all_td_text[15],
                'sc_R': all_td_text[16],
                's3c_R': all_td_text[17],
                'len_amp': all_td_text[27]
            }

            number_of_primerpair = '{} '.format(number_of_primerpair.replace(' ', '_'))
            all_primers_dicts[number_of_primerpair]=primers_dict

            #Далее нужно представить результаты подборки праймеров в удобном формате
            titles = '{:20}{:30}{:20}{:10}{:10}{:10}{:10}{:10}{:25}{:25}'.format(number_of_primerpair, 'Sequence (5\'->3\')', 'Template strand', 'Length', 'Start', 'Stop', 'Tm', 'GC%', 'Self complementarity', 'Self 3\' complementarity')
            f_primer = '{:20}{:30}{:20}{:10}{:10}{:10}{:10}{:10}{:25}{:25}'.format(primers_dict['direct_F'], primers_dict['seq_F'], primers_dict['template_strand_F'], primers_dict['len_F'], primers_dict['start_F'], primers_dict['stop_F'], primers_dict['tm_F'], primers_dict['gc_F'], primers_dict['sc_F'], primers_dict['s3c_F'])
            r_primer = '{:20}{:30}{:20}{:10}{:10}{:10}{:10}{:10}{:25}{:25}'.format(primers_dict['direct_R'], primers_dict['seq_R'], primers_dict['template_strand_R'], primers_dict['len_R'], primers_dict['start_R'], primers_dict['stop_R'], primers_dict['tm_R'], primers_dict['gc_R'], primers_dict['sc_R'], primers_dict['s3c_R'])
            len_amp = '{:20}{:30}\n'.format('Product length', primers_dict['len_amp'])

            primer_string = '{}\n{}\n{}\n{}'.format(titles, f_primer, r_primer, len_amp)
            all_primer_strings.append(primer_string)

        #Часть по полному принту вывода (сиквенс + таблица праймеров)
        print_everything = print_everything.join(all_primer_strings)
        print(print_everything)
        #Присваивание выше созданному дикту пары ключ-значение таким вот способом(через знак =).
        result_dict['{}'.format(number_of_exon)]=all_primers_dicts
        #Тоже самое, только на уровень ниже, добавляю последовательность экзона в дикт.
        result_dict['{}'.format(number_of_exon)]['exon_sequence']=exon_sequence_for_users
        #Тоже самое, добавляются позиции начала и конца последовательности для поиска праймеров, необходимо это для отбора по SNP ниже.
        result_dict['{}'.format(number_of_exon)]['seq_search_positions']={'seq_search_start': begin_of_search_seq, 'seq_search_end': end_of_search_seq}
        #Добавляется ссылка на праймербласт
        if pb_link_checker == 'checked' and SHOW_PB_LINK == ['checked']:
            result_dict['{}'.format(number_of_exon)]['primer_blast_link']=pb_specific_link
            #print(pb_link_checker)
            #print(pb_specific_link)
        #Добавление чекера на неспицифику праймеров
        non_specific_primers_message = 'Primers may <b>not</b> be specific to the input PCR template'#тут дальше продолжнение
        #Primers may not be specific to the input PCR template as targets were found in selected database:Refseq mRNA (Organism limited to Homo sapiens)
        if non_specific_primers_message in html_primers:
            result_dict['{}'.format(number_of_exon)]['specific_checker']=['non_specific_primers']
        #Нужно для отслеживания того что все праймеры подобранные ПБ были выброшены из-за полиморфизмов и отображения пользователю.
        if 'Primer_pair_1 ' in result_dict['{}'.format(number_of_exon)]:
            result_dict['{}'.format(number_of_exon)]['primers_checker']=['checked']
        #Злесь в случае если праймер бласт нашел совпадение с другим (не нашим геном), то в дикс будет записано сообщения об этом, для уведомления пользователя.
        if SEARCH_SPECIFIC_PRIMER == ['checked']:
            if similar_genes_dict[gene_checkbox_id]['gene'].upper() != gene_name:
                result_dict['{}'.format(number_of_exon)]['similar_gene']=similar_genes_dict[gene_checkbox_id]['gene'].upper()
        #Если время подбора праймеров больше максимальногов времени, то процесс будет прерван и сделана об этом запись в дикт.
        if SEARCH_SPECIFIC_PRIMER == ['checked']:
            if waiting_time_is_over == 'waiting_time_is_over':
                result_dict['{}'.format(number_of_exon)]['waiting_time_is_over']=['waiting_time_is_over']
        #Если поиск праймеров не удовлетворяет заданным параметрам об этом будет записано.
        if SEARCH_SPECIFIC_PRIMER == ['checked']:
            if wrong_params == 'wrong_params':
                result_dict['{}'.format(number_of_exon)]['wrong_params']=['wrong_params']


    #Формирование и запуск разных потоков
    threads = []
    for element in taken_exons_id:
        t = threading.Thread(target=get_get_primers, args=[element, result_dict, statuses_dict, gp_request_id, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, FIVE_SAVE_EXON_DISTANCE_PB, THREE_SAVE_EXON_DISTANCE_PB, info_dict, auto_distances])
        t.start()
        threads.append(t)
        time.sleep(1.1) #Здесь слип необходим потому что между запуском потока нужны интервалы, иначе библиотека Requests жаулется на слишком большок количество запросов за единицу времени, поэтому искусственно замедляем.
        #print('Поток запущен.')
    for t in threads:
        t.join()

    exons_result_dict = []
    for ee in result_dict:
        exons_result_dict.append(ee)


    #Отсев парймеров поп полиморфизмам
    if NO_SNP == ['checked']:
        MAF = float(MAX_MAF)
        #Здесь будет происходить отбор праймероа на отсуствие полиморфизмов в месте их отжига на темплейте.
        #Запросы для форвард праймеров.
        result_dict_deepcopy = copy.deepcopy(result_dict)
        print(result_dict)

        not_primer_pairs = ['exon_sequence', 'seq_search_positions', 'primer_blast_link', 'specific_checker', 'waiting_time_is_over', 'primers_checker', 'similar_gene', 'wrong_params']

        for exon in result_dict:
            #Если праймеры есть, будет создан статус об их проверке на полиморфизмы
            if 'Primer_pair_1 ' in result_dict[exon]:
                state = 'polymorphisms_checking'
                state_statuses(statuses_dict, exon, gp_request_id, state)
            print(exon)
            for primer_pair in result_dict[exon]:
                contain_polymorphisms = False
                if primer_pair not in not_primer_pairs:
                    print('\n---Forward primers---\n')
                    print('\n{}({})'.format(primer_pair, exon))

                    if strand == '1':
                        genome_primer_F_position = '{}:{}-{}'.format(chromosome, result_dict[exon]['seq_search_positions']['seq_search_start'] + (result_dict[exon]['exon_sequence']).lower().find((result_dict[exon][primer_pair]['seq_F']).lower()), result_dict[exon]['seq_search_positions']['seq_search_start'] + (result_dict[exon]['exon_sequence']).lower().find((result_dict[exon][primer_pair]['seq_F']).lower()) + len(result_dict[exon][primer_pair]['seq_F']) - 1)
                    elif strand == '-1':
                        genome_primer_F_position = '{}:{}-{}'.format(chromosome, result_dict[exon]['seq_search_positions']['seq_search_end'] - (result_dict[exon]['exon_sequence']).lower().find((result_dict[exon][primer_pair]['seq_F']).lower()) - len(result_dict[exon][primer_pair]['seq_F']) + 1, result_dict[exon]['seq_search_positions']['seq_search_end'] - (result_dict[exon]['exon_sequence']).lower().find((result_dict[exon][primer_pair]['seq_F']).lower()))
                    print(genome_primer_F_position)

                    #Заопрос возвращающий всю инфу по данному диапазону, включая все варианты.
                    ext = "/overlap/region/{}/{}?feature=variation;".format(SPECIES, genome_primer_F_position)
                    headers={"Content-Type" : "application/json"}
                    #dict_of_variants_info = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                    #dict_of_variants_info = cycreq(requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=None))
                    dict_of_variants_info = cyc_get_req_for_ens(server, ext, headers)

                    if not dict_of_variants_info.ok:
                      dict_of_variants_info.raise_for_status()
                      sys.exit()
                    dict_of_variants_info = dict_of_variants_info.json()
                    #print(dict_of_variants_info)

                    list_of_variants = []
                    for variant in dict_of_variants_info:
                        list_of_variants.append(variant['id'])

                    dict_of_variants = {}
                    dict_of_variants["ids"]=list_of_variants
                    dict_of_variants = str(dict_of_variants).replace('\'', '\"')
                    print(dict_of_variants)


                    #Запрос собирающий инфу по найденным выше вариантам, включая MAF.
                    ext = "/variation/{}".format(SPECIES)
                    headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
                    #variation_info = requests.post(server+ext, headers=headers, data=dict_of_variants)
                    #variation_info = cycreq(requests.post(server+ext, headers=headers, data=dict_of_variants, timeout=None))
                    variation_info = cyc_post_req_for_ens(server, ext, headers, dict_of_variants)

                    if not variation_info.ok:
                      variation_info.raise_for_status()
                      sys.exit()
                    variation_info = variation_info.json()

                    variants_maf = {}
                    for variant in variation_info:
                        variants_maf[variant]={'MAF': variation_info[variant]['MAF'], 'start': variation_info[variant]['mappings'][0]['start'], 'end': variation_info[variant]['mappings'][0]['end']}
                    print(variants_maf)

                    for variant_maf in variants_maf:
                        print(variants_maf[variant_maf]['MAF'])
                        if variants_maf[variant_maf]['MAF'] != None:
                            print('2')
                            if variants_maf[variant_maf]['MAF'] > MAF:
                                print('3')
                                del result_dict_deepcopy[exon][primer_pair]
                                print("Отбор по SNP выполнен")
                                contain_polymorphisms = True
                                break

                    if contain_polymorphisms != True:
                        print('\n---Reverse primers---\n')
                        print('\n{}({})'.format(primer_pair, exon))

                        #Запросы для реверс праймеров.
                        if strand == '1':
                            genome_primer_R_position = '{}:{}-{}'.format(chromosome, result_dict_deepcopy[exon]['seq_search_positions']['seq_search_start'] + (result_dict_deepcopy[exon]['exon_sequence']).lower().find(get_reverse_complement_seq(result_dict_deepcopy[exon][primer_pair]['seq_R']).lower()), result_dict_deepcopy[exon]['seq_search_positions']['seq_search_start'] + (result_dict_deepcopy[exon]['exon_sequence']).lower().find(get_reverse_complement_seq(result_dict_deepcopy[exon][primer_pair]['seq_R']).lower()) + len(result_dict_deepcopy[exon][primer_pair]['seq_F']) - 1)
                        elif strand == '-1':
                            genome_primer_R_position = '{}:{}-{}'.format(chromosome, result_dict_deepcopy[exon]['seq_search_positions']['seq_search_end'] - (result_dict_deepcopy[exon]['exon_sequence']).lower().find(get_reverse_complement_seq(result_dict_deepcopy[exon][primer_pair]['seq_R']).lower()) - len(result_dict_deepcopy[exon][primer_pair]['seq_F']) + 1, result_dict_deepcopy[exon]['seq_search_positions']['seq_search_end'] - (result_dict_deepcopy[exon]['exon_sequence']).lower().find(get_reverse_complement_seq(result_dict_deepcopy[exon][primer_pair]['seq_R']).lower()))
                        print(genome_primer_R_position)


                        #Заопрос возвращающий всю инфу по данному диапазону, включая все варианты.
                        ext = "/overlap/region/{}/{}?feature=variation;".format(SPECIES, genome_primer_R_position)
                        headers={"Content-Type" : "application/json"}
                        #dict_of_variants_info = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
                        #dict_of_variants_info = cycreq(requests.get(server+ext, headers={ "Content-Type" : "application/json"}, timeout=None))
                        dict_of_variants_info = cyc_get_req_for_ens(server, ext, headers)

                        if not dict_of_variants_info.ok:
                          dict_of_variants_info.raise_for_status()
                          sys.exit()
                        dict_of_variants_info = dict_of_variants_info.json()

                        list_of_variants = []
                        for variant in dict_of_variants_info:
                            list_of_variants.append(variant['id'])

                        dict_of_variants = {}
                        dict_of_variants["ids"]=list_of_variants
                        dict_of_variants = str(dict_of_variants).replace('\'', '\"')
                        print(dict_of_variants)


                        #Запрос собирающий инфу по найденным выше вариантам, включая MAF.
                        ext = "/variation/{}".format(SPECIES)
                        headers={ "Content-Type" : "application/json", "Accept" : "application/json"}
                        #variation_info = requests.post(server+ext, headers=headers, data=dict_of_variants)
                        #variation_info = cycreq(requests.post(server+ext, headers=headers, data=dict_of_variants, timeout=None))
                        variation_info = cyc_post_req_for_ens(server, ext, headers, dict_of_variants)

                        if not variation_info.ok:
                          variation_info.raise_for_status()
                          sys.exit()
                        variation_info = variation_info.json()

                        variants_maf = {}
                        for variant in variation_info:
                            variants_maf[variant]={'MAF': variation_info[variant]['MAF'], 'start': variation_info[variant]['mappings'][0]['start'], 'end': variation_info[variant]['mappings'][0]['end']}
                        print(variants_maf)

                        for variant_maf in variants_maf:
                            print(variants_maf[variant_maf]['MAF'])
                            if variants_maf[variant_maf]['MAF'] != None:
                                print('2')
                                if variants_maf[variant_maf]['MAF'] > MAF:
                                    print('3')
                                    del result_dict_deepcopy[exon][primer_pair]
                                    print("Отбор по SNP выполнен")
                                    break

            for number_of_primer_pair in range(1, 11):
                if 'Primer_pair_{} '.format(number_of_primer_pair) in result_dict_deepcopy[exon]:
                    state = 'primers_do_not_contain_polymorphisms'
                    state_statuses(statuses_dict, exon, gp_request_id, state)
                    break
                #else:
                elif 'Primer_pair_{} '.format(number_of_primer_pair) not in result_dict_deepcopy[exon] and 'primers_checker' in result_dict_deepcopy[exon]:
                    state = 'primers_contain_polymorphisms'
                    state_statuses(statuses_dict, exon, gp_request_id, state)

        #Все тоже самое, только для реверс праймеров.
        result_dict_deepcopy_2 = copy.deepcopy(result_dict_deepcopy)

        #for position, element in enumerate(taken_exons_list):
        #Теперь нужно заменить названия ключей в словаре чтобы пары праймеров после удаления некоторых из них шли по порядку.
        result_dict_deepcopy_3 = copy.deepcopy(result_dict_deepcopy_2)
        for exon in result_dict_deepcopy_2:
            for position, primer_pair in enumerate(result_dict_deepcopy_2[exon]):
                if primer_pair not in not_primer_pairs:
                    result_dict_deepcopy_3[exon]['Primer_pair_{} '.format(position + 1)]=result_dict_deepcopy_3[exon].pop(primer_pair)

        result_dict = result_dict_deepcopy_3



    def ordered_exons(result_dict):
        # Функция упорядочивания экзонов в словаре
        ordered_dict = {}
        ordered_list_of_exons = []
        for exon in result_dict:
            exon_number = re.search(r'Exons?_\d*', exon)
            exon_number = exon_number.group(0)
            exon_number = exon_number.replace('s', '')
            ordered_list_of_exons.append(int(exon_number.replace('Exon_', '')))
        ordered_list_of_exons = sorted(ordered_list_of_exons)
        for exon_number in ordered_list_of_exons:
            for exon in result_dict:
                exon_number_in_result_dict = re.search(r'Exons?_\d*', exon)
                exon_number_in_result_dict = exon_number_in_result_dict.group(0)
                exon_number_in_result_dict = exon_number_in_result_dict.replace('s', '')
                exon_number_in_result_dict = exon_number_in_result_dict.replace('Exon_', '')
                if int(exon_number_in_result_dict) == int(exon_number):
                    ordered_dict[exon]=result_dict[exon]

        return ordered_dict

    result_dict = ordered_exons(result_dict)


    result_dict_deepcopy_4 = copy.deepcopy(result_dict)
    scheme_dict = result_dict_deepcopy_4
    def seq_dissection(result_dict, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, scheme_dict, info_dict):
        # Функция для рассчета позиций различных элементов последовательности для построения схем
        for exon in result_dict:
            print(exon, '\n', result_dict[exon]['exon_sequence'])
            for primer_pair in result_dict[exon]:
                if 'Primer_pair' in primer_pair:
                    pre_primer_last_enrty = result_dict[exon]['exon_sequence'].find(result_dict[exon][primer_pair]['seq_F'].lower())
                    f_primer_first_enrty = result_dict[exon]['exon_sequence'].find(result_dict[exon][primer_pair]['seq_F'].lower()) + 1
                    f_primer_last_enrty = result_dict[exon]['exon_sequence'].find(result_dict[exon][primer_pair]['seq_F'].lower()) + len(result_dict[exon][primer_pair]['seq_F'].lower()) - 1
                    first_gap_first_enrty = result_dict[exon]['exon_sequence'].find(result_dict[exon][primer_pair]['seq_F'].lower()) + len(result_dict[exon][primer_pair]['seq_F'].lower())
                    first_gap_last_enrty = info_dict[exon]['search_distances']['F_SEARCH_DISTANCE_PB']
                    exon_first_enrty = info_dict[exon]['search_distances']['F_SEARCH_DISTANCE_PB']
                    exon_last_enrty = len(result_dict[exon]['exon_sequence']) - info_dict[exon]['search_distances']['R_SEARCH_DISTANCE_PB']
                    second_gap_first_enrty = len(result_dict[exon]['exon_sequence']) - info_dict[exon]['search_distances']['R_SEARCH_DISTANCE_PB']
                    second_gap_last_enrty = result_dict[exon]['exon_sequence'].find(get_reverse_complement_seq(result_dict[exon][primer_pair]['seq_R']).lower())
                    r_primer_first_enrty = result_dict[exon]['exon_sequence'].find(get_reverse_complement_seq(result_dict[exon][primer_pair]['seq_R']).lower())
                    r_primer_last_enrty = r_primer_first_enrty + len(result_dict[exon][primer_pair]['seq_R']) - 1
                    post_primer_first_entry = r_primer_last_enrty + 1
                    post_primer_last_entry = len(result_dict[exon]['exon_sequence'])

                    def percentage(element, full_lenght):
                        percentage = (len(element) * 100) / full_lenght
                        percentage = float('{0:.5f}'.format(percentage))
                        return percentage

                    #scheme_dict[primer_pair] = {}
                    scheme_dict[exon][primer_pair]['pre_primer'] = {}
                    scheme_dict[exon][primer_pair]['pre_primer']['seq']=result_dict[exon]['exon_sequence'][0:pre_primer_last_enrty]
                    scheme_dict[exon][primer_pair]['pre_primer']['percentage']=percentage(result_dict[exon]['exon_sequence'][0:pre_primer_last_enrty], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['pre_primer']['partial_percentage']=percentage(result_dict[exon]['exon_sequence'][0:pre_primer_last_enrty], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['pre_primer']['len']=len(result_dict[exon]['exon_sequence'][0:pre_primer_last_enrty])

                    scheme_dict[exon][primer_pair]['f_primer'] = {}
                    scheme_dict[exon][primer_pair]['f_primer']['seq']=result_dict[exon][primer_pair]['seq_F']
                    scheme_dict[exon][primer_pair]['f_primer']['percentage']=percentage(result_dict[exon][primer_pair]['seq_F'], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['f_primer']['partial_percentage']=percentage(result_dict[exon][primer_pair]['seq_F'], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['f_primer']['len']=len(result_dict[exon][primer_pair]['seq_F'])

                    scheme_dict[exon][primer_pair]['first_gap'] = {}
                    scheme_dict[exon][primer_pair]['first_gap']['seq']=result_dict[exon]['exon_sequence'][first_gap_first_enrty:first_gap_last_enrty]
                    scheme_dict[exon][primer_pair]['first_gap']['percentage']=percentage(result_dict[exon]['exon_sequence'][first_gap_first_enrty:first_gap_last_enrty], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['first_gap']['partial_percentage']=percentage(result_dict[exon]['exon_sequence'][first_gap_first_enrty:first_gap_last_enrty], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['first_gap']['len']=len(result_dict[exon]['exon_sequence'][first_gap_first_enrty:first_gap_last_enrty])

                    scheme_dict[exon][primer_pair]['exon'] = {}
                    scheme_dict[exon][primer_pair]['exon']['seq']=result_dict[exon]['exon_sequence'][exon_first_enrty:exon_last_enrty]
                    scheme_dict[exon][primer_pair]['exon']['percentage']=percentage(result_dict[exon]['exon_sequence'][exon_first_enrty:exon_last_enrty], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['exon']['partial_percentage']=percentage(result_dict[exon]['exon_sequence'][exon_first_enrty:exon_last_enrty], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['exon']['len']=len(result_dict[exon]['exon_sequence'][exon_first_enrty:exon_last_enrty])

                    scheme_dict[exon][primer_pair]['second_gap'] = {}
                    scheme_dict[exon][primer_pair]['second_gap']['seq']=result_dict[exon]['exon_sequence'][second_gap_first_enrty:second_gap_last_enrty]
                    scheme_dict[exon][primer_pair]['second_gap']['percentage']=percentage(result_dict[exon]['exon_sequence'][second_gap_first_enrty:second_gap_last_enrty], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['second_gap']['partial_percentage']=percentage(result_dict[exon]['exon_sequence'][second_gap_first_enrty:second_gap_last_enrty], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['second_gap']['len']=len(result_dict[exon]['exon_sequence'][second_gap_first_enrty:second_gap_last_enrty])

                    scheme_dict[exon][primer_pair]['r_primer'] = {}
                    scheme_dict[exon][primer_pair]['r_primer']['seq']=result_dict[exon][primer_pair]['seq_R']
                    scheme_dict[exon][primer_pair]['r_primer']['percentage']=percentage(result_dict[exon][primer_pair]['seq_R'], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['r_primer']['partial_percentage']=percentage(result_dict[exon][primer_pair]['seq_R'], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['r_primer']['len']=len(result_dict[exon][primer_pair]['seq_R'])

                    scheme_dict[exon][primer_pair]['post_primer'] = {}
                    scheme_dict[exon][primer_pair]['post_primer']['seq']=result_dict[exon]['exon_sequence'][post_primer_first_entry:post_primer_last_entry]
                    scheme_dict[exon][primer_pair]['post_primer']['percentage']=percentage(result_dict[exon]['exon_sequence'][post_primer_first_entry:post_primer_last_entry], len(result_dict[exon]['exon_sequence']))
                    scheme_dict[exon][primer_pair]['post_primer']['partial_percentage']=percentage(result_dict[exon]['exon_sequence'][post_primer_first_entry:post_primer_last_entry], len(result_dict[exon]['exon_sequence'][f_primer_first_enrty:r_primer_last_enrty]) + 2)
                    scheme_dict[exon][primer_pair]['post_primer']['len']=len(result_dict[exon]['exon_sequence'][post_primer_first_entry:post_primer_last_entry])

                    template_for_homopolymers_search = (result_dict[exon]['exon_sequence'][f_primer_last_enrty:r_primer_first_enrty]).upper()
                    def homopolymes_search(template_for_homopolymers_search):
                        dict_with_homopolymers = {}
                        homopolymers_lens = list(range(8, 36))
                        homopolymers_lens = homopolymers_lens[::-1]
                        nucleotides = ['A', 'T', 'G', 'C']
                        for nucleotide in nucleotides:
                            for len in homopolymers_lens:
                                if (nucleotide * len) in template_for_homopolymers_search:
                                    dict_with_homopolymers[nucleotide]=len
                                    break

                        print(dict_with_homopolymers)
                        return dict_with_homopolymers

                    dict_with_homopolymers = homopolymes_search(template_for_homopolymers_search)
                    if dict_with_homopolymers != {}:
                        scheme_dict[exon][primer_pair]['dict_with_homopolymers']=dict_with_homopolymers
                        scheme_dict[exon]['homopolymers']={}
                        scheme_dict[exon]['homopolymers']=['checked']

                    regions_list = ['first_gap', 'exon', 'second_gap']
                    for region in regions_list:
                        for nucleotide in dict_with_homopolymers:
                            full_seq = (scheme_dict[exon][primer_pair][region]['seq']).lower()
                            homopolymers_lenght = (nucleotide * dict_with_homopolymers[nucleotide]).lower()
                            if nucleotide == 'A':
                                seq_with_homopymers = '{}{}{}'.format('<span class="border-bottom border-warning border-3">', nucleotide * dict_with_homopolymers[nucleotide], '</span>').lower()
                            elif nucleotide == 'T':
                                seq_with_homopymers = '{}{}{}'.format('<span class="border-bottom border-info border-3">', nucleotide * dict_with_homopolymers[nucleotide], '</span>').lower()
                            elif nucleotide == 'G':
                                seq_with_homopymers = '{}{}{}'.format('<span class="border-bottom border-danger border-3">', nucleotide * dict_with_homopolymers[nucleotide], '</span>').lower()
                            elif nucleotide == 'C':
                                seq_with_homopymers = '{}{}{}'.format('<span class="border-bottom border-primary border-3">', nucleotide * dict_with_homopolymers[nucleotide], '</span>').lower()
                            scheme_dict[exon][primer_pair][region]['seq']=full_seq.replace(homopolymers_lenght, seq_with_homopymers)

        return scheme_dict

    scheme_dict = seq_dissection(result_dict, F_SEARCH_DISTANCE_PB, R_SEARCH_DISTANCE_PB, result_dict_deepcopy_4, info_dict)
    result_dict = scheme_dict

    result_dict_json_name = 'result_dict_{}'.format(gp_request_id)
    with open('result_dicts_dir/{}.json'.format(result_dict_json_name), 'w') as result_dicts_json:
        json.dump(result_dict, result_dicts_json)
    print(result_dict)

    for exon in result_dict:
        state_statuses(statuses_dict, exon, gp_request_id, 'completed')
    print(statuses_dict)

    return result_dict

#get_primers('11', '1', ['5', '6', '11', '12', '13'], ['ENSE00001324401', 'ENSE00003545357', 'ENSE00003558187', 'ENSE00003519127', 'ENSE00000999053', 'ENSE00003663664', 'ENSE00003525485', 'ENSE00003477895', 'ENSE00000824551', 'ENSE00003505813', 'ENSE00003556991', 'ENSE00003646599', 'ENSE00001101378', 'ENSE00002176251'], {'ENSE00001324401': {'start': 44117099, 'end': 44117402}, 'ENSE00003545357': {'start': 44129233, 'end': 44129798}, 'ENSE00003558187': {'start': 44130744, 'end': 44130833}, 'ENSE00003519127': {'start': 44135735, 'end': 44135851}, 'ENSE00000999053': {'start': 44146339, 'end': 44146534}, 'ENSE00003663664': {'start': 44148366, 'end': 44148505}, 'ENSE00003525485': {'start': 44151595, 'end': 44151688}, 'ENSE00003477895': {'start': 44193161, 'end': 44193292}, 'ENSE00000824551': {'start': 44219379, 'end': 44219568}, 'ENSE00003505813': {'start': 44228343, 'end': 44228509}, 'ENSE00003556991': {'start': 44253903, 'end': 44254046}, 'ENSE00003646599': {'start': 44255665, 'end': 44255793}, 'ENSE00001101378': {'start': 44257843, 'end': 44257925}, 'ENSE00002176251': {'start': 44265699, 'end': 44266482}}, ['checked'], [], 'GRCh37')



def get_clear_seq_for_textarea(seq):
    symbols = [' ', '\n', '\r', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', ',', 'Intron', 'Exon', 'ENSE', '.', '/', '|', '\\', '\'', '\"', ';', ':', '<', '>', '!', '@', '#', '$', '%', '^', '&', '*', '_', '(', ')', '+', '=', '							']
    for element in symbols:
        if element in seq:
            seq = seq.replace(element, '')

    return seq


def get_reverse_complement_seq(seq):
    symbols = [' ', '\n', '\r', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', ',', 'Intron', 'Exon', 'ENSE', '.', '/', '|', '\\', '\'', '\"', ';', ':', '<', '>', '!', '@', '#', '$', '%', '^', '&', '*', '_', '(', ')', '+', '=', '							']
    for element in symbols:
        if element in seq:
            seq = seq.replace(element, '')

    #small_nucleotides = []
    #big_nucleotides = []
    #for position, nucleotide in enumerate(seq):
    #    if nucleotide == nucleotide.lower():
    #        small_nucleotides.append(position)
    #    elif nucleotide == nucleotide.upper():
    #        big_nucleotides.append(position)

    reversed_seq = seq[::-1]
    reverse_complement_seq = ''
    for nucleotide in reversed_seq:
        if nucleotide == 'a':
            reverse_complement_seq = reverse_complement_seq + 't'
        elif nucleotide == 'A':
            reverse_complement_seq = reverse_complement_seq + 'T'
        elif nucleotide == 't':
            reverse_complement_seq = reverse_complement_seq + 'a'
        elif nucleotide == 'T':
            reverse_complement_seq = reverse_complement_seq + 'A'
        elif nucleotide == 'g':
            reverse_complement_seq = reverse_complement_seq + 'c'
        elif nucleotide == 'G':
            reverse_complement_seq = reverse_complement_seq + 'C'
        elif nucleotide == 'c':
            reverse_complement_seq = reverse_complement_seq + 'g'
        elif nucleotide == 'C':
            reverse_complement_seq = reverse_complement_seq + 'G'
        elif nucleotide == 'n':
            reverse_complement_seq = reverse_complement_seq + 'n'
        elif nucleotide == 'N':
            reverse_complement_seq = reverse_complement_seq + 'N'

    #reverse_complement_seq_fixed_positions = ''
    #for position, nucleotide in enumerate(reverse_complement_seq):
    #    if position in small_nucleotides:
    #        reverse_complement_seq_fixed_positions += nucleotide.lower()
    #    elif position in big_nucleotides:
    #        reverse_complement_seq_fixed_positions += nucleotide.upper()

    return reverse_complement_seq


def get_seq_len(seq):
    symbols = [' ', '\n', '\r', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', ',', 'Intron', 'Exon', 'ENSE', '.', '/', '|', '\\', '\'', '\"', ';', ':', '<', '>', '!', '@', '#', '$', '%', '^', '&', '*', '_', '(', ')', '+', '=', '							']
    for element in symbols:
        if element in seq:
            seq = seq.replace(element, '')
    len_seq = ''
    if seq:
        len_seq = len(seq)
    else:
        len_seq = 0

    return len_seq


def get_gc_content_f(seq):
    symbols = [' ', '\n', '\r', '0', '1', '2', '3', '4', '5', '6', '7', '8', '9', '-', ',', 'Intron', 'Exon', 'ENSE', '.', '/', '|', '\\', '\'', '\"', ';', ':', '<', '>', '!', '@', '#', '$', '%', '^', '&', '*', '_', '(', ')', '+', '=', '							']
    for element in symbols:
        if element in seq:
            seq = seq.replace(element, '')

    seq = seq.lower()
    gc_len = 0
    for nucleotide in seq:
        if nucleotide == 'g' or nucleotide == 'c':
            gc_len += 1

    at_len = 0
    for nucleotide in seq:
        if nucleotide == 'a' or nucleotide == 't':
            at_len += 1

    general_len = gc_len + at_len

    if seq == '':
        gc_content = 0
    else:
        gc_content = (gc_len / general_len) * 100
        gc_content = '{:.2f} %'.format(gc_content)

    return gc_content


def primers_design_time_counter(taken_exons, SEARCH_SPECIFIC_PRIMER, NO_SNP, pb_server_status):
    import math
    waiting_time_dict = {}
    taken_exons_count = len(taken_exons)
    primers_design_time = 0
    max_primers_design_time = 0
    if SEARCH_SPECIFIC_PRIMER == ['checked'] and NO_SNP == ['checked']:
        primers_design_time = (int(taken_exons_count) * 30) + 210
    elif SEARCH_SPECIFIC_PRIMER == ['checked'] and NO_SNP == []:
        primers_design_time = (int(taken_exons_count) * 2) + 210
    elif SEARCH_SPECIFIC_PRIMER == [] and NO_SNP == []:
        primers_design_time = int(taken_exons_count) * 2
    counted_primers_design_time = math.ceil(primers_design_time / 60)

    if pb_server_status == 'ok':
        max_primers_design_time = math.ceil((primers_design_time + 421 - 210) / 60)
    elif pb_server_status == 'overloaded':
        max_primers_design_time = math.ceil((primers_design_time + 1201 - 210) / 60)

    waiting_time_dict['primers_design_time']=counted_primers_design_time
    waiting_time_dict['max_primers_design_time']=max_primers_design_time

    return waiting_time_dict


def pb_server_status_checker():
    import requests
    #HGD gene ex. 13 (len = 182)
    seq_for_pb = 'GGAACTGCATGAGTGAGTTCATGGGACTCATCCGAGGTCACTATGAGGCAAAGCAAGGTGGGTTCCTGCCAGGGGGAGGGAGTCTACACAGCACAATGACCCCCCATGGACCTGATGCTGACTGCTTTGAGAAGGCCAGCAAGGTCAAGCTGGCACCTGAGAGGATTGCCGATGGCACCATG'
    url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
    data = {
        'INPUT_SEQUENCE': seq_for_pb,
        'PRIMER5_START': '1',
        'PRIMER5_END': '70',
        'PRIMER3_START': '110',
        'PRIMER3_END': '182'
    }

    #r = requests.post(url, data=data)
    r = cyc_post_req(url, data)

    #Здесь условная конструкция, которая если SEARCH_SPECIFIC_PRIMER = [], говорить о том, что выбран неспецифический подбор праймеров
    #а значит нет смысла брать ссылку из первого запроса на второй запрос, который специфически подбирает праймеры. В этом случае второй запрос, который далее
    #необходим для выполнения алгоритма приравнивается к значению первого запроса, просто потому что второй нужен для дальнейшей работы алгоритма, как я уже
    #и написал.
    #В случае же если SEARCH_SPECIFIC_PRIMER = ['checked'], то из хейдера первого запроса будет получена ссылка на страницу специфического подбора праймеров, которая
    #далее будет использована во втором запросе для получения хтмл страницы со специфиечскими праймерами.
    q = 'Making primers specific to your PCR template.'
    pb_server_status = 'ok'
    if q in r.text:
        pb_server_status = 'overloaded'
        print('---\nServer is overloaded\n---')
    else:
        print('---\nServer is OK\n---')

    return pb_server_status


def random_number():
    import random
    random_number = random.randint(1, 2048)
    return random_number
