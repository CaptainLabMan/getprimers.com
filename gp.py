    #!/Library/Frameworks/Python.framework/Versions/3.8/bin/python3


    # От названия гена получить ENSG, от него получить список ENST, выбрать нужный транскрипт и выудить из него ENSE всех экзонов, далее получить их позиции, отнять, прибавить получить последовательность, отправить ее в праймербласт.
    # Часть которая добывает id по названию

def get_transcripts(gene, GCS):
    import requests, sys, time, threading

    # От названия гена получить ENSG, от него получить список ENST, выбрать нужный транскрипт и выудить из него ENSE всех экзонов, далее получить их позиции, отнять, прибавить получить последовательность, отправить ее в праймербласт.
    # Часть которая добывает id по названию

    settings = ''

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"

    genename = gene
    print(genename)

    ext = '/xrefs/symbol/homo_sapiens/{}?'.format(genename)
    #ext = "/xrefs/symbol/homo_sapiens/HGD?"

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    id = decoded[0]
    id = id['id']
    print('\nGene ID:')
    print(id + "\n")




    # Эта чатсь может доставать инфу по транскриптам в зависимости от features
    ext = "/overlap/id/{}?feature=transcript".format(id)
    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    #for element in decoded:
    #    print(element)


    transcripts = []
    for element in decoded:
        transcripts.append(element['transcript_id'])
    #    transcripts.append('{}.{}'.format(element['transcript_id'], element['version']))
    #print(transcripts)

    #Формирую пустой список для добавления туда транскриптов у которых есть РефСек, для дальнейшего объединения списка транскриптов без РефСек и с РефСек.
    transcripts_and_refseq_list = []

    print('Transcripts ID:')

    # Часть по вытаскиванию RefSeq и Часть по выводу на экран инфы по транскриптан с рефсек
    for transcript_id_element in transcripts:
        ext = "/xrefs/id/{}?".format(transcript_id_element)
        r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})
        if not r.ok:
          r.raise_for_status()
          sys.exit()

        decoded = r.json()
    #    for element in decoded:
    #        print(element)
        for dict_element in decoded:
            if dict_element['dbname'] == 'RefSeq_mRNA':
                primary_id = dict_element['primary_id']
                transcripts_and_refseq_list.append('{} - {}'.format(transcript_id_element, primary_id))
    #            print('{} - {}'.format(transcript_id_element, primary_id))
                if transcript_id_element in transcripts:
                    transcripts.remove(transcript_id_element)  # при наличии у транскрипта refseq он будет удалять из списка transcripts чтобы этот транскрипт не дублировался

    #Здесь объединяю списки транскриптов без РефСек и с РефСек.
    #print(transcripts_and_refseq_list)
    for singe_transcript_and_refseq in transcripts_and_refseq_list:
        transcripts.append(singe_transcript_and_refseq)
    #Переворачиваю лист транскриптов чтобы транскрипты с NM были вверху, а не внизу.
    transcripts.reverse()
    # Часть по выводу на экран инфы по транскриптам у которых не refseq

    #print('Transcripts ID:')
#     for element in transcripts:
#        print(element)
#    print('')

    #print('Transcripts ID:')
    with open('params.py', 'a') as t:
        t.write('transcripts = {}\n'.format(transcripts))

    return transcripts


def get_exons(transcript, GCS):
    import requests, sys, time, threading, json


    transcript_id = transcript

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"

    # Эта чатсь может доставать инфу по экзонам в зависимости от features
    ext = "/overlap/id/{}?feature=exon".format(transcript_id)

    r = requests.get(server+ext, headers={ "Content-Type" : "application/json"})

    if not r.ok:
      r.raise_for_status()
      sys.exit()

    decoded = r.json()
    #for element in decoded:
    #    print(element)
    #    print('')

    strand = ''
    chromosome = ''

    exons_id = []
    for element in decoded:
        if element['Parent'] == transcript_id:
            exons_id.append(element['exon_id'])
            strand = element['strand']               # Определение направления цепи
            chromosome = element['seq_region_name']  # Определение [хромосомы]
    #        exons_id.append(element['start'])
    #        exons_id.append(element['end'])
#    print("\nКоличество экзонов: {}\n".format(len(exons_id)))

    #Рассчитывает количество экзонов в транскрипте для передаче html шаблонизатору в exons.html
    exons_count = []
    for nubmer in range(1, len(exons_id) + 1):
        exons_count.append(nubmer)

    #print(exons_id)
    #print(strand)
    #print(chromosome)

    dict_exons = dict()
    for element in decoded:
        if element['Parent'] == transcript_id:
            dict_exons[element['exon_id']] = {'start': element['start'], 'end': element['end']}

    #print(dict_exons)
    #for element in decoded:
    #    print(element)
    #    print('\n')
    #print(decoded)

    #Формирование дикта для сериализации данных в json
    params_json = {'chromosome': chromosome, 'strand': str(strand), 'exons_id': exons_id, 'dict_exons': dict_exons}
    #print(params_json)
    #Запись номера хромосомы и цепи в файл params.json
    with open('params.json', 'w') as params:
        json.dump(params_json, params)

    #Запись номера хромосомы и цепи в файл params.py
    with open('params.py', 'a') as params:
        params.write('exons_id = {}\n'.format(exons_id))
        params.write('chromosome = {}\n'.format(str(chromosome)))
        params.write('strand = {}\n'.format(str(strand)))
        params.write('dict_exons = {}\n'.format(dict_exons))

    #print(exons_count)
    return exons_count



def get_primers(chromosome, strand, taken_exons, exons_id, dict_exons, SEARCH_SPECIFIC_PRIMER, NO_SNP, GCS):
    import requests, sys, time, threading, re, json
    from bs4 import BeautifulSoup

    #print('\nTHIS')
    #print(dict_exons)
    #print('THIS\n')

    # Формирование списка id выбранных экзонов
    taken_exons_id = []
    for element in taken_exons:
        taken_exons_id.append(exons_id[int(element) - 1])
    #print(taken_exons_id)

    if GCS == 'GRCh37':
        server = "https://grch37.rest.ensembl.org"
    elif GCS == 'GRCh38':
        server = "https://rest.ensembl.org"
    # Часть по определению графниц ампликона
#    if settings == 'нет' or settings == 'no':
#        extern_range = input('Введите наружный диапазон отступа от экзона (540):')
#        if extern_range == '':
#            extern_range = 540
#        if extern_range != '':
#            extern_range = int(extern_range)
#
#        intern_range = input('Введите внутренний диапазон отступа от экзона (65):')
#        if intern_range == '':
#            intern_range = 65
#        if intern_range != '':
#            intern_range = int(intern_range)
#    else:
#        extern_range = 540
#        intern_range = 65


    # Часть для получения последовательностей для ПраймерБласта + Часть по рассчетам для вставки в Праймер бласт

    print('\nPrimer design is going.\nPlease wait...\n')

    global result_dict
    result_dict = {}


    def get_get_primers(element):

        #Тут по сути собирается запрос для получения последовательности для поиска праймеров (последовательонсти вместе с отступом для поиска).
        extern_range = 540
        intern_range = 65
        ext = "/sequence/region/human/{}:{}..{}:{}?expand_3prime=0;expand_5prime=0".format(chromosome, int(dict_exons[element]['start']) - extern_range, int(dict_exons[element]['end']) + extern_range, strand)
        r_exon_and_distance = requests.get(server+ext, headers={ "Content-Type" : "text/plain"})
        if not r_exon_and_distance.ok:
          r_exon_and_distance.raise_for_status()
          sys.exit()

    #      return start_of_exon = int(dict_exons[element]['start'])
    #      return end_of_exon = int(dict_exons[element]['end'])

        #Здесь часть похожая на часть сверху, по получению последовательности, но это последовательность строго самого кэзона (+UTR если он есть).
        #Нужна эта часть для того чтобы выделить большими буквами экзон, а маленькими интроны в последовательности.
        ext_exon_only = "/sequence/region/human/{}:{}..{}:{}?expand_3prime=0;expand_5prime=0".format(chromosome, int(dict_exons[element]['start']), int(dict_exons[element]['end']), strand)
        r_exon_only = requests.get(server+ext_exon_only, headers={ "Content-Type" : "text/plain"})
        if not r_exon_only.ok:
          r_exon_only.raise_for_status()
          sys.exit()
        r_exon_only = r_exon_only.text
        #print(r_exon_only.text)

        #Строки ниже необходимы просто для вывода инфы о дистанции для поиска праймеров в консоль.
        distance = '1 - {}\n{} - {}'.format(int(extern_range) - int(intern_range), len(r_exon_and_distance.text) - (int(extern_range) - int(intern_range)), len(r_exon_and_distance.text))
        print_dis = distance
        number_of_exon = exons_id.index(element) + 1    # Номер экзона
        print_exon_number = '\nЭкзон {}'.format(number_of_exon)     #Перенесено в конец с часть по принту
        #if extern_range != 0 and intern_range != 0:         #Перенесено в конец с часть по принту
        #    print(print_dis)

        #print(r.text)
        #Переменная ниже - это последователность на которой будут искаться праймеры вместе с отступами для поиска праймеров.
        #Ниже две новые переменные (print_r_text и seq_for_pb), которые несут в себе одно и то же, просто нужны для разных целей. В принципе можно и одну сдлеать.
        print_r_text = (r_exon_and_distance.text).lower()
        print_r_text = print_r_text.replace(r_exon_only.lower(), r_exon_only.upper())
        exon_sequence = {'sequence': print_r_text}  #Создание словаря, в котром будет храниться последовательность ЭКЗОНА.
        #print(exon_sequence)
        seq_for_pb = '\'{}\''.format(r_exon_and_distance.text)

        #print('')    #Перенесено в конец с часть по принту

        #Часть по рассчету позиций поиска праймеров для вставки в форму PB.
        PRIMER5_END_PB = int(extern_range) - int(intern_range)
        PRIMER3_START_PB = len(r_exon_and_distance.text) - (int(extern_range) - int(intern_range))
        PRIMER3_END_PB = len(r_exon_and_distance.text)



        #Часть по получению тупо страницы HTML с праймерами от PrimerBlast'a
        #Сначала я загоняю последовательность со всеми параметрами в ПБ и он возвращает мне ссылку
        #Датее эту ссылку я вставляю во второй запрос (r2) и он уже дает мне HTML страницу
        #которую я буду парсить на наличие праймеров (естественно со слипом)
        url = 'https://www.ncbi.nlm.nih.gov/tools/primer-blast/primertool.cgi'
        data = {
            'INPUT_SEQUENCE': seq_for_pb,
            'PRIMER5_START': '1',
            'PRIMER5_END': PRIMER5_END_PB,
            'PRIMER3_START': PRIMER3_START_PB,
            'PRIMER3_END': PRIMER3_END_PB,
            'PRIMER_SPECIFICITY_DATABASE': 'PRIMERDB/genome_selected_species',
            'SEARCHMODE': '2',
            'ORGANISM': 'Homo sapiens',
            'TOTAL_PRIMER_SPECIFICITY_MISMATCH': '1',
            'PRIMER_3END_SPECIFICITY_MISMATCH': '1',
            'MISMATCH_REGION_LENGTH': '5',
            'TOTAL_MISMATCH_IGNORE': '6',
            'PRIMER_MISPRIMING_LIBRARY': 'AUTO',
            'HITSIZE': '50000',
            'EVALUE': '30000',
            'WORD_SIZE': '7',
            'LOW_COMPLEXITY_FILTER': 'checked',
            'PRIMER_PRODUCT_MIN': '70',
            'PRIMER_PRODUCT_MAX': '1500',
            'PRIMER_MIN_SIZE': '15',
            'PRIMER_OPT_SIZE': '20',
            'PRIMER_MAX_SIZE': '25',
            'PRIMER_MIN_TM': '57.0',
            'PRIMER_OPT_TM': '60.0',
            'PRIMER_MAX_TM': '63.0',
            'PRIMER_MAX_DIFF_TM': '3',
            'PRIMER_PRODUCT_MAX': '10000',
            'SEARCH_SPECIFIC_PRIMER': SEARCH_SPECIFIC_PRIMER,
            'NO_SNP': NO_SNP
        }

        print(SEARCH_SPECIFIC_PRIMER)
        print(NO_SNP)

        r = requests.post(url, data=data)
        headers_pb = r.headers

        #Здесь условная конструкция, которая если SEARCH_SPECIFIC_PRIMER = [], говорить о том, что выбран неспецифический подбор праймеров
        #а значит нет смысла брать ссылку из первого запроса на второй запрос, который специфически подбирает праймеры. В этом случае второй запрос, который далее
        #необходим для выполнения алгоритма приравнивается к значению первого запроса, просто потому что второй нужен для дальнейшей работы алгоритма, как я уже
        #и написал.
        #В случае же если SEARCH_SPECIFIC_PRIMER = ['checked'], то из хейдера первого запроса будет получена ссылка на страницу специфического подбора праймеров, которая
        #далее будет использована во втором запросе для получения хтмл страницы со специфиечскими праймерами.
        if SEARCH_SPECIFIC_PRIMER == []:
            r2 = r
        elif SEARCH_SPECIFIC_PRIMER == ['checked']:
            pb_specific_link = '{}'.format(headers_pb['NCBI-RCGI-RetryURL'])
            #Вот тут второй запрос и слип
            #Тут нужно грамотно разместить слип чтобы программа могла делать проверки
            #готовности праймеров и в случае не готовности делать слип и потом снова
            #повторять запрос
            r2 = requests.get(pb_specific_link)
            q = 'Making primers specific to your PCR template.'
            while q in r2.text:
            #        print('Request')
                time.sleep(15)
                r2 = requests.get(pb_specific_link)

        #Вот тут второй запрос и слип
        #Тут нужно разместить слип чтобы программа могла делать проверки
        #готовности праймеров и в случае не готовности делать слип и потом снова
        #повторять запрос
        #Здесь (ниже) тоже самое, что и в алгоритме выше, просто пока сотавил.
#        r2 = requests.get(pb_specific_link)
#        q = 'Making primers specific to your PCR template.'
#        while q in r2.text:
    #        print('Request')
#            time.sleep(10)
#            r2 = requests.get(pb_specific_link)

        print('')

        #Тут нужно грамотно разместить слип чтобы программа могла делать проверки
        #готовности праймеров и в случае не готовности делать слип и потом снова
        #повторять запрос
        #После того как q не будет в р2.текст, программа будет идти дальше
        html_primers = r2.text


        #Часть по парсингу и получению только нужных строк

        soup = BeautifulSoup(html_primers, 'html.parser')
        items = soup.find_all('div', class_='prPairInfo')
        #print(items)

        #Тут создается словарь в который будут включены все элементы всех праймеров
        all_primers_dicts = {}

        #Часть по принту того что намного выше
        if extern_range != 0 and intern_range != 0:
            print_print_dis = print_dis
        seq_string = '{}\n{}\n{}\n{}'.format(print_exon_number, print_print_dis, print_r_text, '')


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
            primers_dict = ''
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


#            result_dict[number_of_primerpair]

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
        result_dict['Exon_{}'.format(number_of_exon)]=all_primers_dicts
        #Тоже самое, только на уровень ниже, добавляю последовательность экзона в дикт.
        result_dict['Exon_{}'.format(number_of_exon)]['exon_sequence']=print_r_text


    #Формирование и запуск разных потоков
    threads = []
    for element in taken_exons_id:
        t = threading.Thread(target=get_get_primers, args=[element])
        t.start()
        threads.append(t)
        time.sleep(1) #Здесь слип необходим потому что между запуском потока нужны интервалы, иначе библиотека Requests жаулется на слишком большок количество запросов за единицу времени, поэтому искусственно замедляем.
        print('Поток запущен.')
    for t in threads:
        t.join()

    exons_result_dict = []
    for ee in result_dict:
        exons_result_dict.append(ee)

    #result_dict = json.dumps(result_dict, sort_keys=True)
    #result_dict = json.loads(result_dict)

    with open('params.py', 'a') as params:
        params.write('result_dict = {}\n'.format(result_dict))
        params.write('exons_result_dict = {}\n'.format(exons_result_dict))

    return result_dict
