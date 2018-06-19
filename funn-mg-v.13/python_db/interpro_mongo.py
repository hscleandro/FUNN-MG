# Author: Leandro Correa
# Date: 08.03.2017

import pandas as pd
import sys
import os
from progressbar.progressbar import *

from pymongo import MongoClient

args = sys.argv

# PATH_emg = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/interpro/contig/Canga_emma_hit_proteins_newbler_mgrast_FINAL.tsv'

if '--help' in args:
    os.system('clear')
    print 'Script: Insert a .tsv output from EMG pipeline into the Mongo metagenomic database.'
    print 'Author: Leandro Correa - @hscleandro'
    print 'Date: 08.03.2017\n'
    print 'How to use: python interpro_mongo.py -i INPUT_EMG.TSV -s SAMPLE -p PROJECT -time\n'
    print 'SAMPLE [required]: Sample name.'
    print 'PROJECT [required]: Project name.'
    print 'INPUT_EMG.TSV [required]: File containing the .tsv table of EMG pipeline. The .tsv file must contain 15 fields.\n' \
          'For more details go to: https://github.com/ebi-pf-team/interproscan/wiki/InterProScan5OutputFormats.\n\n'
    print 'time: Graph indicating the total and expected completion time of the execution.\n\n'
    sys.exit('')
else:
    if '-i' in args:
        PATH_emg = args[args.index('-i') + 1]
        split = str.split(PATH_emg, '/')
        PATH = ''
        for s in range(1, len(split[:-1])):
            PATH = PATH + '/' + split[s]
        PATH += '/'
    else:
        sys.exit(
            "\n\nErro: Parameter -i required for script execution. \n\nUse: python interpro_mongo.py --help for details.\n"
        )
    if '-s' in args:
        sample = args[args.index('-s') + 1]
    else:
        sys.exit(
            "\n\nErro: Parameter -s required for script execution. \n\nUse: python interpro_mongo.py --help for details.\n"
        )
    if '-p' in args:
        project = args[args.index('-p') + 1]
    else:
        sys.exit(
            "\n\nErro: Parameter -p required for script execution. \n\nUse: python interpro_mongo.py --help for details.\n"
        )
    if '-time' in args:
        print_time = True
    else:
        print_time = False

    address = os.getcwd()
    collums = ["V1", "V2"]
    config_address = address.replace('src', 'setup.conf')

    config = pd.read_csv(config_address, sep="=", names=collums)

    port = config.iloc[1]["V2"]
    port = port.replace('"', '')
    port = port.replace(' ', '')
    port = int(port)
    host = config.iloc[0]["V2"]
    host = host.replace('"', '')
    host = host.replace(' ', '')

    client = MongoClient(host, port)

    dbname = "interpro"
    interpro_collection = project + "-" + sample
    db = client[dbname]
    collection_interpro = db[interpro_collection]

    dbname = "go"
    go_collection = project + "-" + sample
    db = client[dbname]
    collection_go = db[go_collection]

    dbname = "pathways"
    pathways_collection = project + "-" + sample
    db = client[dbname]
    collection_pathways = db[pathways_collection]

    dbname = "manager"
    dbs = client[dbname]
    collection_sample = dbs.samples

    columns = ['Protein Accession',  # 1
               'Sequence MD5 digest',  # 2
               'Sequence Length',  # 3
               'Analysis',  # 4
               'Signature Accession',  # 5
               'Signature Description',  # 6
               'Start location',  # 7
               'Stop location',  # 8
               'Score',  # 9
               'Status',  # 10
               'Date',  # 11
               'InterPro accession',  # 12
               'InterPro description',  # 13
               'GO annotations',  # 14
               'Pathways annotations']  # 15

    emg_df = pd.read_csv(PATH_emg, sep="\t", names=columns)

    if print_time:
        print '\nLoading. . .\n'
        widgets = ['Update: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(), ' ', FileTransferSpeed()]

        print str(len(emg_df.index)) + ' instances to be inserted in the mongo database.\n\n'
        pbar = ProgressBar(widgets=widgets, maxval=len(emg_df.index) * 1000).start()

    interpro_tool = "loading"
    update = collection_sample.update({"$and":
                                            [{
                                                "sample_name": sample,
                                                "project": project,
                                            }]},
                                        {"$set": {"interpro_tool": interpro_tool}
                                         }, upsert=False)
    if not update.get('updatedExisting'):
        item = {'sample_name': sample.upper(),
                'project': project.upper(),
                'interpro_tool': interpro_tool
                }
        collection_sample.insert(item)

    batch_number = 50000
    batch_documents = []
    batch_go = []
    batch_pathwyas = []
    collection_interpro.drop()
    # i = 0
    for i in range(0, len(emg_df.index)):
        read_id = emg_df.iloc[i]['Protein Accession']
        read_id = read_id.replace(" ", "")
        start_location = emg_df.iloc[i]['Start location']
        stop_location = emg_df.iloc[i]['Stop location']
        signature_accession = emg_df.iloc[i]['Signature Accession']
        signature_description = emg_df.iloc[i]['Signature Description']
        protein_analysis = emg_df.iloc[i]['Analysis']
        Score = emg_df.iloc[i]['Score']
        interpro_accession = emg_df.iloc[i]['InterPro accession']
        interpro_description = emg_df.iloc[i]['InterPro description']

        if type(emg_df.iloc[i]['GO annotations']) is str:
            go_annotations = str.split(emg_df.iloc[i]['GO annotations'], "|")
        else:
            go_annotations = []

        Kegg_Pathways = []
        Metacyc_pathways = []
        Reactome_pathways = []

        if type(emg_df.iloc[i]['Pathways annotations']) is str:

            pathways_annotations = str.split(emg_df.iloc[i]['Pathways annotations'], "|")

            for j in range(0, len(pathways_annotations)):
                if pathways_annotations[j].find('KEGG') == 0:
                    Kegg_Pathways.append(pathways_annotations[j])
                else:
                    if pathways_annotations[j].find('MetaCyc') == 0:
                        Metacyc_pathways.append(pathways_annotations[j])
                    else:
                        if pathways_annotations[j].find('Reactome') == 0:
                            Reactome_pathways.append(pathways_annotations[j])

        document = {'id_seq': read_id,
                    "start_location": str(start_location),
                    "stop_location": str(stop_location),
                    "signature_accession": str(signature_accession),
                    "signature_description": str(signature_description),
                    "protein_analysis": str(protein_analysis),
                    "score": str(Score),
                    "interpro_accession": str(interpro_accession),
                    "interpro_description": str(interpro_description),
                   }
        document_go = {'id_seq': read_id,
                       "go_annotations": str(go_annotations)
                      }
        document_pathways = {'id_seq': read_id,
                             "kegg_Pathways": Kegg_Pathways,
                             "reactome_pathways": Reactome_pathways,
                             "metacyc_pathways": Metacyc_pathways
                             }

        batch_pathwyas.append(document_pathways)
        batch_go.append(document_go)
        batch_documents.append(document)
        if (i + 1) % batch_number == 0:
            collection_interpro.insert(batch_documents)
            collection_go.insert(batch_go)
            collection_pathways.insert(batch_pathwyas)
            batch_documents = []
        if print_time:
            pbar.update(1000 * i + 1)
    if len(batch_documents) > 0:
        collection_interpro.insert(batch_documents)

    if print_time:
        pbar.finish()

    interpro_tool = "OK"
    update = collection_sample.update({"$and":
        [{
            "sample_name": sample,
            "project": project,
        }]},
        {"$set": {"interpro_tool": interpro_tool}
         }, upsert=False)

print "\n\nThe data was successfully stored."
