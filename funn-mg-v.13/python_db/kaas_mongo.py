# Author: Leandro Correa
# Date: 14.03.2017

import pandas as pd
import sys
import os
import re

from pymongo import MongoClient
from progressbar.progressbar import *

# PATH_kaas = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/kaas/CANGA-MG34_MGRAST_N.kaas'
# PATH_metadata = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/kaas/metadata.csv'

args = sys.argv

if '--help' in args:
    os.system('clear')
    print 'Script: Insert a kaas output into the Mongo metagenomic database.'
    print 'Author: Leandro Correa - @hscleandro'
    print 'Date: 14.03.2017\n'
    print 'How to use: python kaas_mongo.py -i INPUT_KAAS -s SAMPLE -p PROJECT --time\n'
    print 'SAMPLE [required]: Sample name.'
    print 'PROJECT [required]: Project name.'
    print 'INPUT_KAAS [required]: File containing the kaas output in .csv format file Containing two fields in the table header:'
    print '1- Sequence ID'
    print '2- Ko familie'
    print 'time: Graph indicating the total and expected completion time of the execution.\n\n'
    sys.exit('')

else:
    if '-i' in args:
        PATH_kaas = args[args.index('-i') + 1]
        split = str.split(PATH_kaas, '/')
        PATH = ''
        for s in range(1, len(split[:-1])):
            PATH = PATH + '/' + split[s]
        PATH += '/'
    else:
        sys.exit(
            "\n\nErro: Parameter -i required for script execution. \n\nUse: python kaas_mongo.py --help for details.\n"
        )
    if '-s' in args:
        sample = args[args.index('-s') + 1]
        sample = sample.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -s required for script execution. \n\nUse: python kaas_mongo.py --help for details.\n"
        )
    if '-p' in args:
        project = args[args.index('-p') + 1]
        project = project.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -p required for script execution. \n\nUse: python kaas_mongo.py --help for details.\n"
        )

    if '--time' in args:
        print_time = True
    else:
        print_time = False

    address = os.getcwd()
    collums = ["V1", "V2"]
    config_address = address.replace('src', 'setup.conf')
    #print "config_address: " + config_address + "\n"
    config = pd.read_csv(config_address, sep="=", names=collums)

    port = config.iloc[1]["V2"]
    port = port.replace('"', '')
    port = port.replace(' ', '')
    port = int(port)
    host = config.iloc[0]["V2"]
    host = host.replace('"', '')
    host = host.replace(' ', '')

    dbname = "kaas"
    kaas_collection = project + "-" + sample
    client = MongoClient(host, port)
    db = client[dbname]
    collection = db[kaas_collection]

    dbname = "manager"
    dbs = client[dbname]
    collection_sample = dbs.samples

    columns = ['read_id',
               'ko'
               ]
    kaas_df = pd.read_csv(PATH_kaas, sep="\t", names=columns)
    if print_time:
        print '\nLoading. . .\n'
        widgets = ['Update: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(), ' ',
                   FileTransferSpeed()]

        print str(len(kaas_df.index)) + ' instances to be inserted in the mongo database.\n\n'
        pbar = ProgressBar(widgets=widgets, maxval=len(kaas_df.index) * 1000).start()

    kaas_tool = "loading"
    update = collection_sample.update({"$and":
                                            [{
                                                "sample_name": sample,
                                                "project": project,
                                            }]},
                                        {"$set": {"kaas_tool": kaas_tool}
                                         }, upsert=False)
    if not update.get('updatedExisting'):
        item = {'sample_name': sample.upper(),
                'project': project.upper(),
                'kaas_tool': kaas_tool
                }
        collection_sample.insert(item)

    # i = 3
    batch_number = 50000
    batch_documents = [];
    collection.drop()
    for i in range(0, len(kaas_df.index)):
        read_id = kaas_df.iloc[i]['read_id']
        read_id = read_id.replace(" ", "")
        ko = kaas_df.iloc[i]['ko']

        document = {'id_seq': read_id,
                    'kegg_ko': re.sub(" ", "", str(ko))
                   }
        batch_documents.append(document)
        if (i + 1) % batch_number == 0:
            collection.insert(batch_documents)
            batch_documents = [];
        if print_time:
            pbar.update(1000 * i + 1)
    if print_time:
        pbar.finish()

    kaas_tool = "OK"
    update = collection_sample.update({"$and":
        [{
            "sample_name": sample,
            "project": project,
        }]},
        {"$set": {"kaas_tool": kaas_tool}
         }, upsert=False)

print "\n\nThe kaas data was successfully stored."
