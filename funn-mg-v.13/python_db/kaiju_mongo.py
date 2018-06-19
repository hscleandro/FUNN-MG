# Author: Leandro Correa
# Date: 08.03.2017

import pandas as pd
import sys
import os
import re
from progressbar.progressbar import *

from pymongo import MongoClient

# PATH_kaiju = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/kaiju/mgm4723928.fnn.kaiju.out.names'
# PATH_metadata = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/kaiju/metadata.csv'

args = sys.argv


if '--help' in args:
    os.system('clear')
    print 'Script: Insert a kaiju output into the Mongo metagenomic database.'
    print 'Author: Leandro Correa - @hscleandro'
    print 'Date: 09.03.2017\n'
    print 'How to use: python kaiju_mongo.py -i INPUT_KAIJU -s SAMPLE -p PROJECT --time\n'
    print 'SAMPLE [required]: Sample name.'
    print 'PROJECT [required]: Project name.'
    print 'INPUT_KAIJU [required]: File containing the kaiju output.\n'
    print 'Input file format: The kaiju output must contain 4 fields:'
    print '1- Tag to indicate whether the sequence was annotated or not (U or C)'
    print '2- Sequence ID'
    print '3- Taxon hierarchy'
    print '4- Taxon name'
    print 'time: Graph indicating the total and expected completion time of the execution.\n\n'
    sys.exit('')

else:
    if '-i' in args:
        PATH_kaiju = args[args.index('-i') + 1]

        split = str.split(PATH_kaiju, '/')
        PATH = ''
        for s in range(1, len(split[:-1])):
            PATH = PATH + '/' + split[s]
        PATH += '/'
    else:
        sys.exit(
            "\n\nErro: Parameter -i required for script execution. \n\nUse: python kaiju_mongo.py --help for details.\n"
        )
    if '-s' in args:
        sample = args[args.index('-s') + 1]
        sample =sample.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -s required for script execution. \n\nUse: python kaiju_mongo.py --help for details.\n"
        )
    if '-p' in args:
        project = args[args.index('-p') + 1]
        project = project.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -p required for script execution. \n\nUse: python kaiju_mongo.py --help for details.\n"
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
    port = port.replace('"','')
    port = port.replace(' ', '')
    port = int(port)
    host = config.iloc[0]["V2"]
    host = host.replace('"', '')
    host = host.replace(' ', '')

    client = MongoClient(host, port)

    dbname = "kaiju"
    taxon_collection = project + "-" + sample
    db = client[dbname]
    collection = db[taxon_collection]

    dbname = "manager"
    dbs = client[dbname]
    collection_sample = dbs.samples

    columns = ['status',
               'read_id',
               'taxon_id',
               'taxon_name']
    kaiju_df = pd.read_csv(PATH_kaiju, sep="\t", names=columns)
    if print_time:
        print '\nLoading. . .\n'
        widgets = ['Update: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(), ' ',
                   FileTransferSpeed()]

        print str(len(kaiju_df.index)) + ' instances to be inserted in the mongo database.\n\n'
        pbar = ProgressBar(widgets=widgets, maxval=len(kaiju_df.index) * 1000).start()

    kaiju_tool = "loading"
    update = collection_sample.update({"$and":
                                            [{
                                                "sample_name": sample,
                                                "project": project,
                                            }]},
                                        {"$set": {"kaiju_tool": kaiju_tool}
                                         }, upsert=False)
    if not update.get('updatedExisting'):

        item = {'sample_name': sample,
                'project': project,
                'kaiju_tool': kaiju_tool
                }
        collection_sample.insert(item)
    # i = 119188
    for i in range(0, len(kaiju_df.index)):
        read_id = kaiju_df.iloc[i]['read_id']
        read_id = read_id.replace(" ", "")
        id_taxon = kaiju_df.iloc[i]['taxon_id']
        name_taxon = kaiju_df.iloc[i]['taxon_name']
        if str(type(name_taxon)) == "<type 'str'>":
            taxons_hierarchy = str.split(name_taxon, ';')[0:7]
        else:
            taxons_hierarchy = ['NA', 'NA', 'NA', 'NA', 'NA', 'NA', 'NA']

        update = collection.update({'id_seq': read_id},
                                           {'$set': {'id_taxon': str(id_taxon).strip(),
                                                     'kingdom': str(taxons_hierarchy[0]).strip(),
                                                     'phylum': str(taxons_hierarchy[1]).strip(),
                                                     'class': str(taxons_hierarchy[2]).strip(),
                                                     'order': str(taxons_hierarchy[3]).strip(),
                                                     'family': str(taxons_hierarchy[4]).strip(),
                                                     'genre': str(taxons_hierarchy[5]).strip(),
                                                     'species': str(taxons_hierarchy[6]).strip()
                                                     },
                                            }, upsert=False)
        if not update.get('updatedExisting'):
            item = {'id_seq': read_id,
                    'id_taxon': str(id_taxon).strip(),
                    'kingdom': str(taxons_hierarchy[0]).strip(),
                    'phylum': str(taxons_hierarchy[1]).strip(),
                    'class': str(taxons_hierarchy[2]).strip(),
                    'order': str(taxons_hierarchy[3]).strip(),
                    'family': str(taxons_hierarchy[4]).strip(),
                    'genre': str(taxons_hierarchy[5]).strip(),
                    'species': str(taxons_hierarchy[6]).strip()
                    }
            ObjectId = collection.insert(item)
        if print_time:
            pbar.update(1000 * i + 1)
    if print_time:
        pbar.finish()

    kaiju_tool = "OK"
    update = collection_sample.update({"$and":
        [{
            "sample_name": sample,
            "project": project,
        }]},
        {"$set": {"kaiju_tool": kaiju_tool}
         }, upsert=False)

print "\n\nThe taxon data was successfully stored."
