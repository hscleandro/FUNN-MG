# Author: Leandro Correa
# Date: 14.03.2017

import re
import pandas as pd
import sys
import os
from progressbar.progressbar import *
import subprocess

from pymongo import MongoClient

# proteomic_file = '/home/leandro/Data/metagenomas/MG_34_Emma/contigs_newbler/hit_protein/Canga_hit_proteins_mgrast_FINAL.fasta'
# PATH_metadata = '/home/leandro/Data/metagenomas/MG_34_Emma/hit_protein/metadata.csv'

args = sys.argv

if '--help' in args:
    os.system('clear')
    print 'Script: Insert a proteomic table output into the Mongo metagenomic database.'
    print 'Author: Leandro Correa - @hscleandro'
    print 'Date: 14.03.2017\n'
    print 'How to use: python proteomic_mongo.py -i INPUT_TABLE -s SAMPLE -p PROJECT --time\n'
    print 'SAMPLE [required]: Sample name.'
    print 'PROJECT [required]: Project name.'
    print 'INPUT_TABLE [required]: File containing the result of proteomic analyse. ' \
          'The proteomic output must be the fasta format contain the header and the aminoacids' \
          ' of each sequence.'
    print 'time: Graph indicating the total and expected completion time of the execution.\n\n'
    sys.exit('')
else:
    if '-i' in args:
        proteomic_file = args[args.index('-i') + 1]

        split = str.split(proteomic_file, '/')
        PATH = ''
        for s in range(1, len(split[:-1])):
            PATH = PATH + '/' + split[s]
        PATH += '/'
    else:
        sys.exit(
            "\n\nErro: Parameter -i required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
        )
    if '-s' in args:
        sample = args[args.index('-s') + 1]
        sample = sample.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -s required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
        )
    if '-p' in args:
        project = args[args.index('-p') + 1]
        project = project.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -p required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
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

    client = MongoClient(host, port)

    dbname = "proteomics"
    proteomics_collection = project + "-" + sample
    db = client[dbname]
    collection = db[proteomics_collection]

    dbname = "manager"
    dbs = client[dbname]
    collection_sample = dbs.samples

    if print_time:
        grep = "grep -c '^>' "
        command = grep + proteomic_file
        count = subprocess.check_output(command, shell=True)
        count = count[:-1]
        print '\nLoading. . .\n'
        widgets = ['Update: ', Percentage(), ' ', Bar(marker=RotatingMarker()), ' ', ETA(), ' ',
                   FileTransferSpeed()]

        pbar = ProgressBar(widgets=widgets, maxval=int(count) * 10).start()


    proteomic_tool = "loading"
    update = collection_sample.update({"$and":
        [{
            "sample_name": sample,
            "project": project,
        }]},
        {"$set": {"proteomic_tool": proteomic_tool}
         }, upsert=False)
    if not update.get('updatedExisting'):
        item = {'sample_name': sample,
                'project': project,
                'proteomic_tool': proteomic_tool
                }
        collection_sample.insert(item)

    i = 1
    for line in open(proteomic_file, 'r'):
        if i % 2 == 1:
            #line = ">contig00003_10126_11573_-"
            read_id = line[1:]
            read_id = re.sub("\n", "", read_id)
            #sequence = str.split(read_id, "_")[0]
            #"""
            update = collection.update({'id_seq': read_id},
                                       {'$set': {'proteomics': "true"
                                                  },
                                        }, upsert=False)
            #print read_id + "\t" + str(update.get('updatedExisting'))
            if not update.get('updatedExisting'):
                item = {'id_seq': read_id,
                        'proteomics': "true"
                        }
                ObjectId = collection.insert(item)

            if print_time:
                pbar.update(i + 1)
        if print_time:
            pbar.finish()
        i += 1
    proteomic_tool = "OK"
    update = collection_sample.update({"$and":
        [{
            "sample_name": sample,
            "project": project,
        }]},
        {"$set": {"proteomic_tool": proteomic_tool}
         }, upsert=False)

print "\n\nThe proteomic data was successfully stored."
