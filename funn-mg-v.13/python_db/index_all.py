from pymongo import MongoClient
import pymongo as mongo
import pandas as pd
import sys
import os

args = sys.argv

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

if '-s' in args:
    sample = args[args.index('-s') + 1]
    sample = sample.upper()
else:
    sys.exit(
        "\n\nErro: Parameter -s required for script execution. \n"
    )
if '-p' in args:
    project = args[args.index('-p') + 1]
    project = project.upper()
else:
    sys.exit(
        "\n\nErro: Parameter -p required for script execution. \n"
    )

taxon_collection = project + "-" + sample
kaas_collection = project + "-" + sample
funn_collection = project + "-" + sample
proteomics_collection = project + "-" + sample
sample = "manager"
taxon = "kaiju"
functional = "kaas"
metabolic = "funn"
proteomics = "proteomics"


client = MongoClient(host, port)

db = client[sample]
all_collection = db.collection_names()

if 'samples' not in all_collection:
    collection = db['samples']
    atributes = ["project", "sample_name"]

    for id in atributes:
        collection.ensure_index([(id, mongo.ASCENDING)])

db = client[taxon]
all_collection = db.collection_names()

if taxon_collection not in all_collection:
    collection = db[taxon_collection]
    atributes = ["id_seq", "kingdom", "phylum", "species",
                 "family", "genre", "class", "id_taxon", "order"]

    for id in atributes:
        collection.ensure_index([(id, mongo.ASCENDING)])

db = client[functional]
all_collection = db.collection_names()

if kaas_collection not in all_collection:
    collection = db[kaas_collection]
    atributes = ["id_seq", "kegg_ko"]

    for id in atributes:
        collection.ensure_index([(id, mongo.ASCENDING)])

db = client[metabolic]
all_collection = db.collection_names()

if funn_collection not in all_collection:
    collection = db[funn_collection]
    atributes = ["betweenness_centrality",
                 "degree","load_roume","load_rahman","choke_point",
                 "ko_functions","kegg_ko","paths","pathways",
                 "class","subclass","Function","genes_relation_sample",
                 "genes_relation_noted","coverage","p_value","q_value",
                 "betweenness_centrality","kos"]

    for id in atributes:
        collection.ensure_index([(id, mongo.ASCENDING)])

db = client[proteomics]
all_collection = db.collection_names()

if proteomics_collection not in all_collection:
    collection = db[proteomics_collection]
    atributes = ["id_seq", "proteomics"]

    for id in atributes:
        collection.ensure_index([(id, mongo.ASCENDING)])

print "\n\nIndexing completed successfully."