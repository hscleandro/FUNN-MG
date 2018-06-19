from pymongo import MongoClient
import sys
import os
import pandas as pd

args = sys.argv

if '--help' in args:
    os.system('clear')
    print 'Script: Delet a sample of the Mongo metagenomic database.'
    print 'Author: Leandro Correa - @hscleandro'
    print 'Date: 19.06.2017\n'
    print 'How to use: python remove_items.py -s SAMPLE -p PROJECT -t TOOL\n'
    print 'SAMPLE [required]: Sample name.'
    print 'PROJECT [required]: Project name.'
    print 'Tool [required]: Base tool of the results of entry into the database:\n' \
          '* Blast\n' \
          '* Interpro\n' \
          '* Kaas\n' \
          '* Kaiju\n' \
          '* Metadata\n' \
          '* Proteomic\n' \
          '* All\n'
    sys.exit('')
else:
    if '-s' in args:
        sample = args[args.index('-s') + 1]
    else:
        sys.exit(
            "\n\nErro: Parameter -s required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
        )
    if '-p' in args:
        project = args[args.index('-p') + 1]
    else:
        sys.exit(
            "\n\nErro: Parameter -p required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
        )
    if '-t' in args:
        tool = args[args.index('-t') + 1]
        tool = tool.upper()
    else:
        sys.exit(
            "\n\nErro: Parameter -t required for script execution. \n\nUse: python proteomic-hit_mongo.py --help for details.\n"
        )

    address = os.getcwd()
    collums = ["V1", "V2"]
    config_address = address.replace('python_db', '') + '/setup.conf'

    config = pd.read_csv(config_address, sep="=", names=collums)

    port = config.iloc[1]["V2"]
    port = port.replace('"', '')
    port = port.replace(' ', '')
    port = int(port)
    host = config.iloc[0]["V2"]
    host = host.replace('"', '')
    host = host.replace(' ', '')

    client = MongoClient(host, port)
    db = client.local
    collection = db.sequences
    collection_sample = db.samples

    if tool == "KAIJU":
        collection_sample.update({"project": project, "sample_name": sample
                                 },
                                 {"$set": {"kaiju_tool": "excluding"}
                                 }
                                 )

        collection.update({"project": project,
                           "id_sample": sample
                          },
                          {"$unset": {"kingdom": 1,
                                      "phylum": 1,
                                      "class": 1,
                                      "order": 1,
                                      "family": 1,
                                      "genre": 1,
                                      "species": 1,
                                      "id_taxon": 1
                                      }
                          }, multi=True)

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                 },
                                 {"$unset": {"kaiju_tool": 1}}
                                 )
    elif tool == "KAAS":
        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$set": {"kaas_tool": "excluding"}
                                  }
                                 )

        collection.update({"project": project,
                           "id_sample": sample
                           },
                          {"$unset": {"kegg_ko": 1
                                      }
                           }, multi=True)

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$unset": {"kaas_tool": 1}}
                                 )
    elif tool == "BLAST":
        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$set": {"blast_tool": "excluding"}
                                  }
                                 )

        collection.update({"project": project,
                           "id_sample": sample
                           },
                          {"$unset": {"blast_id": 1,
                                      "blast_hit": 1,
                                      "blast_score": 1,
                                      "blast_evalue": 1
                                      }
                           }, multi=True)

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$unset": {"blast_tool": 1}}
                                 )
    elif tool == "INTERPRO":
        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$set": {"interpro_tool": "excluding"}
                                  }
                                 )
        collection.update({"project": project,
                           "id_sample": sample
                           },
                          {"$unset": {"orfs_inf": 1
                                      }
                           }, multi=True)

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$unset": {"interpro_tool": 1}}
                                 )
    elif tool == "PROTEOMICS":
        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$set": {"proteomic_tool": "excluding"}
                                  }
                                 )
        collection.update({"project": project,
                           "id_sample": sample
                           },
                          {"$unset": {"proteomics": 1
                                      }
                           }, multi=True)

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$unset": {"proteomic_tool": 1}}
                                 )
    elif tool == "METADATA":
        list = ["project", "sample_name", "kaas_tool", "kaiju_tool", "blast_tool", "interpro_tool"
                "proteomic_tool", "_id"]

        for doc in collection_sample.find({"project": project, "sample_name": sample}):
            items = doc

        items_key = items.keys()
        remove = [item for item in items_key if item not in list]
        element = {}
        for item in remove:
            element[str(item)] = 1

            collection_sample.update({"project": project,
                                      "sample_name": sample
                                     },
                                     {"$unset": element}
                                    )
    elif tool == "ALL":
        list = ["kaas_tool", "kaiju_tool", "kaas_tool", "blast_tool", "interpro_tool", "proteomic_tool"]

        collection_sample.update({"project": project,
                                  "sample_name": sample
                                  },
                                 {"$set": {"sample_status": "excluding"}}
                                 )

        for doc in collection_sample.find({"project": project, "sample_name": sample}):
            items = doc

        items_key = items.keys()
        remove_tools = [item for item in items_key if item in list]

        element = {}
        for item in remove_tools:
            element[str(item)] = "excluding"

            collection_sample.update({"project": project,
                                      "sample_name": sample
                                      },
                                     {"$set": element}
                                     )

        collection.remove({"project": project,
                           "id_sample": sample
                           }
                          )

        collection_sample.remove({"project": project,
                                  "sample_name": sample
                                 }
                                )



print "\nThe sample was deleted successfully."
