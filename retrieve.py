#! /usr/bin/env python

import uniprot
import sys

def retrieve_uniprot_meta(uniprotIDs_file):
    with open(uniprotIDs_file, 'r') as f:
        id_list = [uni_id.strip() for uni_id in f]
        
    meta_dict = uniprot.batch_uniprot_metadata(id_list)
    return meta_dict

def print_uniprot_meta(meta_dict):

    for uni_id in meta_dict:
        data = meta_dict[uni_id]
        print "%s\t%s\t%s\t%s" % (uni_id,
                                  data.get('id', None),
                                  data.get('description', None),
                                  data.get('go', None),)

if __name__ == '__main__':

    meta_dict = retrieve_uniprot_meta(sys.argv[1])
    print_uniprot_meta(meta_dict)
    
    

