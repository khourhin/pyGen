#! /usr/bin/env python

import xml.etree.ElementTree as ET


def get_xml_entries(xml_f, child_l):
    """
    Get xml entries
    """
    cono = ET.parse(xml_f)
    root = cono.getroot()

    for entry in root:
        items = []
        for child in entry:

            if child.tag in child_l:
                items.append(child.text)

        yield items

def get_cono_nuc(cono_xml):
    """
    Print out a fasta file from conoserver_nucleic.xml
    """
    for name, seq in get_xml_entries(cono_xml,["name","sequence"]):
        print ">" + name.replace(" ", "_")
        print seq
#MAIN

import sys
get_cono_nuc(sys.argv[1])

