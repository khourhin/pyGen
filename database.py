#! /usr/bin/python

import os
import argparse
import datetime
import sqlite3
import sys
import basics_fasta as bfas
import work_on_b2go as wb2go

#-------------------------------------------------------------------------------
def createDB(dbName):

    with sqlite3.connect(dbName) as conn:
        c = conn.cursor()
        c.execute('CREATE TABLE jobs (jid INTEGER PRIMARY KEY AUTOINCREMENT, jname, date)')
        c.execute('CREATE TABLE annots (jid, seqid, genecode, EC, GOs, KEGGs)')
        c.execute('CREATE TABLE seqs (jid, seqid, seq, len)')
        conn.commit()

#-------------------------------------------------------------------------------
def createNewJob(jname, dbName):
    
        with sqlite3.connect(dbName) as conn:
            c = conn.cursor()
            
            c.execute('SELECT * FROM jobs WHERE jname=?', jname)
            res = c.fetchone()
            if res != None:
                raise IOError("The jobname '{}' already exists in the database".format(jname))

            now = datetime.datetime.now()
            c.execute('INSERT INTO jobs (jname, date) VALUES (?,?)', (jname, now))
            
#-------------------------------------------------------------------------------
def loadFasta(fas_f, jid, dbName):
    
    fas_d = bfas.fasta_2_dict(fas_f)

    with sqlite3.connect(dbName) as conn:
        c = conn.cursor()

        for seqid in fas_d:
            seq = fas_d[seqid]
            c.execute('INSERT INTO seqs (jid, seqid, seq, len) VALUES (?,?,?,?)',
                      (jid, seqid, seq, len(seq)))

#-------------------------------------------------------------------------------
def loadEnsemblAnnot(ensembl_f, jid, dbName):
    pass
    

#-------------------------------------------------------------------------------
def loadB2go(b2go_f, jid, db):
    """ NOT WORKING YET
    NOTE FOR LATER: have one tabke of each types of annots (ex: 1 GO in each row)
    """
    b2go_d = wb2go.parse_b2go(b2go_f)

    with sqlite3.connect(db) as conn:
        c = conn.cursor()

        for seqid in b2go_d:
            GOs =  ",".join(b2go_d[seqid]["GOs"])
            c.execute('INSERT INTO annots(seqid, GOs) VALUES(?, ?)', (seqid, GOs))

        conn.commit()

#-------------------------------------------------------------------------------
# MAIN
#-------------------------------------------------------------------------------
if __name__ == "__main__":

    DBNAME = "test.db"
    NEWJOB = "job1"
    FASTA = "/home/tiennou/bin/biopy_lib/demo_data/AP1.fas"

#    createDB( DBNAME )
#    createNewJob( NEWJOB, DBNAME )
#    load_fasta( FASTA, NEWJOB, DBNAME )

    parser = argparse.ArgumentParser("For creating and browsing NGS databases")
    parser.add_argument("db", help="The path to the db to create/update")
    parser.add_argument("-j","--jobname", type=str, help="An identifiers to identify various datasets in a same database")
    parser.add_argument("-f","--fasta", help="A simple multi-fasta file to feed the database with")
    parser.add_argument("-a", "--annots", help="NOT WORKING YET !! An annotation file with on seqid by lines.")

    args = parser.parse_args()

    if os.path.isfile(args.db):
        print "Database already exists"
    else:
        print "Creating database {}".format(args.db)
        createDB( args.db )

    if args.jobname:
        print "Creating job {0} in db {1}".format(args.jobname, args.db)
        createNewJob( args.jobname, args.db)

        if args.fasta:
            print "Adding sequences from {0} in db {1} with jobname {2}".format(args.fasta, args.db, args.jobname)
            loadFasta(args.fasta, args.jobname, args.db)
        if args.annots:
            # Only work with b2go outputs so far:
            loadB2go(args.annots, args.jobname, args.db)

