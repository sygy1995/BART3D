# Time-stamp: <2017-08-10>
'''
Copyright (c) 2017, 2018 Chongzhi Zang, Zhenjia Wang <zhenjia@virginia.edu>

This code is free software; you can redistribute it and/or modify it 
under the terms of the BSD License.

@status: release candidate
@version: $Id$
@author: Chongzhi Zang, Zhenjia Wang
@contact: zhenjia@virginia.edu

'''

import configparser
import os

script_dir = os.path.dirname(os.path.realpath(__file__))

def conf_validate():
    '''
    Read user provided path from 'bart.conf' config file
    '''
    config = configparser.ConfigParser()
    config_path = os.path.join(script_dir, 'bart.conf')
    if not os.path.exists(config_path):
        sys.stderr.write("CRITICAL: bart.conf does not exist in {}!\n".format(script_dir))
        sys.exit(1)

    config.read(config_path)
    return config

def opt_validate(options):
    '''
    Validate input options and specify used data.
    '''
    config = conf_validate()

    if not options.outdir or options.outdir=='BART3D_output/':
        options.outdir = os.path.join(os.getcwd(), 'BART3D_output') # create output directory at current working directory

    # === hg38 ===
    if options.species == 'hg38':   
        if config['path']['hg38_library_dir']:
            data_dir = os.path.join(config['path']['hg38_library_dir'], 'hg38_library')
        else:
            data_dir = os.path.join(script_dir, 'hg38_library')
        print("Library directory:" + data_dir)

        # file for standardization

        options.normfile = data_dir+os.sep+'bart2_hg38_H3K27ac.dat'

    # === mm10 ===
    elif options.species == 'mm10': 
        if config['path']['mm10_library_dir']:
            data_dir = os.path.join(config['path']['mm10_library_dir'], 'mm10_library')
        else:
            data_dir = os.path.join(script_dir, 'mm10_library')
        print("Library directory:" + data_dir)
        
        # file for standardization

        options.normfile = os.path.join(data_dir, 'bart2_mm10_H3K27ac.dat')

    # used for generating cis-regulatory profile
    options.rp = data_dir+os.sep+'bart2_'+options.species+'_RP.h5'
    options.rpkm = data_dir+os.sep+'bart2_'+options.species+'_UDHS_H3K27ac.h5'
    options.tss = data_dir+os.sep+'bart2_'+options.species+'_refseq_TSS.bed'
    options.dhsfile = data_dir+os.sep+'bart2_'+options.species+'_UDHS.bed'
    options.desc = data_dir+os.sep+'bart2_'+options.species+'_H3K27ac_Description.dat'

    # used for BART AUC calculation
    options.tffile = data_dir+os.sep+'bart2_'+options.species+'_TF_file.json'
    options.tfoverlap = data_dir+os.sep+'bart2_'+options.species+'_TF_overlap.json'

    return options

