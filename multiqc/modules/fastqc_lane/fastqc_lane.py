#!/usr/bin/env python

""" MultiQC module to parse output from FastQC
"""

######################################################
#### LOOKING FOR AN EXAMPLE OF HOW MULTIQC WORKS? ####
######################################################
#### Stop! This module is huge and complicated.   ####
#### Have a look at Bowtie or STAR for a simpler  ####
#### example. CONTRIBUTING.md has documentation.  ####
######################################################

from __future__ import print_function
from collections import defaultdict, OrderedDict
import json
import logging
import os
import re
import zipfile
import re

#from pprint import pprint

from multiqc import config, plots
from multiqc.modules import fastqc

# Initialise the logger
log = logging.getLogger(__name__)

class MultiqcModule(fastqc.MultiqcModule):

    def __init__(self):

        # Initialise the parent object
        super(fastqc.MultiqcModule, self).__init__(name='FastQC by Lane', anchor='fastqc-by-lane', 
        href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/", 
        info="groups FastQC reports by lane.")

        self.fastqc_data = defaultdict(lambda: defaultdict(lambda: defaultdict( lambda: defaultdict() ) ) )
        self.fastqc_stats = dict()
        self.fastqc_statuses = defaultdict(lambda: defaultdict())
        
        # Find and parse unzipped FastQC reports
        for f in self.find_log_files(config.sp['fastqc']['data']):
            s_name = self.clean_s_name(os.path.basename(f['root']), os.path.dirname(f['root']))
            self.parse_fastqc_report(f['f'], s_name, f)

        # Find and parse zipped FastQC reports
        s_names = []
        for f in self.find_log_files(config.sp['fastqc']['zip'], filecontents=False):
            s_name = f['fn']
            if s_name.endswith('_fastqc.zip'):
                s_name = s_name[:-11]
            try:
                fqc_zip = zipfile.ZipFile(os.path.join(f['root'], f['fn']))
            except zipfile.BadZipfile:
                continue
            # FastQC zip files should have just one directory inside, containing report
            d_name = fqc_zip.namelist()[0]
            try:
                with fqc_zip.open(os.path.join(d_name, 'fastqc_data.txt')) as fh:
                    r_data = fh.read().decode('utf8')
                    self.parse_fastqc_report(r_data, s_name, f)
            except KeyError:
                log.warning("Error - can't find fastqc_raw_data.txt in {}".format(f))
            s_names.append(s_name)

        if len(self.fastqc_stats) == 0:
            log.debug("Could not find any reports in {}".format(config.analysis_dir))
            raise UserWarning

        log.info("Found {} reports".format(len(self.fastqc_stats)))

        self.status_colours = {
            'pass': '#5cb85c',
            'warn': '#f0ad4e',
            'fail': '#d9534f',
            'default': '#999'
        }

        self.group_by_lane(s_names)
        self.add_lane_statuses()

        #print("Lane statuses")
       # pprint(self.lane_statuses)
        #print("Lane statuses colors")
        #pprint(self.lane_statuses_colors)

        # Add to self.css and self.js to be included in template
#        self.css = {
#            'assets/css/multiqc_fastqc.css' : os.path.realpath(os.path.join(os.path.dirname(__file__), 'assets', 'css', 'multiqc_fastqc.css'))
#        }
#        self.js = {
#            'assets/js/multiqc_fastqc.js' : os.path.realpath(os.path.join(os.path.dirname(__file__), 'assets', 'js', 'multiqc_fastqc.js'))
#        }
#        
#        # Add to the general statistics table
#        self.fastqc_stats_table()
#        
#        # Write the basic stats table data to a file
        self.write_data_file(self.fastqc_statuses, 'multiqc_fastqc_lane')
        self.sections = list()

        self.sequence_quality_plot()
        self.per_seq_quality_plot()
        self.sequence_content_plot()
        self.gc_content_plot()
        self.n_content_plot()
        self.seq_length_dist_plot()
        self.seq_dup_levels_plot()
        self.adapter_content_plot()


    def group_by_lane(self, s_names):
        self.lanes = defaultdict(list)
        for s_name in s_names:
            m = re.match('.*_(L00[1-8])_R[12]_[0-9]{3}$', s_name)
            self.lanes[m.group(1)].append(s_name)


    # dummy function, fill in future if needed
    def add_lane_data(self):
        for metric, values in self.fastqc_data.items():
            print(metric)
            for lane, s_names in self.lanes.items():
                print(lane)
                for s_name in s_names:
                    print(s_name)
                    pprint(values[s_name])
    

    # dummy function, fill in future if needed
    def add_lane_stats(self):
        for lane, s_names in self.lanes.items():
            print(lane)
            for s_name in s_names:
                self.fastqc_stats[s_name]
                # this won't work because there is some additional data in key string in some metrics"
                

    def add_lane_statuses(self):
        self.lane_statuses = defaultdict(lambda: defaultdict(lambda: {'pass': 0, 'warn': 0, 'fail': 0}))
        self.lane_statuses_colors = defaultdict(lambda: defaultdict())

        for metric, data in self.fastqc_statuses.items():
            for lane, s_names in self.lanes.items():
                for s_name in s_names:
                    self.lane_statuses[metric][lane][data[s_name]] += 1
                        
                self.lane_statuses_colors[metric][data[s_name]] = {
                    'color': self.status_colours[data[s_name]],
                    'name': data[s_name],
                }

    def statuses_plot(self, metric, title):
        colors = self.lane_statuses_colors[metric]
        data = self.lane_statuses[metric]

        pconfig = {
            'title': title,
            'cpswitch': False,
        }

        self.sections.append({
            'name': title,
            'anchor': metric,
            'content': plots.bargraph.plot(data, colors, pconfig)
        })

    def sequence_quality_plot(self):
        self.statuses_plot('sequence_quality', 'Sequence Quality Histograms by Lane')

    def per_seq_quality_plot(self):
        self.statuses_plot('per_seq_quality', 'Per Sequence Quality Scores by Lane')

    def sequence_content_plot(self):
        self.statuses_plot('sequence_content', 'Per Base Sequence Content by Lane')

    def gc_content_plot(self):
        self.statuses_plot('gc_content', 'Per Sequence GC Content by Lane')
    
    def n_content_plot(self):
        self.statuses_plot('n_content', 'Per Base N Content by Lane')
    
    def seq_length_dist_plot(self):
        self.statuses_plot('seq_length_dist', 'Sequence Length Distribution by Lane')

    def seq_dup_levels_plot(self):
        self.statuses_plot('seq_dup_levels', 'Sequence Duplication Levels by Lane')

    def adapter_content_plot(self):
        self.statuses_plot('adapter_content', 'Adapter Content by Lane')
