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

        self.lanes = defaultdict(list)
        self.lane_data = defaultdict(lambda: defaultdict(lambda: defaultdict(int)))
        self.lane_statuses = defaultdict(lambda: defaultdict(lambda: {'pass': 0, 'warn': 0, 'fail': 0}))
        self.lane_reads = defaultdict(int)
        self.lane_stats = defaultdict(lambda: defaultdict(float))
        self.undetermined_check = defaultdict(lambda: defaultdict(int))
        
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

        self.status_colours = {
            'pass': '#5cb85c',
            'warn': '#f0ad4e',
            'fail': '#d9534f',
            'default': '#999'
        }

        log.info("Found {} reports".format(len(self.fastqc_stats)))

        self.group_by_lane(s_names)
        self.total_lane_reads()
        self.aggregate_lane_statuses()
        self.aggregate_lane_data_n_content()
        self.aggregate_lane_data_sequence_quality()
        self.aggregate_lane_data_per_seq_quality()
        self.aggregate_lane_data_seq_dup_levels()
        self.aggregate_lane_stats()
        self.aggregate_undetermined_check()
        self.add_undetermined_check()

        self.write_data_file(self.lane_statuses, 'multiqc_fastqc_lane')
        self.sections = list()

        self.lane_stats_table()
        self.sequence_quality_plot()
        self.per_seq_quality_plot()
        self.sequence_content_plot()
        self.gc_content_plot()
        self.n_content_plot()
        self.seq_length_dist_plot()
        self.seq_dup_levels_plot()
        self.adapter_content_plot()


    def group_by_lane(self, s_names):
        for s_name in s_names:
            m = re.match('.*_(L00[1-8])_R[12]_[0-9]{3}$', s_name)
            self.lanes[m.group(1)].append(s_name)


    def total_lane_reads(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                self.lane_reads[lane] += self.fastqc_stats[s_name]['total_sequences']
    

    def aggregate_lane_data_n_content(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                for pos, value in self.fastqc_data['n_content'][s_name].items():
                    self.lane_data['n_content'][lane][pos] += value * \
                        self.fastqc_stats[s_name]['total_sequences']

            for pos in self.lane_data['n_content'][lane]:
                self.lane_data['n_content'][lane][pos] /= self.lane_reads[lane]


    def aggregate_lane_data_sequence_quality(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                for pos, value in self.fastqc_data['sequence_quality']['mean'][s_name].items():
                    self.lane_data['sequence_quality'][lane][pos] += value * \
                        self.fastqc_stats[s_name]['total_sequences']

            for pos in self.lane_data['sequence_quality'][lane]:
                self.lane_data['sequence_quality'][lane][pos] /= self.lane_reads[lane]


    def aggregate_lane_data_per_seq_quality(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                for qual, value in self.fastqc_data['per_seq_quality'][s_name].items():
                    self.lane_data['per_seq_quality'][lane][qual] += value

            for qual in self.lane_data['per_seq_quality'][lane]:
                self.lane_data['per_seq_quality'][lane][qual] = \
                    int(self.lane_data['per_seq_quality'][lane][qual] / \
                    len(self.lanes[lane]))


    def aggregate_lane_data_seq_dup_levels(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                for pos, value in self.fastqc_data['seq_dup_levels'][s_name].items():
                    self.lane_data['seq_dup_levels'][lane][pos] += value * \
                        self.fastqc_stats[s_name]['total_sequences']

            for pos in self.lane_data['seq_dup_levels'][lane]:
                self.lane_data['seq_dup_levels'][lane][pos] /= self.lane_reads[lane]


    def aggregate_lane_stats(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                self.lane_stats[lane]['percent_duplicates'] += \
                    self.fastqc_stats[s_name]['percent_duplicates'] * \
                    self.fastqc_stats[s_name]['total_sequences']
                self.lane_stats[lane]['percent_gc'] += \
                    self.fastqc_stats[s_name]['percent_gc'] * \
                    self.fastqc_stats[s_name]['total_sequences']
                self.lane_stats[lane]['avg_sequence_length'] += \
                    self.fastqc_stats[s_name]['avg_sequence_length'] * \
                    self.fastqc_stats[s_name]['total_sequences']

            self.lane_stats[lane]['total_sequences'] = self.lane_reads[lane]
            self.lane_stats[lane]['percent_duplicates'] /= self.lane_reads[lane]
            self.lane_stats[lane]['percent_gc'] /= self.lane_reads[lane]
            self.lane_stats[lane]['avg_sequence_length'] /= self.lane_reads[lane]
                

    def aggregate_lane_statuses(self):
        for metric, data in self.fastqc_statuses.items():
            for lane, s_names in self.lanes.items():
                for s_name in s_names:
                    self.lane_statuses[metric][lane][data[s_name]] += 1


    def lane_stats_table(self):
        """ Add some single-number stats to the basic statistics
        table at the top of the report """
        
        # Are sequence lengths interesting?
        seq_lengths = [x['avg_sequence_length'] for x in self.lane_stats.values()]
        seq_length_range = max(seq_lengths) - min(seq_lengths)
        
        headers = OrderedDict()
        headers['percent_duplicates'] = {
            'title': '% Dups',
            'description': '% Duplicate',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'RdYlGn-rev',
            'format': '{:.1f}%'
        }
        headers['percent_gc'] = {
            'title': '% GC',
            'description': 'Average % GC Content',
            'max': 100,
            'min': 0,
            'suffix': '%',
            'scale': 'Set1',
            'format': '{:.0f}%'
        }
        headers['avg_sequence_length'] = {
            'title': 'Length',
            'description': 'Average Sequence Length (bp)',
            'min': 0,
            'suffix': 'bp',
            'scale': 'RdYlGn',
            'format': '{:.0f}',
            'hidden': False if seq_length_range > 10 else True
        }
        headers['total_sequences'] = {
            'title': 'M Seqs',
            'description': 'Total Sequences in Lane (millions)',
            'min': 0,
            'scale': 'Blues',
            'modify': lambda x: x / 1000000,
            'shared_key': 'read_count'
        }

        pconfig = { 
            'id': 'general_stats_table_by_lane',
            'table_title': 'general statistics by Lane',
            'save_file': True,
            'raw_data_fn':'multiqc_general_stats_by_lane'
        }   

        self.sections.append({
            'name': 'General Statistics by Lane',
            'anchor': 'fastqc_by_lane_general_stats',
            'content': plots.table.plot(self.lane_stats, headers, pconfig)
        })


    def statuses_plot(self, metric, title):
        data = self.lane_statuses[metric]

        colors = {
            'pass': {
                'color': '#5cb85c',
                'name': 'pass',
            },
            'warn': {
                'color': '#f0ad4e',
                'name': 'warn',
            },
            'fail': {
                'color': '#d9534f',
                'name': 'fail',
            },
        }

        pconfig = {
            'title': title,
            'cpswitch': False,
        }

        return plots.bargraph.plot(data, colors, pconfig)


    def sequence_quality_plot(self):
        name = 'Sequence Quality Histograms by Lane'
        metric = 'sequence_quality'

        if metric not in self.fastqc_data or not self.fastqc_data[metric]:
            log.debug('sequence_quality not found in FastQC reports')
            return None    
        
        pconfig = {
            'id': 'fastqc_by_lane_sequence_quality_plot',
            'title': 'Mean Quality Scores by Lane',
            'ylab': 'Phred Score',
            'xlab': 'Position (bp)',
            'ymin': 0,
            'xDecimals': False,
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}',
            'colors': self.get_status_cols('sequence_quality'),
            'yPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.sections.append({
            'name': name,
            'anchor': 'fastqc_by_lane_' + metric,
            'content': self.statuses_plot(metric, name) + \
                       plots.linegraph.plot(self.lane_data[metric], pconfig)
        })


    def per_seq_quality_plot(self):
        name = 'Per Sequence Quality Scores by Lane'
        metric = 'per_seq_quality'

        if metric not in self.fastqc_data or not self.fastqc_data[metric]:
            log.debug('per_seq_quality not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_by_lane_per_seq_quality_plot',
            'title': name,
            'ylab': 'Count',
            'xlab': 'Mean Sequence Quality (Phred Score)',
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols(metric),
            'tt_label': '<b>Phred {point.x}</b>: {point.y} reads',
            'xPlotBands': [
                {'from': 28, 'to': 100, 'color': '#c3e6c3'},
                {'from': 20, 'to': 28, 'color': '#e6dcc3'},
                {'from': 0, 'to': 20, 'color': '#e6c3c3'},
            ]
        }
        self.sections.append({
            'name': name,
            'anchor': 'fastqc_by_lane_' + metric,
            'content': self.statuses_plot(metric, name) + \
                       plots.linegraph.plot(self.lane_data[metric], pconfig)
        })


    def sequence_content_plot(self):
        name = 'Per Base Sequence Content by Lane'
        metric = 'sequence_content'

        self.sections.append({
            'name': name,
            'anchor': metric,
            'content': self.statuses_plot(metric, name)
        })


    def gc_content_plot(self):
        name = 'Per Sequence GC Content by Lane'
        metric = 'sequence_content'
        
        self.sections.append({
            'name': name,
            'anchor': metric,
            'content': self.statuses_plot(metric, name)
        })
    

    def n_content_plot(self):
        name = 'Per Base N Content by Lane'
        metric = 'n_content'

        if 'n_content' not in self.lane_data or not self.lane_data['n_content']:
            log.debug('n_content not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_by_lane_n_content_plot',
            'title': 'Per Base N Content by Lane',
            'ylab': 'Percentage N-Count',
            'xlab': 'Position in Read (bp)',
            'yCeiling': 100,
            'yMinRange': 5,
            'ymin': 0,
            'xmin': 0,
            'xDecimals': False,
            'colors': self.get_status_cols('n_content'),
            'tt_label': '<b>Base {point.x}</b>: {point.y:.2f}%',
            'yPlotBands': [
                {'from': 20, 'to': 100, 'color': '#e6c3c3'},
                {'from': 5, 'to': 20, 'color': '#e6dcc3'},
                {'from': 0, 'to': 5, 'color': '#c3e6c3'},
            ]
        }

        self.sections.append({
            'name': name,
            'anchor': 'fastqc_by_lane_' + metric,
            'content': self.statuses_plot(metric, name) + \
                       plots.linegraph.plot(self.lane_data[metric], pconfig)
        })
    

    def seq_length_dist_plot(self):
        name = 'Sequence Length Distribution by Lane'
        metric = 'seq_length_dist'

        self.sections.append({
            'name': name,
            'anchor': metric,
            'content': self.statuses_plot(metric, name)
        })


    def seq_dup_levels_plot(self):
        name = 'Sequence Duplication Levels by Lane'
        metric = 'seq_dup_levels'

        if metric not in self.fastqc_data or not self.fastqc_data[metric]:
            log.debug('seq_dup_levels not found in FastQC reports')
            return None
        
        pconfig = {
            'id': 'fastqc_by_lane_seq_dup_levels_plot',
            'title': name,
            'categories': True,
            'ylab': '% of Library',
            'xlab': 'Sequence Duplication Level',
            'ymax': 100,
            'ymin': 0,
            'yMinTickInterval': 0.1,
            'colors': self.get_status_cols(metric),
            'tt_label': '<b>{point.x}</b>: {point.y:.1f}%',
        }
    
        self.sections.append({
            'name': name,
            'anchor': 'fastqc_by_lane_' + metric,
            'content': self.statuses_plot(metric, name) + \
                       plots.linegraph.plot(self.lane_data[metric], pconfig)
        })


    def adapter_content_plot(self):
        name = 'Adapter Content by Lane'
        metric = 'adapter_content'

        self.sections.append({
            'name': name,
            'anchor': metric,
            'content': self.statuses_plot(metric, name)
        })


    def get_status_cols(self, metric):
        """ Helper function - returns a list of colours according to the FastQC
        status of this module for each sample. """
        colours = dict()
        highest = 0
        final_status = 'pass'

        for lane, statuses in self.lane_statuses[metric].items():
            for status, count in statuses.items():
                if count > highest:
                    highest = count
                    final_status = status

            colours[lane] = self.status_colours.get(final_status, self.status_colours['default'])

        return colours


    def aggregate_undetermined_check(self):
        for lane, s_names in self.lanes.items():
            for s_name in s_names:
                if s_name.startswith('Undetermined_'):
                    undet_name = s_name

            for s_name in s_names:
                if s_name == undet_name:
                    self.undetermined_check[s_name]['undetermined_check'] = 0
                    continue

                if self.fastqc_stats[s_name]['total_sequences'] < \
                        self.fastqc_stats[undet_name]['total_sequences']:
                    self.undetermined_check[s_name]['undetermined_check'] = 2
                else:
                    self.undetermined_check[s_name]['undetermined_check'] = 1

                if self.fastqc_stats[s_name]['total_sequences'] < 1000000:
                    self.undetermined_check[s_name]['undetermined_check'] += 2


    def add_undetermined_check(self):
        headers = OrderedDict()
        headers['undetermined_check'] = {
            'title': 'Seqs < Undet',
            'description': 'Number of sequences lower than number of undetermined sequences of this lane: 1 = OK, 2 = seq < undet, 3 = seq < 1.000.000, 4 = 2 & 3',
            'min': 0,
            'max': 4,
            'scale': 'Paired',
        }

        self.general_stats_addcols(self.undetermined_check, headers)
