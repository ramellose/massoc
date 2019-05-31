"""
The process panel allows users to format files. They can split data by cluster ID or sample features,
or set a prevalence filter.
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import wx
from wx.lib.pubsub import pub
import numpy as np
import biom
from scipy.sparse import csr_matrix
from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.figure import Figure
import sys
import logging
import os
import logging.handlers

logger = logging.getLogger()
logger.setLevel(logging.INFO)

# handler to sys.stdout
sh = logging.StreamHandler(sys.stdout)
sh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
sh.setFormatter(formatter)
logger.addHandler(sh)


class ProcessPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.frame = parent
        # subscribe to input tab
        pub.subscribe(self.load_settings, 'load_settings')
        pub.subscribe(self.enable_tax, 'receive_tax')
        pub.subscribe(self.toggle_network, 'toggle_network')
        pub.subscribe(self.set_settings, 'input_settings')
        pub.subscribe(self.set_meta, 'input_metadata')


        btnsize = (300, -1)
        btnmargin = 10
        boxsize = (300, 50)
        # adds columns

        self.settings = dict()
        self.name = None
        self.prev = 20
        self.split = None
        self.cluster = None
        self.agglom = None
        self.nclust = None
        self.rar = None
        self.min = None
        self.filelist = None
        self.meta = None
        # metadata only required for setting split list
        self.topsizer = wx.BoxSizer(wx.HORIZONTAL)

        # defines columns
        self.leftsizer = wx.BoxSizer(wx.VERTICAL)
        self.topleftsizer = wx.BoxSizer(wx.VERTICAL)
        self.bottomleftsizer = wx.BoxSizer(wx.VERTICAL)

        # include a column to report errors or check file formats
        self.rightsizer = wx.BoxSizer(wx.VERTICAL)

        # set savefile prefix
        self.prefix_txt = wx.StaticText(self, label='Prefix for saved files')
        self.prefix = wx.TextCtrl(self, value='', size=btnsize)
        self.prefix.Bind(wx.EVT_MOTION, self.update_help)
        self.prefix.Bind(wx.EVT_TEXT, self.get_filename)

        # diagnostics button
        self.dg_button = wx.StaticText(self, label='Show data properties', size=btnsize)
        self.file_list = wx.ListBox(self, choices=['Select files first'], size=btnsize)
        self.file_list.Bind(wx.EVT_LISTBOX, self.register_figures)
        self.file_list.Bind(wx.EVT_MOTION, self.update_help)
        self.figure1 = Figure(figsize=(3, 2))
        self.prevfig = self.figure1.add_subplot(211)
        self.prevfig.set_xlabel('Prevalence')
        self.prevfig.set_title('Taxon prevalence')
        self.prevfig.set_ylabel('Number of taxa')
        self.figure1.set_tight_layout(True)
        self.figure2 = Figure(figsize=(3, 2))
        self.rarfig = self.figure2.add_subplot(211)
        self.rarfig.set_xlabel('Count number')
        self.rarfig.set_title('Sample counts')
        self.rarfig.set_ylabel('Number of samples')
        self.figure2.set_tight_layout(True)
        self.canvas1 = FigureCanvas(self, -1, self.figure1)
        self.canvas2 = FigureCanvas(self, -1, self.figure2)

        # set minimum count
        self.min_txt = wx.StaticText(self, label='Remove taxa with mean count below:')
        self.min_number = wx.TextCtrl(self, value='0', size=btnsize)
        self.min_number.Bind(wx.EVT_TEXT, self.update_min)
        self.min_number.Bind(wx.EVT_MOTION, self.update_help)

        # set rarefaction
        self.rar_txt = wx.StaticText(self, label='Rarefy:')
        self.rar_choice = wx.ListBox(self, choices=['To even depth', 'To count number'], size=(300, 50))
        self.rar_choice.Bind(wx.EVT_LISTBOX, self.get_rarefaction)
        self.rar_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.rar_number = wx.TextCtrl(self, value='', size=btnsize)
        self.rar_number.Bind(wx.EVT_TEXT, self.update_rarefaction)
        self.rar_number.Hide()

        # set prevalence filter
        self.prev_txt = wx.StaticText(self, label='Set prevalence filter in %')
        self.prev_val = wx.TextCtrl(self, value='20', size=btnsize)
        self.prev_slider = wx.Slider(self, value=20, minValue=0, maxValue=100, size=btnsize, style=wx.SL_AUTOTICKS)
        self.prev_slider.SetTickFreq(10)
        self.prev_val.Bind(wx.EVT_TEXT, self.update_slider)
        self.prev_val.Bind(wx.EVT_MOTION, self.update_help)
        self.prev_slider.Bind(wx.EVT_SLIDER, self.update_text)
        self.prev_slider.Bind(wx.EVT_MOTION, self.update_help)

        # split samples by
        self.split_txt = wx.StaticText(self, label='Split samples by metadata variable')
        self.split_list = wx.ListBox(self, choices=['Select a BIOM or metadata file first'], size=(300, 100))
        self.split_list.Bind(wx.EVT_MOTION, self.update_help)
        self.split_list.Bind(wx.EVT_LISTBOX, self.split_files)

        # cluster samples
        self.cluster_btn = wx.Button(self, label='Show clustering dialog', size=btnsize)
        self.cluster_btn.Bind(wx.EVT_MOTION, self.update_help)
        self.cluster_btn.Bind(wx.EVT_BUTTON, self.show_dialog)

        # select taxonomic levels
        self.tax_txt = wx.StaticText(self, label='Taxonomic levels to analyze')
        self.tax_choice = wx.ListBox(self, choices=['OTU', 'Species', 'Genus',
                                     'Family', 'Order', 'Class', 'Phylum'],
                                     size=(boxsize[0], 150), style=wx.LB_MULTIPLE)
        self.tax_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.tax_choice.Bind(wx.EVT_LISTBOX, self.get_levels)
        self.tax_choice.Enable(False)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.file_list, 0, wx.EXPAND | wx.ALL, 10)

        self.topleftsizer.Add(self.prefix_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.prefix, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(10)
        self.topleftsizer.Add(self.dg_button, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.file_list, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(10)
        self.topleftsizer.Add(self.canvas1, 0, wx.EXPAND | wx.ALL, 10)
        self.topleftsizer.Add(self.canvas2, 0, wx.EXPAND | wx.ALL, 10)
        self.topleftsizer.AddSpacer(30)

        self.rightsizer.Add(self.min_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.min_number, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(10)
        self.rightsizer.Add(self.rar_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.rar_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.rar_number, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(10)
        self.rightsizer.Add(self.prev_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.prev_val, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.prev_slider, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(10)
        self.rightsizer.Add(self.split_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.split_list, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(10)
        self.rightsizer.Add(self.tax_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.tax_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(10)
        self.rightsizer.Add(self.cluster_btn, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.leftsizer.Add(self.topleftsizer, wx.ALL, 20)
        self.leftsizer.AddSpacer(10)
        self.rightsizer.Add(self.bottomleftsizer, wx.ALL, 20)
        self.bottomleftsizer.ShowItems(show=False)
        self.topsizer.Add(self.leftsizer, 0, wx.ALL, 20)
        self.topsizer.Add(self.rightsizer, 0, wx.ALL, 20)
        self.SetSizerAndFit(self.topsizer)
        self.Fit()

        # help strings for buttons
        self.buttons = {self.prev_val: 'With a prevalence filter of 20%, only taxa present in 20% of '
                                       'samples are retained. ',
                        self.prefix: 'This prefix is used to generate filenames. \n'
                                     'If you have more than one BIOM file, separate by ; .',
                        self.prev_slider: 'With a prevalence filter of 0.2, only taxa present in 20% '
                                       'of samples are retained. ',
                        self.split_list: 'Separate a file according to sample metadata and run network analysis on '
                                         'both of the split files. ',
                        self.cluster_btn: 'Cluster samples according to different algorithms; either use cluster '
                                          'identities to split files, or include them as additional metadata. ',
                        self.rar_choice: 'Specify whether samples should be rarefied to equal depth or a '
                                         'specific value. Samples below a specified threshold are removed.',
                        self.min_number: 'Taxa with low mean counts can be removed if there are too many '
                                         'high-prevalence, '
                                         'low-abundance taxa.',
                        self.file_list: 'Show taxon prevalence and counts per sample for imported count tables.',
                        self.tax_choice: 'Select which taxonomic level to agglomerate to before network inference.'}

    def update_help(self, event):
        btn = event.GetEventObject()
        if btn in self.buttons:
            status = self.buttons[btn]
            pub.sendMessage('change_statusbar', msg=status)

    def get_filename(self, event):
        text = self.prefix.GetValue()
        self.name = text.split(';')
        self.name = [x.strip() for x in self.name]
        self.send_settings()

    def get_rarefaction(self, event):
        rar = self.rar_choice.GetSelection()
        text = self.rar_choice.GetString(rar)
        if text == 'To even depth':
            self.rar = 'True'
        if text == 'To count number':
            self.rar_number.Show()
            self.Layout()
        self.send_settings()

    def update_rarefaction(self, event):
        value = self.rar_number.GetValue()
        self.rar = value

    def update_min(self, event):
        value = self.min_number.GetValue()
        if len(value) > 0:
            self.min = int(value)

    def enable_tax(self, msg):
        self.tax_choice.Enable(True)

    def split_files(self, event):
        name = self.split_list.GetSelection()
        text = self.split_list.GetString(name)
        self.split = text
        self.send_settings()

    def get_levels(self, event):
        text = list()
        try:
            ids = self.tax_choice.GetSelections()
            for i in ids:
                text.append(self.tax_choice.GetString(i))
            text = [x.lower() for x in text]
            self.agglom = text
        except KeyError:
            pass
        self.send_settings()

    def update_slider(self, event):
        value = self.prev_val.GetValue()
        try:
            value = float(value)
            self.prev_slider.SetValue(value)
            self.prev = float(value)
            self.send_settings()
        except KeyError:
            pass

    def update_text(self, event):
        value = self.prev_slider.GetValue()
        try:
            value = str(value)
            self.prev_val.SetValue(value)
            self.prev = float(value)
        except KeyError:
            pass

    def show_dialog(self, event):
        dlg = ClusterDialog(self.filelist)
        dlg.ShowModal()

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'prev': self.prev, 'split': self.split, 'name': self.name,
                    'levels': self.agglom, 'nclust': self.nclust,
                    'rar': self.rar, 'min': self.min, 'cluster': self.cluster}
        pub.sendMessage('process_settings', msg=settings)

    def set_settings(self, msg):
        """
        Stores settings file as tab property so it can be read by save_settings.
        """
        self.settings = msg
        self.filelist = list()
        try:
            if self.settings['biom_file'] is not None:
                for name in self.settings['biom_file']:
                    self.filelist.append(name)
            if self.settings['otu_table'] is not None:
                for name in self.settings['otu_table']:
                    self.filelist.append(name)
            if len(self.filelist) > 0:
                self.file_list.Set(self.filelist)
                self.file_list.SetSelection(0)
                self.generate_figures()
        except Exception:
            logger.error("Failed to save settings. ", exc_info=True)

    def set_meta(self, msg):
        """Stores a dictionary of BIOM filenames and metadata vars."""
        self.meta = msg[0]
        self.split = msg[1]
        self.generate_figures()

    def load_settings(self, msg):
        """
        Listener function that changes input values
        to values specified in settings file.
        """
        # try
        if msg['min'] is not None:
            self.min = msg['min']
            self.min_number.SetValue(str(msg['min']))
        else:
            self.min = None
            self.min_number.SetValue('')
        if msg['name'] is not None:
            self.name = msg['name']
            self.prefix.SetValue("; ".join(msg['name']))
        else:
            self.name = None
            self.prefix.SetValue('')
            self.file_list.Set(['Select files first'])
            self.figure1.clf()
            self.figure2.clf()
        if msg['rar'] is not None:
            self.rar = msg['rar']
            rar = msg['rar']
            try:
                rar = int(rar)
                self.rar_number.SetValue(str(rar))
                self.rar_choice.SetSelection(1)
            except ValueError:
                self.rar_choice.SetSelection(0)
        else:
            self.rar = None
            choice = self.rar_choice.GetSelection()
            self.rar_choice.Deselect(choice)
        if msg['prev'] is not None:
            self.prev = msg['prev']
            self.prev_slider.SetTick(msg['prev'])
            self.prev_val.SetValue(str(msg['prev']))
        else:
            self.prev = 20
            self.prev_slider.SetTick(20)
            self.prev_val.SetValue(str(20))
        if msg['split'] is not None:
            self.split = msg['split']
        else:
            self.split = None
            self.split_list.Set(['Select a BIOM or metadata file first'])
        if msg['levels'] is not None:
            self.agglom = msg['levels']
            agglomdict = {'otu': 0, 'species': 1, 'genus': 2,
                          'family': 3, 'order': 4, 'class': 5,
                          'phylum': 6}
            for tax in msg['levels']:
                self.tax_choice.SetSelection(agglomdict[tax])
        else:
            self.agglom = None
            choice = self.tax_choice.GetSelections()
            for selection in choice:
                self.tax_choice.Deselect(selection)
        if msg['cluster'] is not None:
            self.cluster = msg['cluster']
            self.nclust = msg['nclust']
        else:
            self.cluster = None
            self.nclust = None
        # except Exception:
        #     logger.error("Unable to load settings", exc_info=True)
        self.send_settings()

    def register_figures(self, event):
        """Registers listbox event and calls generate_figures."""
        self.generate_figures()

    def generate_figures(self):
        """Generates figures for diagnostics canvas.
        Also sets the split file params. """
        file = self.file_list.GetSelection()
        if file != -1:
            file = self.file_list.GetString(file)
            biomfile = biom.load_table(file)
            if biomfile.metadata(axis='sample'):
                varlist = list(biomfile.metadata_to_dataframe(axis='sample').columns)
                varlist.sort()
                self.split_list.Set(varlist)
            else:
                if self.meta:
                    if file in self.meta:
                        varlist = self.meta[file]
                        varlist.sort()
                        self.split_list.Set(varlist)
            if self.split:
                split = self.split_list.FindString(self.split)
                self.split_list.SetSelection(split)
            data = biomfile.matrix_data
            data = csr_matrix.todense(data)
            fracs = np.count_nonzero(data, axis=1)
            nsamples = data.shape[1]
            fracs = fracs / nsamples
            self.prevfig.clear()
            self.prevfig.hist(fracs, bins=20)
            self.prevfig.set_xlabel('Prevalence')
            self.prevfig.set_title('Taxon prevalence')
            self.prevfig.set_ylabel('Number of taxa')
            sample_sums = np.transpose(np.count_nonzero(data, axis=0))
            self.rarfig.clear()
            self.rarfig.hist(sample_sums, bins=40)
            self.rarfig.set_xlabel('Count number')
            self.rarfig.set_title('Sample counts')
            self.rarfig.set_ylabel('Number of samples')
            self.canvas1.draw()
            self.canvas2.draw()



    def toggle_network(self, msg):
        if msg == 'Yes':
            self.min_txt.Enable(True)
            self.min_number.Enable(True)
            self.rar_txt.Enable(True)
            self.rar_choice.Enable(True)
            self.prev_txt.Enable(True)
            self.prev_val.Enable(True)
            self.prev_slider.Enable(True)
            self.split_txt.Enable(True)
            self.split_list.Enable(True)
            self.cluster_btn.Enable(True)
            self.tax_txt.Enable(True)
            self.tax_choice.Enable(True)
        if msg == 'No':
            self.min_txt.Enable(False)
            self.min_number.Enable(False)
            self.rar_txt.Enable(False)
            self.rar_choice.Enable(False)
            self.prev_txt.Enable(False)
            self.prev_val.Enable(False)
            self.prev_slider.Enable(False)
            self.split_txt.Enable(False)
            self.split_list.Enable(False)
            self.cluster_btn.Enable(False)
            self.tax_txt.Enable(False)
            self.tax_choice.Enable(False)


class ClusterDialog(wx.Dialog):
    def __init__(self, files):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Cluster samples", size=(700,500))
        self.cluster = None
        self.split = None
        pub.subscribe(self.load_settings, 'load_settings')

        self.alg_txt = wx.StaticText(self, label='Clustering algorithm')
        self.cluster_choice = wx.ListBox(self, choices=['K-means', 'DBSCAN', 'Gaussian', 'Spectral',
                                                        'Affinity Propagation', 'None'])
        self.act_txt = wx.StaticText(self, label='Use cluster ID to: ')
        self.cluster_proc = wx.ListBox(self, choices=['Add to metadata', 'Split files'])

        # cluster plots
        self.dg_button = wx.StaticText(self, label='Show preview of sample clustering')
        self.file_list = wx.ListBox(self, choices=files)
        self.file_list.SetSelection(0)
        self.figure1 = Figure(figsize=(4,3))
        self.prev = self.figure1.add_subplot(111)
        self.prev.set_xlabel('PCA axis 1')
        self.prev.set_title('PCA plot of CLR-transformed samples')
        self.prev.set_ylabel('PCA axis 2')
        self.figure1.set_tight_layout(True)
        self.canvas1 = FigureCanvas(self, -1, self.figure1)

        # choose to keep clusters
        self.ok_txt = wx.StaticText(self, label='Proceed with clustering?')
        self.ok_btn = wx.Button(self, label='OK')
        self.ok_btn.Bind(wx.EVT_BUTTON, self.cluster_func)
        self.no_btn = wx.Button(self, wx.ID_CANCEL, label='Cancel')

        self.topsizer = wx.BoxSizer(wx.HORIZONTAL)
        self.leftsizer = wx.BoxSizer(wx.VERTICAL)
        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.alg_txt, wx.ALIGN_CENTER_HORIZONTAL)
        self.leftsizer.Add(self.cluster_choice, wx.EXPAND | wx.ALL, 20)
        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.act_txt, wx.ALIGN_CENTER_HORIZONTAL)
        self.leftsizer.Add(self.cluster_proc, wx.EXPAND | wx.ALL, 20)
        self.leftsizer.AddSpacer(20)
        self.rightsizer = wx.BoxSizer(wx.VERTICAL)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.dg_button, wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.file_list, wx.EXPAND | wx.ALL, 20)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.canvas1, wx.EXPAND | wx.ALL, 20)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.ok_txt, wx.ALIGN_CENTER_HORIZONTAL)
        self.minisizer = wx.BoxSizer(wx.HORIZONTAL)
        self.minisizer.Add(self.ok_btn, wx.EXPAND | wx.ALL, 5)
        self.minisizer.Add(self.no_btn, wx.EXPAND | wx.ALL, 5)
        self.rightsizer.Add(self.minisizer)
        self.rightsizer.AddSpacer(20)
        self.topsizer.AddSpacer(20)
        self.topsizer.Add(self.leftsizer, wx.ALL, 20)
        self.topsizer.AddSpacer(10)
        self.topsizer.Add(self.rightsizer, wx.ALL, 20)
        self.topsizer.AddSpacer(20)
        self.SetSizer(self.topsizer)

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'split': self.split, 'cluster': self.cluster,
                    }
        pub.sendMessage('update_settings', msg=settings)

    def set_settings(self, msg):
        """
        Stores settings file as tab property so it can be read by save_settings.
        """
        self.settings = msg
        self.filelist = list()
        try:
            if self.settings['biom_file'] is not None:
                for name in self.settings['biom_file']:
                    self.filelist.append(name)
            if self.settings['otu_table'] is not None:
                for name in self.settings['otu_table']:
                    self.filelist.append(name)
            if len(self.filelist) > 0:
                self.file_list.Set(self.filelist)
                self.file_list.SetSelection(0)
                self.generate_figures()
        except Exception:
            logger.error("Failed to save settings. ", exc_info=True)

    def cluster_func(self, event):
        """
        This function specifies whether the cluster ID should be used to split files.
        """
        try:
            clus = self.cluster_choice.GetSelection()
            text = self.cluster_choice.GetString(clus)
            self.cluster = list()
            self.cluster.append(text)
        except KeyError:
            pass
        try:
            proc = self.cluster_proc.GetSelections()
            text = list()
            for i in proc:
                text.append(self.cluster_proc.GetString(i))
            if text is 'Split files':
                self.split = list('TRUE')
        except KeyError:
            pass
        self.send_settings()


    def load_settings(self, msg):
        """
        Listener function that changes input values
        to values specified in settings file.
        """
        try:
            if msg['cluster'] is not None:
                self.cluster = msg['cluster']
                clus = self.cluster_choice.FindString(msg['cluster'][0])
                self.cluster_choice.SetSelection(clus)
                if msg['split'] is 'TRUE':
                    self.cluster_proc.SetSelection(1)
                else:
                    self.cluster_proc.SetSelection(0)
        except Exception:
            logger.error("Unable to load settings. ", exc_info=True)

    def register_figures(self, event):
        """Registers listbox event and calls generate_figures."""
        self.generate_cluster_figures()

    def generate_cluster_figures(self):
        """Generates figures for diagnostics canvas."""
        from massoc.scripts.batch import Batch
        from sklearn.cluster import KMeans, DBSCAN, SpectralClustering, AffinityPropagation
        from sklearn.mixture import GaussianMixture
        from sklearn.metrics import silhouette_score
        from sklearn.decomposition import PCA
        nums = list(range(2, 5))
        try:
            file = self.file_list.GetSelection()
            file = self.file_list.GetString(file)
            x = 'init'
            biomfile = {x: biom.load_table(file)}
            algo = self.cluster_choice.GetSelection()
            algo = self.cluster_choice.GetString(algo)
            inputs = {'biom_file': [file],
                      'cluster': [algo]}
            normbatch = Batch(biomfile, inputs)
            normbatch = normbatch.normalize_transform(mode='clr')
            norm_table = normbatch.otu[x]
            topscore = 0
            bestcluster = [1] * len(norm_table.ids())
            data = csr_matrix.todense(norm_table.matrix_data)
            data = np.matrix.transpose(data)
            data = PCA(n_components=2).fit_transform(data)
            randomclust = np.random.randint(2, size=len(data))
            sh_score = [silhouette_score(data, randomclust)]
            # K-means clustering, tests 2-4 clusters
            if inputs['cluster'][0] == 'K-means':
                for i in nums:
                    clusters = KMeans(i).fit_predict(data)
                    silhouette_avg = silhouette_score(data, clusters)
                    sh_score.append(silhouette_avg)
                topscore = int(np.argmax(sh_score) + 1)
                bestcluster = KMeans(topscore).fit_predict(data)
            # DBSCAN clustering, automatically finds optimal cluster size
            if inputs['cluster'][0] == 'DBSCAN':
                bestcluster = DBSCAN().fit_predict(data)
                topscore = len(set(bestcluster)) - (1 if -1 in bestcluster else 0)
            # Gaussian Mixture Model (gmm) probability distribution
            if inputs['cluster'][0] == 'Gaussian':
                for i in nums:
                    fit = GaussianMixture(i).fit(data)
                    clusters = fit.predict(data)
                    silhouette_avg = silhouette_score(data, clusters)
                    sh_score.append(silhouette_avg)
                topscore = int(np.argmax(sh_score) + 1)
                bestfit = GaussianMixture(topscore).fit(data)
                bestcluster = bestfit.predict(data)
            # Spectral Clustering
            if inputs['cluster'][0] == 'Spectral':
                for i in nums:
                    clusters = SpectralClustering(i).fit_predict(data)
                    silhouette_avg = silhouette_score(data, clusters)
                    sh_score.append(silhouette_avg)
                topscore = int(np.argmax(sh_score) + 1)
                bestcluster = SpectralClustering(topscore).fit_predict(data)
            # Affinity Propagation clustering
            if inputs['cluster'] == 'Affinity':
                bestcluster = AffinityPropagation().fit_predict(data)
                topscore = len(set(bestcluster)) - (1 if -1 in bestcluster else 0)
            if max(sh_score) < 0.25:
                raise ValueError("Silhouette score too low: please try a different algorithm. "
                                 "Your data may not be suitable for clustering.")
            for i in range(topscore):
                mask, = np.where(bestcluster == i)
                for j in mask:
                    norm_table._sample_metadata[j]['cluster'] = inputs['cluster'][0] + '_' + str(i)
            x, y = zip(*data)
            self.prev.scatter(x, y, bestcluster)
            self.canvas1.draw()
        except Exception:
            logger.error("Failed to generate figures. ", exc_info=True)


