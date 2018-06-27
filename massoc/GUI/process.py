"""
The process panel allows users to format files. They can split data by cluster ID or sample features,
or set a prevalence filter.
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'BSD'

import wx
from wx.lib.pubsub import pub
import numpy as np
import biom
from scipy.sparse import csr_matrix
import matplotlib
matplotlib.use('WXAgg')

from matplotlib.backends.backend_wxagg import FigureCanvasWxAgg as FigureCanvas
from matplotlib.backends.backend_wx import NavigationToolbar2Wx
from matplotlib.figure import Figure


class ProcessPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.frame = parent
        # subscribe to input tab
        pub.subscribe(self.receive_meta, 'receive_metadata')
        pub.subscribe(self.get_directory, 'update_settings')
        pub.subscribe(self.clear_settings_proc, 'clear_settings')
        pub.subscribe(self.load_settings_proc, 'load_settings')
        pub.subscribe(self.enable_tax, 'receive_tax')

        btnsize = (300, -1)
        btnmargin = 10
        boxsize = (300, 50)
        # adds columns

        self.settings = {'biom_file': None, 'otu_table': None, 'tax_table': None, 'sample_data': None,
                         'otu_meta': None, 'cluster': None, 'split': None, 'prev': None, 'fp': None,
                         'levels': None, 'tools': None, 'spiec': None, 'conet': None, 'spar_pval': None,
                         'spar_boot': None, 'nclust': None, 'name': None, 'cores': None, 'rar': None}

        self.name = None
        self.dir = None
        self.prev = ['20']
        self.split = None
        self.cluster = None
        self.agglom = None
        self.nclust = None
        self.rar = None
        self.min = None
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

        self.dg_button = wx.Button(self, label='Show data properties', size=btnsize)
        self.dg_button.Bind(wx.EVT_MOTION, self.update_help)
        self.dg_button.Bind(wx.EVT_BUTTON, self.show_diagnostics)

        # set minimum count
        self.min_txt = wx.StaticText(self, label='Remove taxa with mean count below:')
        self.min_number = wx.TextCtrl(self, value='0', size=btnsize)
        self.min_number.Bind(wx.EVT_TEXT, self.update_min)
        self.min_number.Bind(wx.EVT_MOTION, self.update_help)

        # set rarefaction
        self.rar_txt = wx.StaticText(self, label='Rarefy:')
        self.rar_choice = wx.ListBox(self, choices=['To even depth', 'To count number'], size=(300, 40))
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
        self.alg_txt = wx.StaticText(self, label='Clustering algorithm')
        self.cluster_choice = wx.ListBox(self, choices=['K-means', 'DBSCAN', 'Gaussian', 'Spectral',
                                                        'Affinity Propagation', 'None'], size=(boxsize[0], 100))
        self.cluster_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.cluster_choice.Bind(wx.EVT_LISTBOX, self.cluster_alg)
        self.act_txt = wx.StaticText(self, label='Use cluster ID to: ')
        self.cluster_proc = wx.ListBox(self, choices=['Add to metadata', 'Split files'], size=(boxsize[0], 40))
        self.cluster_proc.Bind(wx.EVT_LISTBOX, self.cluster_func)

        # select taxonomic levels
        self.tax_txt = wx.StaticText(self, label='Taxonomic levels to analyze')
        self.tax_choice = wx.CheckListBox(self, choices=['OTU', 'Species', 'Genus',
                                                         'Family', 'Order', 'Class', 'Phylum'],
                                          size=(boxsize[0], 130), style=wx.LB_MULTIPLE)
        self.tax_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.tax_choice.Bind(wx.EVT_CHECKLISTBOX, self.get_levels)
        self.tax_choice.Enable(False)

        self.topleftsizer.Add(self.prefix_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.prefix, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(40)
        self.topleftsizer.Add(self.dg_button, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(30)
        self.topleftsizer.Add(self.min_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.min_number, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(30)
        self.topleftsizer.Add(self.rar_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.rar_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.rar_number, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(30)
        self.topleftsizer.Add(self.prev_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.prev_val, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.prev_slider, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(40)
        self.topleftsizer.Add(self.split_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.split_list, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.rightsizer.Add(self.tax_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.tax_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(40)
        self.rightsizer.Add(self.cluster_btn, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.bottomleftsizer.Add(self.alg_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.bottomleftsizer.Add(self.cluster_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.bottomleftsizer.Add(self.act_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.bottomleftsizer.Add(self.cluster_proc, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.leftsizer.Add(self.topleftsizer, wx.ALL, 20)
        self.leftsizer.AddSpacer(10)
        self.rightsizer.Add(self.bottomleftsizer, wx.ALL, 20)
        self.bottomleftsizer.ShowItems(show=False)
        self.topsizer.Add(self.leftsizer, 0, wx.ALL, 20)
        self.topsizer.Add(self.rightsizer, 0, wx.ALL, 20)
        self.SetSizerAndFit(self.topsizer)

        # help strings for buttons
        self.buttons = {self.prev_val: 'With a prevalence filter of 20%, only taxa present in 20% of '
                                       'samples are retained. ',
                        self.prefix: 'Outputs generated by massoc follow a standard format. '
                                     'This prefix is used to generate filenames. ',
                        self.prev_slider: 'With a prevalence filter of 0.2, only taxa present in 20% '
                                       'of samples are retained. ',
                        self.split_list: 'Separate a file according to sample metadata and run network analysis on '
                                         'both of the split files. ',
                        self.cluster_btn: 'Cluster samples according to different algorithms; either use cluster '
                                          'identities to split files, or include them as additional metadata. ',
                        self.cluster_choice: 'For a detailed explanation of the clustering algorithms, please visit '
                                             'http://scikit-learn.org/stable/modules/'
                                             'classes.html#module-sklearn.cluster',
                        self.rar_choice: 'Specify whether samples should be rarefied to equal depth or a '
                                         'specific value. Samples below a specified threshold are removed.',
                        self.min_number: 'Taxa with low mean counts can be removed if there are too many '
                                         'high-prevalence, '
                                         'low-abundance taxa.',
                        self.dg_button: 'Show taxon prevalence and counts per sample for imported count tables.'}

    def update_help(self, event):
        btn = event.GetEventObject()
        if btn in self.buttons:
            status = self.buttons[btn]
            pub.sendMessage('change_statusbar', msg=status)

    def get_filename(self, event):
        text = self.prefix.GetValue()
        if self.name is None:
            self.name = list()
            self.name.append(text)
        else:
            self.name[0] = text
        self.send_settings()

    def get_rarefaction(self, event):
        rar = self.rar_choice.GetSelection()
        text = self.rar_choice.GetString(rar)
        self.rar = list()
        if text == 'To even depth':
            self.rar.append('True')
        if text == 'To count number':
            self.rar_number.Show()
            self.Layout()
        self.send_settings()

    def update_rarefaction(self, event):
        value = self.rar_number.GetValue()
        self.rar = list()
        try:
            self.rar.append(value)
            self.send_settings()
        except KeyError:
            pass

    def update_min(self, event):
        value = self.min_number.GetValue()
        self.min = list()
        try:
            self.min.append(value)
            self.send_settings()
        except KeyError:
            pass

    def enable_tax(self, msg):
        self.tax_choice.Enable(True)

    def split_files(self, event):
        self.split = list()
        name = self.split_list.GetSelection()
        text = self.split_list.GetString(name)
        try:
            self.split.append(text)
            self.send_settings()
        except KeyError:
            pass

    def cluster_alg(self, event):
        try:
            clus = self.cluster_choice.GetSelection()
            text = self.cluster_choice.GetString(clus)
            self.cluster = list()
            self.cluster.append(text)
        except KeyError:
            pass
        self.send_settings()

    def cluster_func(self, event):
        """
        This function specifies whether the cluster ID should be used to split files.
        """
        try:
            proc = self.cluster_proc.GetSelections()
            text = self.cluster_proc.GetString(proc)
            if text is 'Split files':
                if self.split is None:
                    self.split = list('TRUE')
                else:
                    self.split[0] = 'TRUE'
        except KeyError:
            pass
        self.send_settings()

    def show_diagnostics(self, event):
        """
        Starts diagnostics window to show prevalence + counts for imported tables.
        """
        dlg = Diagnostics(self.settings)
        dlg.ShowModal()

    def get_levels(self, event):
        text = list()
        try:
            text = self.tax_choice.GetCheckedStrings()
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
            self.prev = [str(value)]
            self.send_settings()
        except KeyError:
            pass

    def update_text(self, event):
        value = self.prev_slider.GetValue()
        try:
            value = str(value)
            self.prev_val.SetValue(value)
            self.prev = [str(value)]
        except KeyError:
            pass

    def show_dialog(self, event):
        self.bottomleftsizer.ShowItems(show=True)
        self.Layout()

    def get_directory(self, msg):
        try:
            self.dir = msg['dir']
        except KeyError:
            pass
        try:
            for key in msg:
                self.settings[key] = msg[key]
        except:
            pass

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'prev': self.prev, 'split': self.split, 'name': self.name,
                    'cluster': self.cluster, 'levels': self.agglom, 'nclust': self.nclust,
                    'rar': self.rar, 'min': self.min}
        pub.sendMessage('update_settings', msg=settings)

    def receive_meta(self, msg):
        """
        Listener function that registers whether a BIOM file with metadata
        or separate metadata file has been supplied.
        It opens the file to check whether there are variables present,
        and updates the split_list listbox to include these variables.
        """
        self.split_list.Set(msg)

    def clear_settings_proc(self, msg):
        """
        Listener function that clears all input boxes.
        The new "empty" settings are not sent to
        the mainframe or network tabs,
        because those also receive this message.
        """
        if msg == 'CLEAR':
            self.prefix.SetValue("")
            self.min_number.SetValue("0")
            self.rar_number.SetValue("0")
            for cb in np.nditer(self.rar_choice.GetSelections):
                self.rar_choice.Deselect(cb)
            self.prev_slider.SetTick(20)
            self.prev_val.SetValue("")
            for cb in np.nditer(self.split_list.GetSelections):
                self.split_list.Deselect(cb)
            for cb in np.nditer(self.cluster_choice.GetSelections):
                self.cluster_choice.Deselect(cb)
            for cb in np.nditer(self.cluster_proc.GetSelections):
                self.cluster_proc.Deselect(cb)
            for cb in self.tax_choice.Checked:
                self.tax_choice.Check(cb, False)

    def load_settings_proc(self, msg):
        """
        Listener function that changes input values
        to values specified in settings file.
        """
        if msg['min'] is not None:
            self.min = msg['min']
            self.min_number.SetValue(msg['min'][0])
        if msg['fp'] is not None:
            self.dir = msg['fp']
        if msg['name'] is not None:
            self.prefix.SetValue(msg['name'][0])
            self.name = msg['name']
        if msg['rar'] is not None:
            self.rar = msg['rar']
            rar = msg['rar'][0]
            try:
                rar = int(rar)
                self.rar_number.SetValue(rar)
            except ValueError:
                self.rar_choice.SetSelection(0)
        if msg['prev'] is not None:
            self.prev = msg['prev']
            self.prev_slider.SetTick(int(float((msg['prev'][0]))))
            self.prev_val.SetValue(msg['prev'][0])
        if msg['split'] is not None:
            self.split = msg['split']
            if msg['split'][0] is not 'TRUE':
                split = self.split_list.FindString(msg['split'][0])
                self.split_list.SetSelection(split)
        if msg['cluster'] is not None:
            self.bottomleftsizer.ShowItems(show=True)
            self.Layout()
            self.cluster = msg['cluster']
            clus = self.cluster_choice.FindString(msg['cluster'][0])
            self.cluster_choice.SetSelection(clus)
            if msg['split'] is 'TRUE':
                self.cluster_proc.SetSelection(1)
            else:
                self.cluster_proc.SetSelection(0)
        if msg['levels'] is not None:
            self.agglom = msg['levels']
            agglomdict = {'otu': 0, 'species': 1, 'genus': 2,
                          'family': 3, 'order': 4, 'class': 5,
                          'phylum': 6}
            for tax in msg['levels']:
                self.tax_choice.Check(agglomdict[tax], True)

class Diagnostics(wx.Dialog):
    def __init__(self, settings):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Data properties")
        filelist = list()
        if settings['biom_file'] is not None:
            for name in settings['biom_file']:
                filelist.append(name)
        if settings['otu_table'] is not None:
            for name in settings['otu_table']:
                filelist.append(name)
        self.file_list = wx.ListBox(self, choices=filelist)
        self.file_list.Bind(wx.EVT_LISTBOX, self.register_figures)
        self.file_list.SetSelection(0)
        self.figure1 = Figure(figsize=(3,2))
        self.prev = self.figure1.add_subplot(211)
        self.prev.set_xlabel('Prevalence')
        self.prev.set_title('Histogram of taxon prevalence')
        self.prev.set_ylabel('Number of taxa')
        self.figure1.tight_layout()
        self.figure2 = Figure(figsize=(3,2))
        self.rar = self.figure2.add_subplot(211)
        self.rar.set_xlabel('Count number')
        self.rar.set_title('Histogram of sample counts')
        self.rar.set_ylabel('Number of samples')
        self.figure2.tight_layout()
        self.canvas1 = FigureCanvas(self, -1, self.figure1)
        self.canvas2 = FigureCanvas(self, -1, self.figure2)

        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.file_list, 0, wx.EXPAND | wx.ALL, 10)
        sizer.AddSpacer(10)
        sizer.Add(self.canvas1, 0, wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.canvas2, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(sizer)
        self.Fit()
        self.generate_figures()

    def register_figures(self, event):
        """Registers listbox event and calls generate_figures."""
        self.generate_figures()

    def generate_figures(self):
        """Generates figures for diagnostics canvas."""
        file = self.file_list.GetSelection()
        file = self.file_list.GetString(file)
        biomfile = biom.load_table(file)
        data = biomfile.matrix_data
        data = csr_matrix.todense(data)
        fracs = np.count_nonzero(data, axis=1)
        nsamples = data.shape[1]
        fracs = fracs / nsamples
        self.prev.hist(fracs, bins=20)
        sample_sums = np.transpose(np.count_nonzero(data, axis=0))
        self.rar.hist(sample_sums, bins=40)
        self.canvas1.draw()
        self.canvas2.draw()





