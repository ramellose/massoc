"""
The network panel allows users to select network inference tools to run and to adjust the settings of these tools.
It also provides an overview of current settings, and can execute massoc with these settings.
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import os
from threading import Thread

import numpy as np
import wx
from biom import load_table
from biom.parse import MetadataMap
from massoc.scripts.batch import Batch
from massoc.scripts.netwrap import Nets
from wx.lib.pubsub import pub

import massoc
from massoc.scripts.main import run_parallel


class NetworkPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # subscribe to inputs from tabwindow
        self.frame = parent

        btnsize = (300, -1)
        btnmargin = 10
        boxsize = (300, 50)
        # adds columns

        self.tools = None
        self.spar_boot = None
        self.spar_pval = None
        self.cores = None
        self.conet = None
        self.spar = None

        # while mainframe contains the main settings file, the NetworkPanel also needs it to generate a command call
        self.settings = {'biom_file': None, 'otu_table': None, 'tax_table': None, 'sample_data': None,
                         'otu_meta': None, 'cluster': None, 'split': None, 'prev': None, 'fp': None,
                         'levels': None, 'tools': None, 'spiec': None, 'conet': None, 'spar': None, 'spar_pval': None,
                         'spar_boot': None, 'nclust': None, 'name': None, 'cores': None, 'rar': None, 'min': None}

        pub.subscribe(self.review_settings, 'show_settings')
        pub.subscribe(self.format_settings, 'update_settings')
        pub.subscribe(self.clear_settings_net, 'clear_settings')
        pub.subscribe(self.load_settings_net, 'load_settings')

        self.topsizer = wx.BoxSizer(wx.HORIZONTAL)

        # defines columns
        self.leftsizer = wx.BoxSizer(wx.VERTICAL)
        self.topleftsizer = wx.BoxSizer(wx.VERTICAL)
        self.bottomleftsizer = wx.BoxSizer(wx.VERTICAL)

        # include a column to report errors or check file formats
        self.rightsizer = wx.BoxSizer(wx.VERTICAL)

        # select tools
        self.tool_txt = wx.StaticText(self, label='Select network inference tools to run')
        self.tool_box = wx.CheckListBox(self, choices=['SparCC', 'CoNet', 'SPIEC-EASI'],
                                        size=(boxsize[0], 110), style=wx.LB_MULTIPLE)
        self.tool_box.Bind(wx.EVT_MOTION, self.update_help)
        self.tool_box.Bind(wx.EVT_CHECKLISTBOX, self.list_tools)

        # select CoNet.jar
        self.conet_button = wx.Button(self, label="Select CoNet3 folder", size=btnsize)
        self.conet_button.Bind(wx.EVT_BUTTON, self.open_conet)
        self.conet_button.Bind(wx.EVT_MOTION, self.update_help)
        self.conet_txt = wx.TextCtrl(self, size=btnsize)

        # select SparCC folder
        self.spar_button = wx.Button(self, label="Select SparCC folder", size=btnsize)
        self.spar_button.Bind(wx.EVT_BUTTON, self.open_spar)
        self.spar_button.Bind(wx.EVT_MOTION, self.update_help)
        self.spar_txt = wx.TextCtrl(self, size=btnsize)

        # show & adjust tool settings
        self.settings_txt = wx.StaticText(self, label="Change settings for: ")
        self.settings_choice = wx.ListBox(self, choices=['SparCC', 'CoNet', 'SPIEC-EASI'], size=boxsize)
        self.settings_choice.Bind(wx.EVT_LISTBOX, self.show_settings)

        # SparCC settings
        self.sparbox = wx.BoxSizer(wx.VERTICAL)
        self.spar_txt1 = wx.StaticText(self, label='Number of SparCC bootstrap iterations')
        self.spar_boot_txt = wx.TextCtrl(self, value='100', size=btnsize)
        self.spar_boot_txt.Bind(wx.EVT_TEXT, self.send_spar)
        self.spar_txt2 = wx.StaticText(self, label='Threshold for pseudo p-values')
        self.spar_pval_txt = wx.TextCtrl(self, value='0.001', size=btnsize)
        self.spar_pval_txt.Bind(wx.EVT_TEXT, self.send_spar)
        self.sparbox.Add(self.spar_txt1, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.sparbox.Add(self.spar_boot_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.sparbox.Add(self.spar_txt2, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.sparbox.Add(self.spar_pval_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.sparbox.ShowItems(show=False)

        # CoNet settings
        self.conetbox = wx.BoxSizer(wx.VERTICAL)
        self.conet_sets = wx.StaticText(self, label='To adjust CoNet settings, replace the CoNet BASH call \n'
                                                   'in this file by new ones generated with the original \n'
                                                   'CoNet application. ')
        self.conet_file = wx.GenericDirCtrl(self, dir=(os.path.dirname(massoc.__file__) + '\\execs\\CoNet.sh'),
                                            size=(boxsize[0], 300))
        self.conetbox.Add(self.conet_sets, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.conetbox.AddSpacer(5)
        self.conetbox.Add(self.conet_file, flag=wx.ALIGN_LEFT)
        self.conetbox.ShowItems(show=False)

        # SPIEC-EASI settings
        self.spiecbox = wx.BoxSizer(wx.VERTICAL)
        self.spiec_txt1 = wx.StaticText(self, label='SPIEC-EASI algorithm')
        self.spiec_alg = wx.ListBox(self, choices=['Meinshausen-Buhlmann', 'Graphical Lasso'], size=(boxsize[0], 40))
        self.spiec_alg.Bind(wx.EVT_LISTBOX, self.change_spiec_alg)
        self.spiec_txt2 = wx.StaticText(self, label='Number of StARS repetitions')
        self.spiec_star = wx.TextCtrl(self, value='50', size=btnsize)
        self.spiec_star.Bind(wx.EVT_TEXT, self.change_spiec_stars)
        self.spiecbox.Add(self.spiec_txt1, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.spiecbox.Add(self.spiec_alg, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.spiecbox.Add(self.spiec_txt2, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.spiecbox.Add(self.spiec_star, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.spiecbox.ShowItems(show=False)

        # multiprocessing dialog
        self.jobs_text = wx.StaticText(self, label='Number of jobs to run: ')
        self.jobs_cores = wx.StaticText(self, label='Number of processes for network inference:')
        self.jobs_choice = wx.TextCtrl(self, value='4', size=btnsize)
        self.jobs_choice.Bind(wx.EVT_TEXT, self.send_cores)
        self.jobs_choice.Bind(wx.EVT_MOTION, self.update_help)

        # review settings
        self.rev_text = wx.StaticText(self, label='Current settings')
        self.review = wx.TextCtrl(self, value='', size=(boxsize[0], 200), style=wx.TE_READONLY | wx.TE_MULTILINE)
        self.call = wx.Button(self, label='Export as command line call', size=btnsize)
        self.call.Bind(wx.EVT_BUTTON, self.generate_call)
        self.call.Bind(wx.EVT_MOTION, self.update_help)
        self.call.SetFont(wx.Font(16, wx.DECORATIVE, wx.NORMAL, wx.BOLD))

        # Run button
        self.go = wx.Button(self, label='Run network inference', size=(btnsize[0], 40))
        self.go.Bind(wx.EVT_BUTTON, self.run_network)
        self.go.SetFont(wx.Font(18, wx.DECORATIVE, wx.NORMAL, wx.BOLD))
        self.go.SetBackgroundColour(wx.Colour(0, 153, 51))

        self.rightsizer.Add(self.jobs_text, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.jobs_cores, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.jobs_choice, flag=wx.ALIGN_LEFT)
        self.rightsizer.Add(self.rev_text, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.review, flag=wx.ALIGN_LEFT)
        self.rightsizer.Add(self.call, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.go, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.topleftsizer.Add(self.tool_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.AddSpacer(5)
        self.topleftsizer.Add(self.tool_box, flag=wx.ALIGN_LEFT)
        self.topleftsizer.AddSpacer(40)
        self.topleftsizer.Add(self.conet_button, flag=wx.ALIGN_LEFT)
        self.topleftsizer.Add(self.conet_txt, flag=wx.ALIGN_LEFT)
        self.topleftsizer.AddSpacer(5)
        self.topleftsizer.Add(self.spar_button, flag=wx.ALIGN_LEFT)
        self.topleftsizer.Add(self.spar_txt, flag=wx.ALIGN_LEFT)
        self.topleftsizer.AddSpacer(40)
        self.topleftsizer.Add(self.settings_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.topleftsizer.Add(self.settings_choice, flag=wx.ALIGN_CENTER_HORIZONTAL)

        self.bottomleftsizer.Add(self.sparbox, wx.ALL, 20)
        self.bottomleftsizer.Add(self.conetbox, wx.ALL, 20)
        self.bottomleftsizer.Add(self.spiecbox, wx.ALL, 20)

        self.leftsizer.Add(self.topleftsizer, wx.ALL, 20)
        self.leftsizer.AddSpacer(10)
        self.leftsizer.Add(self.bottomleftsizer, wx.ALL, 20)
        self.topsizer.Add(self.leftsizer, 0, wx.ALL, 20)
        self.topsizer.Add(self.rightsizer, 0, wx.ALL, 20)
        self.SetSizerAndFit(self.topsizer)

        # help strings for buttons
        self.buttons = {self.tool_box: 'Run SparCC, CoNet or SPIEC-EASI. Check the help files for more information.',
                        self.call: 'Generate a command line call to run this pipeline.',
                        self.jobs_choice: 'Distributing jobs over multiple processes = lower runtime.',
                        self.conet_button: 'Select the location of your CoNet.jar file. ',
                        self.spar_button: 'Select the location of your SparCC folder. '}

    def update_help(self, event):
        btn = event.GetEventObject()
        if btn in self.buttons:
            status = self.buttons[btn]
            pub.sendMessage('change_statusbar', msg=status)

    def show_settings(self, event):
        tool = self.settings_choice.GetStringSelection()
        if tool == 'SparCC':
            self.sparbox.ShowItems(show=True)
            self.conetbox.ShowItems(show=False)
            self.spiecbox.ShowItems(show=False)
        if tool == 'CoNet':
            self.sparbox.ShowItems(show=False)
            self.conetbox.ShowItems(show=True)
            self.spiecbox.ShowItems(show=False)
        if tool == 'SPIEC-EASI':
            self.sparbox.ShowItems(show=False)
            self.conetbox.ShowItems(show=False)
            self.spiecbox.ShowItems(show=True)
        self.Layout()

    def open_conet(self, event):
        """
        Create file dialog and show it.
        """
        dlg = wx.DirDialog(self, "Select CoNet3 directory", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.conet = list()
            self.conet.append(dlg.GetPath())
        self.conet_txt.SetValue(self.conet[0])
        self.send_settings()
        dlg.Destroy()

    def open_spar(self, event):
        dlg = wx.DirDialog(self, "Select SparCC directory", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.spar = list()
            self.spar.append(dlg.GetPath())
        self.spar_txt.SetValue(self.spar[0])
        self.send_settings()
        dlg.Destroy()

    def list_tools(self, event):
        text = list(self.tool_box.GetCheckedStrings())
        self.tools = [x.lower() for x in text]
        self.send_settings()

    def send_spar(self, event):
        btn = event.GetEventObject()
        if btn is self.spar_pval_txt:
            self.spar_pval = list()
            self.spar_pval.append(self.spar_pval_txt.GetValue())
        if btn is self.spar_boot_txt:
            self.spar_boot = list()
            self.spar_boot.append(self.spar_boot_txt.GetValue())
        self.send_settings()

    def send_cores(self, event):
        self.cores = list()
        self.cores.append(self.jobs_choice.GetValue())
        self.send_settings()

    def change_spiec_alg(self, event):
        """Takes settings input from GUI and uses it to rewrite the SPIEC-EASI script in execs."""
        path = os.path.dirname(massoc.__file__) + '\\execs\\spieceasi.r'
        path = path.replace('\\', '/')
        spiec_script = open(path, 'r')
        commands = spiec_script.read()
        spiec_script.close()
        alg = self.spiec_alg.GetSelection()
        text = self.spiec_alg.GetString(alg)
        if text == 'Meinshausen-Buhlmann':
            commands = commands.replace('method = "glasso"', 'method = "mb"')
        if text == 'Graphical Lasso':
            commands = commands.replace('method = "mb"', 'method = "glasso"')
        spiec_script = open(path, 'w')
        spiec_script.write(commands)
        spiec_script.close()

    def change_spiec_stars(self, event):
        """Takes settings input from GUI and uses it to rewrite the SPIEC-EASI script in execs."""
        path = os.path.dirname(massoc.__file__) + '\\execs\\spieceasi.r'
        path = path.replace('\\', '/')
        spiec_script = open(path, 'r')
        commands = spiec_script.read()
        spiec_script.close()
        new_stars = 'rep.num=' + self.spiec_star.GetValue() + '))'
        old_stars = commands.find('rep.num=')
        old_stars = commands[old_stars:(old_stars+20)]  # need to make sure there is a \n in there
        old_stars = old_stars.split('\n')[0]
        commands = commands.replace(old_stars, new_stars)
        spiec_script = open(path, 'w')
        spiec_script.write(commands)
        spiec_script.close()

    def generate_call(self, event):
        # otu_meta, nclust are not integrated in GUI yet
        command = list()
        command.append('main.py')
        command_dict = {'biom_file': '-biom', 'otu_table': '-otu', 'tax_table': '-tax', 'sample_data': '-s',
                        'otu_meta': '-od', 'cluster': '-cl', 'split': '-split', 'nclust': '-nclust', 'cores': '-cores',
                        'prev': '-prev', 'fp': '-o', 'levels': '-levels', 'tools': '-tools', 'spiec': '-spiec',
                        'conet': '-conet', 'spar': 'spar', 'spar_pval': '-spar_pval', 'spar_boot': '-spar_boot',
                        'name': '-n', 'rar': '-rar', 'min': '-min'}
        for i in command_dict:
            if self.settings[i] is not None:
                if len(self.settings[i]) > 0:
                    command.append(command_dict[i])
                    if len(self.settings[i]) > 1:
                        subcommand = ' '.join(self.settings[i])
                        command.append(subcommand)
                    else:
                        command.append(self.settings[i][0])
        command = ' '.join(command)
        dlg = wx.TextEntryDialog(None, "", "Command line call", command)
        dlg.ShowModal()

    def run_network(self, event):
        eg = Thread(target=massoc_worker, args=(self.settings,))
        eg.start()
        dlg = LoadingBar(self.settings)
        dlg.ShowModal()

    def review_settings(self, msg):
        """
        Listener function for settings, adjusts review pane.
        """
        reviewtext = list()
        command_dict = {'biom_file': 'BIOM file:', 'otu_table': 'OTU table:', 'tax_table': 'Taxonomy table:',
                        'sample_data': 'Sample data:',
                        'otu_meta': 'Taxa metadata:', 'cluster': 'Cluster algorithm:', 'split': 'Split files by:',
                        'nclust': 'Maximum number of clusters to evaluate:',
                        'prev': 'Prevalence filter:', 'levels': 'Taxonomic levels:',
                        'tools': 'Network inference tools to run:', 'spiec': 'Settings for SPIEC-EASI:',
                        'conet': 'Location of CoNet executable:', 'spar': 'Location of SparCC folder',
                        'spar_pval': 'SparCC p-value:',
                        'spar_boot': 'SparCC bootstrap number:', 'name': 'Prefix for output files:',
                        'fp': 'Location for output files:', 'cores': 'Number of processes:', 'rar': 'Rarefaction:',
                        'min': 'Minimal mean count:'}
        for i in command_dict:
            subcommand = list()
            if self.settings[i] is not None:
                if len(self.settings[i]) > 0:
                    subcommand.append(command_dict[i])
                    if len(self.settings[i]) > 1:
                        multicommand = ' '.join(self.settings[i])
                        subcommand.append(multicommand)
                    else:
                        subcommand.append(self.settings[i][0])
                    reviewtext.append(' '.join(subcommand))
        reviewtext = '\n'.join(reviewtext)
        self.review.SetValue(reviewtext)

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'tools': self.tools, 'spar_pval': self.spar_pval, 'spar_boot': self.spar_boot,
                    'cores': self.cores, 'spar': self.spar, 'conet': self.conet}
        pub.sendMessage('update_settings', msg=settings)

    def format_settings(self, msg):
        """
        Listener function for settings from tabs in notebook.
        """
        try:
            for key in msg:
                self.settings[key] = msg[key]
        except KeyError:
            pass
        if key in msg is 'cluster':  # at the moment, the maximum # of clusters to evaluate cannot be selected
            self.settings['nclust'] = [4]
        if self.settings['tools'] is not None and self.settings['levels'] is not None:
            njobs = len(self.settings['tools']) * len(self.settings['levels'])
            njobs = 'Number of jobs to run: ' + str(njobs)
            self.jobs_text.SetLabel(label=njobs)

    def clear_settings_net(self, msg):
        """
        Listener function that clears all input boxes.
        The new "empty" settings are not sent to
        the mainframe or network tabs,
        because those also receive this message.
        """
        if msg == 'CLEAR':
            self.jobs_choice.SetValue("")
            self.review.SetValue("")
            for cb in self.tool_box.Checked:
                self.tool_box.Check(cb, False)
            for cb in np.nditer(self.settings_choice.GetSelections):
                self.settings_choice.Deselect(cb)
            self.spar_boot.SetValue("100")
            self.spar_pval.SetValue("0.001")
            for cb in np.nditer(self.spiec_alg.GetSelections):
                self.spiec_alg.Deselect(cb)

            self.spiec_star.SetValue("50")
            self.settings = {'biom_file': None, 'otu_table': None, 'tax_table': None, 'sample_data': None,
                             'otu_meta': None, 'cluster': None, 'split': None, 'prev': None, 'fp': None,
                             'levels': None, 'tools': None, 'spiec': None, 'conet': None, 'spar': None,
                             'spar_pval': None, 'min': None,
                             'spar_boot': None, 'nclust': None, 'name': None, 'cores': None, 'rar': None}
            self.send_settings()

    def load_settings_net(self, msg):
        """
        Listener function that changes input values
        to values specified in settings file.
        """
        self.settings = msg
        if msg['tools'] is not None:
            self.tools = msg['tools']
            tooldict = {'sparcc': 0, 'conet': 1, 'spiec-easi': 2}
            for tool in msg['tools']:
                self.tool_box.Check(tooldict[tool], True)
        if msg['spar_boot'] is not None:
            self.sparbox.ShowItems(show=True)
            self.Layout()
            self.spar_boot = msg['spar_boot']
            self.spar_boot_txt.SetValue(msg['spar_boot'][0])
        if msg['spar_pval'] is not None:
            self.sparbox.ShowItems(show=True)
            self.Layout()
            self.spar_pval = msg['spar_pval']
            self.spar_pval_txt.SetValue(msg['spar_pval'][0])
        if msg['cores'] is not None:
            self.cores = msg['cores']
            self.jobs_choice.SetValue(msg['cores'][0])
        self.send_settings()


class LoadingBar(wx.Dialog):
    def __init__(self, settings):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Progress")
        if settings['tools'] is not None and settings['levels'] is not None:
            self.number_networks = (len(settings['tools']) * len(settings['levels'])) + 1

        self.count = 0
        self.text = wx.StaticText(self, label="Starting...")
        self.progress = wx.Gauge(self, range=self.number_networks)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.text, 0, wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.progress, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(sizer)
        pub.subscribe(self.get_progress, "update")

    def get_progress(self, msg):
        """This is not working yet! Find a better way to track progress of the thread... """
        if msg is 'Writing to disk...' or 'Starting network inference. This may take some time!':
            self.count += 1
        if msg is 'Starting network inference. This may take some time!':
            self.timer.Start()
        if msg is 'Finished running network inference!':
            self.timer.Stop()
            print(self.timer.GetInterval())
            self.Destroy()
        self.progress.SetValue(self.count)
        print(msg)
        self.text.SetLabel(msg)


def massoc_worker(inputs):
    """
    Alternative version of massoc's main pipe.
    Uses publisher to send messages instead of sys.stdout.write.
    """
    filestore = {}
    if inputs['biom_file'] is None:
        if inputs['otu_table'] is None:
            raise ValueError("Please supply either a biom file "
                             "or a tab-delimited OTU table!")
    i = 0
    if inputs['biom_file'] is not None:
        for x in inputs['biom_file']:
            biomtab = load_table(x)
            filestore[inputs['name'][i]] = biomtab
            i += 1
    if inputs['otu_table'] is not None:
        j = 0  # j is used to match sample + tax data to OTU data
        for x in inputs['otu_table']:
            input_fp = x
            sample_metadata_fp = None
            observation_metadata_fp = None
            obs_data = None
            sample_data = None
            biomtab = load_table(input_fp)
            try:
                sample_metadata_fp = inputs['sample_data'][j]
                observation_metadata_fp = inputs['tax_table'][j]
            except KeyError:
                pass
            if sample_metadata_fp is not None:
                sample_f = open(sample_metadata_fp, 'r')
                sample_data = MetadataMap.from_file(sample_f)
                sample_f.close()
                biomtab.add_metadata(sample_data, axis='sample')
            if observation_metadata_fp is not None:
                obs_f = open(observation_metadata_fp, 'r')
                obs_data = MetadataMap.from_file(obs_f)
                obs_f.close()
                # for taxonomy collapsing,
                # metadata variable needs to be a complete list
                # not separate entries for each tax level
                for i in list(obs_data):
                    tax = list()
                    for j in list(obs_data[i]):
                        tax.append(obs_data[i][j])
                        obs_data[i].pop(j, None)
                    obs_data[i]['taxonomy'] = tax
                biomtab.add_metadata(obs_data, axis='observation')
            filestore[inputs['name'][i]] = biomtab
            i += 1
            j += 1
    bioms = Batch(filestore, inputs)
    if inputs['cluster'] is not None:
        pub.sendMessage('update', msg='Clustering BIOM files...')
        bioms.cluster_biom()
    if inputs['split'] is not None and inputs['split'] is not 'TRUE':
        bioms.split_biom()
    if inputs['min'] is not None:
        pub.sendMessage('update', msg='Setting minimum mean abundance...')
        bioms.prev_filter(mode='min')
    if inputs['rar'] is not None:
        pub.sendMessage('update', msg='Rarefying counts...')
        bioms.rarefy()
    if inputs['prev'] is not None:
        pub.sendMessage('update', msg='Setting prevalence filter...')
        bioms.prev_filter(mode='prev')
    nets = Nets(bioms)
    pub.sendMessage('update', msg='Starting network inference. This may take some time!')
    nets = run_parallel(nets)
    nets.write_networks()
    pub.sendMessage('update', msg="Finished running network inference!")

