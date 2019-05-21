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

import wx
from wx.lib.pubsub import pub
import massoc
from run_massoc import run_network, get_input
from massoc.scripts.batch import read_settings
from time import sleep
import sys
import logging
from copy import deepcopy
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

# handler to file
# only handler with 'w' mode, rest is 'a'
# once this handler is started, the file writing is cleared
# other handlers append to the file
logpath = "\\".join(os.getcwd().split("\\")[:-1]) + '\\massoc.log'
# filelog path is one folder above massoc
# pyinstaller creates a temporary folder, so log would be deleted
fh = logging.handlers.RotatingFileHandler (maxBytes=500,
                                      filename=logpath, mode='a')
fh.setLevel(logging.INFO)
formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
fh.setFormatter(formatter)
logger.addHandler(fh)


class NetworkPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # subscribe to inputs from tabwindow
        self.frame = parent

        btnsize = (300, -1)
        boxsize = (300, 50)
        # adds columns

        self.tools = None
        self.spar_boot = None
        self.spar_pval = None
        self.cores = 4
        self.conet = None
        self.conet_bash = None
        self.spar = None

        # while mainframe contains the main settings file, the NetworkPanel also needs it to generate a command call
        self.settings = dict()
        self.settings['procbioms'] = None
        self.settings['network'] = None

        pub.subscribe(self.review_settings, 'show_settings')
        pub.subscribe(self.format_settings, 'input_settings')
        pub.subscribe(self.format_settings, 'process_settings')
        pub.subscribe(self.format_settings, 'network_settings')
        pub.subscribe(self.load_settings, 'load_settings')
        pub.subscribe(self.toggle_network, 'toggle_network')

        self.topsizer = wx.BoxSizer(wx.HORIZONTAL)

        # defines columns
        self.leftsizer = wx.BoxSizer(wx.VERTICAL)
        self.topleftsizer = wx.BoxSizer(wx.VERTICAL)
        self.bottomleftsizer = wx.BoxSizer(wx.VERTICAL)

        # include a column to report errors or check file formats
        self.rightsizer = wx.BoxSizer(wx.VERTICAL)

        # select tools
        self.tool_txt = wx.StaticText(self, label='Select network inference tools to run')
        self.tool_box = wx.ListBox(self, choices=['SparCC', 'CoNet', 'SPIEC-EASI'],
                                        size=(boxsize[0], 110), style=wx.LB_MULTIPLE)
        self.tool_box.Bind(wx.EVT_MOTION, self.update_help)
        self.tool_box.Bind(wx.EVT_LISTBOX, self.list_tools)

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
        self.conet_sets = wx.StaticText(self, label='To adjust CoNet settings, provide an alternative CoNet BASH call.')
        self.conet_bash_txt = wx.TextCtrl(self, value='', size=btnsize)
        self.conetbox.Add(self.conet_sets, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.conetbox.AddSpacer(5)
        self.conetbox.Add(self.conet_bash_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
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

        # Run button
        self.go = wx.Button(self, label='Run network inference', size=(btnsize[0], 40))
        self.go.Bind(wx.EVT_BUTTON, self.run_network)
        self.go.SetFont(wx.Font(16, wx.DECORATIVE, wx.NORMAL, wx.BOLD))
        self.go.SetBackgroundColour(wx.Colour(0, 153, 51))

        self.rightsizer.Add(self.jobs_text, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.jobs_cores, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.jobs_choice, flag=wx.ALIGN_LEFT)
        self.rightsizer.AddSpacer(40)
        self.rightsizer.Add(self.rev_text, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.review, flag=wx.ALIGN_LEFT)
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
            dlg = wx.FileDialog(
                self, message="Select count tables",
                defaultFile="",
                style=wx.FD_OPEN | wx.FD_MULTIPLE | wx.FD_CHANGE_DIR
            )
            if dlg.ShowModal() == wx.ID_OK:
                path = dlg.GetPath()
                path = path.replace('\\', '/')
                self.conet_bash = path
                if len(path) > 0:
                    self.conet_bash_txt.SetValue(self.conet_bash)
            self.send_settings()
            dlg.Destroy()
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
            self.conet = dlg.GetPath()
        self.conet_txt.SetValue(self.conet)
        self.send_settings()
        dlg.Destroy()

    def open_spar(self, event):
        dlg = wx.DirDialog(self, "Select SparCC directory", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.spar = dlg.GetPath()
        self.spar_txt.SetValue(self.spar)
        self.send_settings()
        dlg.Destroy()

    def list_tools(self, event):
        ids = self.tool_box.GetSelections()
        text = list()
        for i in ids:
            text.append(self.tool_box.GetString(i))
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
        self.cores = int(self.jobs_choice.GetValue())
        self.send_settings()

    def change_spiec_alg(self, event):
        """Takes settings input from GUI and uses it to rewrite the SPIEC-EASI script in execs."""
        try:
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
        except Exception:
            logger.error("Failed to change SPIEC-EASI algorithm, ", exc_info=True)

    def change_spiec_stars(self, event):
        """Takes settings input from GUI and uses it to rewrite the SPIEC-EASI script in execs."""
        try:
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
        except Exception:
            logger.error("Failed to change number of repetitions. ", exc_info=True)

    def run_network(self, event):
        self.settings['network'] = list()
        self.settings['procbioms'] = list()
        network_names = list()
        biom_names = list()
        for tool in self.settings['tools']:
            for level in self.settings['levels']:
                for name in self.settings['name']:
                    filename = self.settings['fp'] + '/' + tool + '_' + name + '_' + level + '.txt'
                    biomname = self.settings['fp'] + '/' + name + '_' + level + '.hdf5'
                    biom_names.append(biomname)
                    network_names.append(filename)
        self.settings['network'] = network_names
        self.settings['procbioms'] = biom_names
        try:
            eg = Thread(target=massoc_worker, args=(self.settings,))
            eg.start()
            dlg = LoadingBar(self.settings)
            dlg.ShowModal()
            eg.join()  # necessary for ubuntu thread to quit crashing
        except Exception:
            logger.error("Failed to start worker thread. ", exc_info=True)
        self.send_settings()


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
            if i in self.settings and self.settings[i] is not None:
                subcommand.append(command_dict[i])
                if type(self.settings[i]) == list:
                    multicommand = ' '.join(self.settings[i])
                    subcommand.append(multicommand)
                else:
                    subcommand.append(str(self.settings[i]))
                reviewtext.append(' '.join(subcommand))
        reviewtext = '\n'.join(reviewtext)
        self.review.SetValue(reviewtext)

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'tools': self.tools, 'spar_pval': self.spar_pval, 'spar_boot': self.spar_boot,
                    'cores': self.cores, 'spar': self.spar, 'conet': self.conet, 'conet_bash': self.conet_bash,
                    'spiec': None,  # it is not possible to specify R file from GUI
                    'procbioms': self.settings['procbioms'], 'network': self.settings['network']}
        pub.sendMessage('network_settings', msg=settings)

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
        if 'tools' in self.settings and 'levels' in self.settings and \
                self.settings['tools'] is not None and self.settings['levels'] is not None:
            njobs = len(self.settings['tools']) * len(self.settings['levels'])
            njobs = 'Number of jobs to run: ' + str(njobs)
            self.jobs_text.SetLabel(label=njobs)

    def load_settings(self, msg):
        """
        Listener function that changes input values
        to values specified in settings file.
        """
        self.settings = deepcopy(msg)
        if msg['tools'] is not None:
            self.tools = msg['tools']
            tooldict = {'sparcc': 0, 'conet': 1, 'spiec-easi': 2}
            for tool in msg['tools']:
                self.tool_box.SetSelection(tooldict[tool])
        else:
            self.tools = None
            choice = self.tool_box.GetSelections()
            for selection in choice:
                self.tool_box.Deselect(selection)
        if msg['spar_boot'] is not None:
            self.sparbox.ShowItems(show=True)
            self.Layout()
            self.spar_boot = msg['spar_boot']
            self.spar_boot_txt.SetValue(msg['spar_boot'])
        else:
            self.sparbox.ShowItems(show=False)
            self.spar_boot_txt.SetValue('')
            self.spar_boot = None
        if msg['spar_pval'] is not None:
            self.sparbox.ShowItems(show=True)
            self.Layout()
            self.spar_pval = msg['spar_pval']
            self.spar_pval_txt.SetValue(msg['spar_pval'])
        else:
            self.sparbox.ShowItems(show=False)
            self.spar_pval_txt.SetValue('')
            self.spar_pval = None
        if msg['cores'] is not None:
            self.cores = msg['cores']
            self.jobs_choice.SetValue(str(msg['cores']))
        else:
            self.cores = None
            self.jobs_choice.SetValue('')
        if msg['conet'] is not None:
            self.conet = msg['conet']
            self.conet_txt.SetValue(self.conet)
        if msg['conet_bash'] is not None:
            self.conet_bash = msg['conet_bash']
            self.conet_bash_txt.SetValue(msg['conet_bash'])
        else:
            self.conet_bash = None
            self.conet_bash_txt.SetValue('')
        if msg['spar'] is not None:
            self.spar = msg['spar']
            self.spar_txt.SetValue(msg['spar'])
        else:
            self.spar = None
            self.spar_txt.SetValue('')
        self.send_settings()

    def toggle_network(self, msg):
        if msg == 'Yes':
            self.tool_txt.Enable(True)
            self.jobs_choice.Enable(True)
            self.go.Enable(True)
            self.tool_box.Enable(True)
            self.conet_button.Enable(True)
            self.spar_button.Enable(True)
            self.settings_txt.Enable(True)
            self.settings_choice.Enable(True)
            self.jobs_cores.Enable(True)
            self.jobs_text.Enable(True)
        if msg == 'No':
            self.tool_txt.Enable(False)
            self.jobs_choice.Enable(False)
            self.go.Enable(False)
            self.tool_box.Enable(False)
            self.conet_button.Enable(False)
            self.spar_button.Enable(False)
            self.settings_txt.Enable(False)
            self.settings_choice.Enable(False)
            self.jobs_cores.Enable(False)
            self.jobs_text.Enable(False)


class LoadingBar(wx.Dialog):
    def __init__(self, settings):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Progress")
        if settings['tools'] is not None and settings['levels'] is not None:
            self.number_networks = (len(settings['tools']) * len(settings['levels'])) + 1
        self.count = 0
        self.text = wx.StaticText(self, label="Starting...")
        self.progress = wx.Gauge(self, range=self.number_networks+2)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.text, 0, wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.progress, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(sizer)
        pub.subscribe(self.get_progress, "update")

    def get_progress(self, msg):
        """Progress bar appears to work, but end /start is not calculated correctly.
        Issue is trivial but may be nice to fix. """
        if msg is 'Writing to disk...' or 'Starting network inference. This may take some time!':
            self.count += 1
        self.progress.SetValue(self.count)
        self.text.SetLabel(msg)
        if msg == 'Finished running network inference!':
            sleep(3)
            self.Destroy()


def massoc_worker(inputs):
    """
    Alternative version of massoc's main pipe.
    Uses publisher to send messages instead of sys.stdout.write.
    """
    get_input(inputs, publish=True)
    run_network(inputs, publish=True)



