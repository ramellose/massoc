"""
The network panel allows users to select network inference tools to run and to adjust the settings of these tools.
It also provides an overview of current settings, and can execute massoc with these settings.
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

from threading import Thread
import wx
from wx.lib.pubsub import pub
from massoc.scripts.main import general_settings
from massoc.scripts.netbase import ImportDriver
from massoc.scripts.netstats import Driver
import webbrowser
from biom import load_table
import networkx as nx
from networkx import NetworkXError
from subprocess import Popen
from massoc.scripts.main import resource_path
from psutil import Process
from time import sleep
from platform import system
import logging
import logging.handlers as handlers
logger = logging.getLogger()
hdlr = logging.FileHandler(resource_path("massoc.log"))
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)

class DataPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        # subscribe to inputs from tabwindow
        self.frame = parent

        btnsize = (300, -1)
        boxsize = (300, 50)
        # adds columns
        pub.subscribe(self.enable_tax, 'receive_tax')
        pub.subscribe(self.receive_meta, 'receive_metadata')
        pub.subscribe(self.check_networks, 'update_settings')
        pub.subscribe(self.check_database, 'database_log')
        pub.subscribe(self.update_pid, 'pid')

        self.settings = general_settings
        self.agglom = None
        self.agglom_weight = None
        self.logic = None
        self.assoc = None
        self.networks = None
        self.export = None
        self.checks = str()
        self.address = ['bolt://localhost:7687']
        self.username = ['neo4j']
        self.password = ['neo4j']
        self.neo4j = None
        self.process = None
        self.gml_name = ['network']

        # defines columns
        self.leftsizer = wx.BoxSizer(wx.VERTICAL)
        self.rightsizer = wx.BoxSizer(wx.VERTICAL)
        self.topsizer = wx.BoxSizer(wx.HORIZONTAL)

        # Opening neo4j folder
        self.neo_btn = wx.Button(self, label="Select Neo4j folder", size=btnsize)
        self.neo_btn.Bind(wx.EVT_BUTTON, self.open_neo)
        self.neo_btn.Bind(wx.EVT_MOTION, self.update_help)

        self.address_txt = wx.StaticText(self, label='Neo4j database address')
        self.address_box = wx.TextCtrl(self, value='bolt://localhost:7687', size=btnsize)
        self.username_txt = wx.StaticText(self, label='Neo4j username & password')
        self.username_box = wx.TextCtrl(self, value='neo4j', size=btnsize)
        self.pass_box = wx.TextCtrl(self, value='neo4j', size=btnsize)
        self.address_box.Bind(wx.EVT_TEXT, self.update_address)
        self.address_txt.Bind(wx.EVT_MOTION, self.update_help)
        self.address_box.Bind(wx.EVT_MOTION, self.update_help)
        self.username_box.Bind(wx.EVT_TEXT, self.update_username)
        self.username_box.Bind(wx.EVT_MOTION, self.update_help)
        self.pass_box.Bind(wx.EVT_TEXT, self.update_pass)
        self.pass_box.Bind(wx.EVT_MOTION, self.update_help)

        # set up database
        self.data_button = wx.Button(self, label='Launch database', size=btnsize)
        self.data_button.Bind(wx.EVT_MOTION, self.update_help)
        self.data_button.Bind(wx.EVT_BUTTON, self.start_database)

        # open database in browser
        self.data_browser = wx.Button(self, label='Open database in browser', size=btnsize)
        self.data_browser.Bind(wx.EVT_MOTION, self.update_help)
        self.data_browser.Bind(wx.EVT_BUTTON, self.open_browser)
        self.data_browser.Enable(False)

        # close database
        self.close_button = wx.Button(self, label='Close database', size=btnsize)
        self.close_button.Bind(wx.EVT_MOTION, self.update_help)
        self.close_button.Bind(wx.EVT_BUTTON, self.close_database)

        # select taxonomic levels
        self.tax_txt = wx.StaticText(self, label='Agglomerate edges to: ')
        self.tax_choice = wx.CheckListBox(self, choices=['Species', 'Genus',
                                                         'Family', 'Order', 'Class', 'Phylum'],
                                          size=(boxsize[0], 110), style=wx.LB_MULTIPLE)
        self.tax_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.tax_choice.Bind(wx.EVT_CHECKLISTBOX, self.get_levels)
        self.tax_choice.Enable(False)
        self.tax_txt.Enable(False)

        # button for agglomeration
        self.weight_txt = wx.StaticText(self, label='During network agglomeration:')
        self.tax_weight = wx.ListBox(self, choices=['Take weight into account', 'Ignore weight'], size=(boxsize[0], 40))
        self.tax_weight.Bind(wx.EVT_MOTION, self.update_help)
        self.tax_weight.Bind(wx.EVT_LISTBOX, self.weight_agglomeration)
        self.tax_weight.Enable(False)
        self.weight_txt.Enable(False)


        # button for sample association
        self.assoc_txt = wx.StaticText(self, label='Associate taxa to: ' )
        self.assoc_box = wx.CheckListBox(self, choices=['Select a BIOM or metadata file first'], size=(300, 90))
        self.assoc_box.Bind(wx.EVT_MOTION, self.update_help)
        self.assoc_box.Bind(wx.EVT_CHECKLISTBOX, self.run_association)
        self.assoc_txt.Enable(False)
        self.assoc_box.Enable(False)

        # logic operations
        self.logic_txt = wx.StaticText(self, label='For networks, perform the following logic operations...')
        self.logic_choice = wx.ListBox(self, choices=['None', 'Union', 'Intersection', 'Difference'],
                                          size=(boxsize[0], 70))
        self.logic_choice.Bind(wx.EVT_MOTION, self.update_help)
        self.logic_choice.Bind(wx.EVT_LISTBOX, self.get_logic)
        self.logic_choice.Enable(False)
        self.logic_txt.Enable(False)

        # review pane
        # review settings
        self.rev_text = wx.StaticText(self, label='Current operations')
        self.review = wx.TextCtrl(self, value='', size=(boxsize[0], 300), style=wx.TE_READONLY | wx.TE_MULTILINE)

        # Run button
        self.go = wx.Button(self, label='Run database operations', size=(btnsize[0], 40))
        self.go.Bind(wx.EVT_BUTTON, self.run_database)
        self.go.SetFont(wx.Font(16, wx.DECORATIVE, wx.NORMAL, wx.BOLD))
        self.go.SetBackgroundColour(wx.Colour(0, 153, 51))

        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.address_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.leftsizer.Add(self.address_box, flag=wx.ALIGN_LEFT)
        self.leftsizer.AddSpacer(10)
        self.leftsizer.Add(self.username_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.leftsizer.Add(self.username_box, flag=wx.ALIGN_LEFT)
        self.leftsizer.Add(self.pass_box, flag=wx.ALIGN_LEFT)
        self.leftsizer.Add(self.neo_btn, flag=wx.ALIGN_LEFT)
        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.data_button, flag=wx.ALIGN_LEFT)
        self.leftsizer.Add(self.close_button, flag=wx.ALIGN_LEFT)
        self.leftsizer.Add(self.data_browser, flag=wx.ALIGN_LEFT)
        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.rev_text, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.leftsizer.Add(self.review, flag=wx.ALIGN_CENTER)


        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.tax_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.tax_choice, flag=wx.ALIGN_LEFT)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.weight_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.tax_weight, flag=wx.ALIGN_LEFT)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.assoc_txt, flag=wx.ALIGN_CENTER_HORIZONTAL)
        self.rightsizer.Add(self.assoc_box, flag=wx.ALIGN_LEFT)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.logic_txt)
        self.rightsizer.Add(self.logic_choice)
        self.rightsizer.AddSpacer(20)
        self.rightsizer.Add(self.go, flag=wx.ALIGN_CENTER)
        self.rightsizer.AddSpacer(20)

        self.topsizer.AddSpacer(10)
        self.topsizer.Add(self.leftsizer)
        self.topsizer.AddSpacer(10)
        self.topsizer.Add(self.rightsizer)

        self.SetSizerAndFit(self.topsizer)
        self.Fit()

        # help strings for buttons
        self.buttons = {self.pass_box: 'Supply password for Neo4j database.'
                                       'For details on configuring your database, check the Neo4j manual.',
                        self.address_txt: 'Supply address of Neo4j database.'
                                          'For details on configuring your database, check the Neo4j manual.',
                        self.address_box: 'Supply address of Neo4j database.'
                                          'For details on configuring your database, check the Neo4j manual.',
                        self.username_txt: 'Supply username for Neo4j database.'
                                           'For details on configuring your database, check the Neo4j manual.',
                        self.username_box: 'Supply username for Neo4j database.'
                                           'For details on configuring your database, check the Neo4j manual.',
                        self.data_button: 'Launch local Neo4j database.',
                        self.close_button: 'Shut down local Neo4j database.',
                        self.neo_btn: 'Location of your Neo4j folder.',
                        self.data_browser: 'Open Neo4j Browser and explore your data.',
                        self.tax_choice: 'Associations that are taxonomically similar at the specified level are '
                                         'combined into agglomerated associations. ',
                        self.tax_weight: 'If selected, only edges with matching weight are agglomerated. ',
                        self.assoc_box: "Taxa are linked to categorical variables through a hypergeometric test"
                                           " and to continous variables through Spearman's rank correlation.",
                        self.logic_choice: 'Find associations that are present in only one or all of your networks.',
                        self.go: 'Run the selected operations and export a GraphML file.'
                        }


    def update_help(self, event):
        btn = event.GetEventObject()
        if btn in self.buttons:
            status = self.buttons[btn]
            pub.sendMessage('change_statusbar', msg=status)

    def open_neo(self, event):
        dlg = wx.DirDialog(self, "Select Neo4j directory", style=wx.DD_DEFAULT_STYLE)
        if dlg.ShowModal() == wx.ID_OK:
            self.neo4j = list()
            self.neo4j.append(dlg.GetPath())
        self.send_settings()
        dlg.Destroy()

    def update_address(self, event):
        text = self.address_box.GetValue()
        if self.address is None:
            self.address = list()
            self.address.append(text)
        else:
            self.address[0] = text
        self.send_settings()

    def update_username(self, event):
        text = self.username_box.GetValue()
        if self.username is None:
            self.username = list()
            self.username.append(text)
        else:
            self.username[0] = text
        self.send_settings()

    def update_pass(self, event):
        text = self.pass_box.GetValue()
        if self.password is None:
            self.password = list()
            self.password.append(text)
        else:
            self.password[0] = text
        self.send_settings()


    def update_pid(self, msg):
        """Listener for Neo4j PID."""
        self.process = msg

    def start_database(self, event):
        checks = str()
        try:
            eg = Thread(target=data_starter, args=(self.settings,))
            eg.start()
            eg.join()
        except Exception:
            logger.error("Failed to initiate database", exc_info=True)
        # removed dlg.LoadingBar() dlg.ShowModal()
        self.tax_txt.Enable(True)
        self.tax_choice.Enable(True)
        self.tax_weight.Enable(True)
        self.data_browser.Enable(True)
        self.assoc_txt.Enable(True)
        self.assoc_box.Enable(True)
        self.logic_txt.Enable(True)
        self.logic_choice.Enable(True)
        self.weight_txt.Enable(True)

    def close_database(self, event):
        try:
            # there is a lingering Java process that places a lock on the database.
            # terminating the subprocess does NOT terminate the Java process,
            # so the store lock has to be deleted manually.
            # This is different for Linux & Windows machines and may not be trivial
            # however, PID solution may be platform-independent
            # CURRENT SOLUTION:
            # get parent PID of subprocess
            # use psutil to get child PIDs
            # kill child PIDs too
            parent_pid = self.process
            parent = Process(parent_pid)
            children = parent.children(recursive=True)
            for child in children:
                child.kill()
            # apparently killing the children also kills the parent
        except Exception:
            logger.error("Failed to close database", exc_info=True)

    def open_browser(self, event):
        url = "http://localhost:7474/browser/"
        webbrowser.open(url)

    def enable_tax(self, msg):
        self.tax_choice.Enable(True)

    def get_levels(self, event):
        text = self.tax_choice.GetCheckedStrings()
        self.agglom = list(text)
        self.send_settings()

    def weight_agglomeration(self, event):
        self.agglom_weight = list()
        name = self.tax_weight.GetSelection()
        text = self.tax_weight.GetString(name)
        self.agglom_weight.append(text)
        self.send_settings()

    def run_association(self, event):
        self.assoc = list()
        text = list()
        ids = self.assoc_box.GetSelections()
        for i in ids:
            text.append(self.assoc_box.GetString(i))
        self.assoc = text
        self.send_settings()

    def get_logic(self, event):
        self.logic = list()
        name = self.logic_choice.GetSelection()
        text = self.logic_choice.GetString(name)
        self.logic.append(text)
        self.send_settings()

    def receive_meta(self, msg):
        """
        Listener function that registers whether a BIOM file with metadata
        or separate metadata file has been supplied.
        It opens the file to check whether there are variables present,
        and updates the split_list listbox to include these variables.
        """
        self.assoc_box.Set(msg)

    def send_settings(self):
        """
        Publisher function for settings
        """
        settings = {'assoc': self.assoc, 'agglom': self.agglom,
                    'logic': self.logic, 'agglom_weight': self.agglom_weight,
                    'export': self.export,
                    'address': self.address, 'username': self.username,
                    'password': self.password, 'neo4j': self.neo4j, 'gml_name': self.gml_name}
        pub.sendMessage('update_settings', msg=settings)

    def run_database(self, event):
        try:
            eg = Thread(target=data_worker, args=(self.settings,))
            eg.start()
        except Exception:
            logger.error("Failed to start database worker", exc_info=True)
        # removed LoadingBar()


    def check_networks(self, msg):
        # define how files should be checked for, it is important that import functions work!
        try:
            if msg['network']:
                for file in msg['network']:
                    try:
                        network = nx.read_weighted_edgelist(file)
                        self.checks += "Loaded network from " + file + ". \n\n"
                        nodes = len(network.nodes)
                        edges = len(network.edges)
                        self.checks += "This network has " + str(nodes) + \
                                       " nodes and " + str(edges) + " edges. \n\n"
                    except (TypeError, NetworkXError):
                        wx.LogError("Could not import edge list '%s'." % file)
                        logger.error("Could not import edge list", exc_info=True)
                    try:
                        weight = nx.get_edge_attributes(network, 'weight')
                        if len(weight) > 0:
                            self.checks += 'This is a weighted network. \n\n'
                        else:
                            self.checks += 'This is an unweighted network. \n\n'
                    except (TypeError, NetworkXError):
                        wx.LogError("Could not access edge metadata '%s'." % file)
                        logger.error("Could not access edge metadata", exc_info=True)
                    # need to check if network IDs match IDs in the BIOM files.
                    allbioms = list()
                    allbioms.extend(msg['procbioms'])
                    allbioms.extend(msg['otu_table'])
                    allbioms.extend(msg['biom_file'])
                    match = 0
                    taxa = None
                    for biomfile in allbioms:
                        try:
                            biomtab = load_table(biomfile)
                            taxa = biomtab.ids(axis='observation')
                        except TypeError:
                            wx.LogError("Could not access source BIOM file '%s'." % file)
                            logger.error("Could not access source BIOM file", exc_info=True)
                        if taxa:
                            nodes = list(network.nodes)
                            if all(elem in taxa for elem in nodes):
                                match += 1
                                self.checks += 'Node identifiers in ' + biomfile + \
                                               ' matched node identifiers in ' + file + '. \n\n'
                    if match == 0:
                        wx.LogError("No BIOM file matched network nodes!")
                        logger.error("No BIOM file matched network nodes!", exc_info=True)
                self.review.SetValue(self.checks)
        except KeyError:
            pass

    def check_database(self, msg):
        self.checks += msg
        self.review.SetValue(self.checks)


class LoadingBar(wx.Dialog):
    def __init__(self):
        """Constructor"""
        wx.Dialog.__init__(self, None, title="Progress")
        self.count = 0
        self.text = wx.StaticText(self, label="Starting...")
        self.progress = wx.Gauge(self, range=4)
        sizer = wx.BoxSizer(wx.VERTICAL)
        sizer.Add(self.text, 0, wx.EXPAND | wx.ALL, 10)
        sizer.Add(self.progress, 0, wx.EXPAND | wx.ALL, 10)
        self.SetSizer(sizer)
        pub.subscribe(self.get_progress, "update")

    def get_progress(self, msg):
        """This is not working yet! Find a better way to track progress of the thread... """
        testlist = ["Agglomerating edges...", "Associating samples...",
                    "Uploading BIOM files...", "Uploading network files...",
                    "Processing associations...", "Exporting network..."]
        if msg in testlist:
            self.count += 1
        self.progress.SetValue(self.count)
        self.text.SetLabel(msg)
        if msg == 'Completed database operations!':
            self.progress.SetValue(4)
            sleep(3)
            self.Destroy()


def data_worker(inputs):
    """
    Carries out operations on database as specified by user.
    """
    try:
        pub.sendMessage('update', msg='Starting database drivers.')
        importdriver = ImportDriver(user=inputs['username'][0],
                                    password=inputs['password'][0],
                                    uri=inputs['address'][0])
        statdriver = Driver(user=inputs['username'][0],
                            password=inputs['password'][0],
                            uri=inputs['address'][0])
    except Exception:
        logger.error("Failed to start database worker. ", exc_info=True)
    try:
        # write operations here
        if inputs['agglom']:
            tax_list = ['Species', 'Genus', 'Family', 'Order', 'Class', 'Phylum', 'Kingdom']
            level_id = tax_list.index(inputs['agglom'][0])
            if inputs['agglom_weight'][0]:
                mode = inputs['agglom_weight'][0]
            else:
                mode = 'Ignore weight'
            for level in range(0, level_id+1):
                pub.sendMessage('update', msg="Agglomerating edges...")
                statdriver.agglomerate_network(level=tax_list[level], mode=mode)
        if inputs['assoc']:
            pub.sendMessage('update', msg="Associating samples...")
            statdriver.associate_samples(properties=inputs['assoc'])
        pairlist = None
        if inputs['logic']:
            if inputs['logic'][0] == 'Union':
                pairlist = statdriver.graph_union()
            if inputs['logic'][0] == 'Intersection':
                pairlist = statdriver.graph_intersection()
            if inputs['logic'][0] == 'Difference':
                pairlist = statdriver.graph_difference()
        filename = inputs['fp'][0] + '/' + inputs['gml_name'][0] + '.graphml'
        pub.sendMessage('update', msg="Exporting network...")
        importdriver.export_network(path=filename, pairlist=pairlist, mode='basic')
        pub.sendMessage('update', msg="Completed database operations!")
    except Exception:
        logger.error("Failed to run database worker", exc_info=True)
    try:
        logfile = open(resource_path("massoc.log"), 'r')
        logtext = logfile.read()
        logfile.close()
        dump = open(inputs['fp'][0], 'w')
        dump.write(logtext)
        dump.close()
    except Exception:
        pass


def data_starter(inputs):
    """Starts up database and uploads specified files. """
    checks = str()
    try:
        pub.sendMessage('update', msg='Starting database.')
        # run database
        if system() == 'Windows':
            filepath = inputs['neo4j'][0] + '/bin/neo4j.bat console'
        else:
            filepath = inputs['neo4j'][0] + '/bin/neo4j console'
        filepath = filepath.replace("\\", "/")
        if system() == 'Windows':
            p = Popen(filepath, shell=True)
        else:
            p = Popen(["gnome-terminal", "-e", filepath])  # x-term compatible alternative terminal
        pub.sendMessage('pid', msg=p.pid)
    except Exception:
        logger.error("Failed to initiate database", exc_info=True)
    i = 0
    importdriver = None
    while not importdriver and i < 10:
        sleep(12)
        importdriver = ImportDriver(user=inputs['username'][0],
                                    password=inputs['password'][0],
                                    uri=inputs['address'][0])
        i += 1
    if i == 10:
        logger.error("Unable to access Neo4j database.", exc_info=True)
        pub.sendMessage('update', msg='Could not access Neo4j database!')
    importdriver.clear_database()
    try:
        pub.sendMessage('update', msg='Uploading BIOM files...')
        itemlist = list()
        for item in inputs['procbioms']:
            biomfile = load_table(item)
            importdriver.convert_biom(biomfile=biomfile, exp_id=item)
            itemlist.append(item)
    except Exception:
        checks += 'Unable to upload BIOM files to Neo4j database.'
        logger.error("Unable to upload BIOM files to Neo4j database.", exc_info=True)
    try:
        pub.sendMessage('update', msg='Uploading network files...')
        for item in inputs['network']:
            network = nx.read_weighted_edgelist(item)
            importdriver.convert_networkx(network=network, network_id=item, mode='weight')
            itemlist.append(item)
        checks += 'Successfully uploaded the following items and networks to the database: \n'
        for item in itemlist:
            checks += (item + '\n')
        checks += '\n'
    except Exception:
        checks += 'Unable to upload network files to Neo4j database.'
        logger.error("Unable to upload network files to Neo4j database.", exc_info=True)
    pub.sendMessage('update', msg='Completed database operations!')
    pub.sendMessage('database_log', msg=checks)


if __name__ == "__main__":
    app = wx.App(False)
    app.MainLoop()
