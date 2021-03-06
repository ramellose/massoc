#!/usr/bin/env python

"""
This interface covers all of massoc's main features.
It allows users to select appropriate settings and export these as the appropriate command line call.
Moreover, it incorporates checks to make sure supplied files are correct.
By visualizing input and output, it provides interactive feedback to users that helps them in their
decision-making process.

To do: write a safe close button that also terminates R script processes.
Right now, the R scripts keep running even if you quit running massoc.
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import wx
import massoc
from massoc.scripts.main import resource_path
from wx.lib.pubsub import pub
from massoc.GUI.intro import IntroPanel
from massoc.GUI.input import InputPanel
from massoc.GUI.process import ProcessPanel
from massoc.GUI.network import NetworkPanel
from massoc.GUI.database import DataPanel
from massoc.GUI.analysis import AnalysisPanel
import multiprocessing
# source: https://stackoverflow.com/questions/4004353/logging-strategy-for-gui-program

general_settings = {"biom_file": None,
                     "otu_table": None,
                     "tax_table": None,
                     "sample_data": None,
                     "otu_meta": None,
                     "cluster": None,
                     "split": None,
                     "prev": 20,
                     "fp": None,
                     "levels": None,
                     "tools": None,
                     "spiec": None,
                     "conet": None,
                     "conet_bash": None,
                     "spiec_settings": None,
                     "spar": None,
                     "spar_pval": None,
                     "spar_boot": None,
                     "nclust": None,
                     "name": None,
                     "cores": None,
                     "rar": None,
                     "min": None,
                     "network": None,
                     "assoc": None,
                     "agglom": None,
                     "logic": None,
                     "agglom_weight": None,
                     "export": None,
                     "neo4j": None,
                     "procbioms": None,
                     "address": "bolt://localhost:7687",
                     "username": "neo4j",
                     "password": "neo4j",
                     "variable": None,
                     "weight": None,
                     "networks": None,
                     "output": None,
                     "add": None}


class BuildFrame(wx.Frame):
    """Constructor"""
    def __init__(self):
        wx.Frame.__init__(self, None, title='massoc', size=(800, 700))

        ico = wx.Icon(resource_path("massoc.png"), wx.BITMAP_TYPE_PNG)
        self.SetIcon(ico)

        p = wx.Panel(self)
        self.nb = wx.Notebook(p)
        self.tab1 = IntroPanel(self.nb)
        self.tab2 = InputPanel(self.nb)
        self.tab3 = ProcessPanel(self.nb)
        self.tab4 = NetworkPanel(self.nb)
        self.tab5 = DataPanel(self.nb)
        self.tab6 = AnalysisPanel(self.nb)

        self.nb.AddPage(self.tab1, "Start")
        self.nb.AddPage(self.tab2, "Input files")
        self.nb.AddPage(self.tab3, "Preprocessing")
        self.nb.AddPage(self.tab4, "Network inference")
        self.nb.AddPage(self.tab5, "Network database")
        self.nb.AddPage(self.tab6, "Network analysis")

        self.settings = general_settings

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        p.SetSizer(sizer)

        # listens to help messages from uncoupled tab files
        self.CreateStatusBar()
        pub.subscribe(self.change_statusbar, 'change_statusbar')
        self.Show()
        pub.subscribe(self.format_settings, 'input_settings')
        pub.subscribe(self.format_settings, 'process_settings')
        pub.subscribe(self.format_settings, 'network_settings')
        pub.subscribe(self.format_settings, 'data_settings')
        pub.subscribe(self.format_settings, 'analysis_settings')
        pub.subscribe(self.load_settings, 'load_settings')

    def format_settings(self, msg):
        """
        Listener function for settings from tabs in notebook.
        """
        try:
            for key in msg:
                self.settings[key] = msg[key]
        except:
            pass
        pub.sendMessage('show_settings', msg=self.settings)

    def load_settings(self, msg):
        try:
            for key in msg:
                self.settings[key] = msg[key]
        except:
            pass
        pub.sendMessage('show_settings', msg=self.settings)

    def change_statusbar(self, msg):
        self.SetStatusText(msg)


if __name__ == "__main__":
    multiprocessing.freeze_support()
    app = wx.App(False)
    frame = BuildFrame()
    app.MainLoop()
