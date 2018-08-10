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
from massoc.scripts.main import resource_path
from wx.lib.pubsub import pub
from massoc.GUI.intro import IntroPanel
from massoc.GUI.input import InputPanel
from massoc.GUI.process import ProcessPanel
from massoc.GUI.network import NetworkPanel
from massoc.GUI.database import DataPanel
from massoc.scripts.main import general_settings

# source: https://stackoverflow.com/questions/4004353/logging-strategy-for-gui-program
import logging
import logging.handlers as handlers
logger = logging.getLogger()
hdlr = logging.FileHandler(resource_path("massoc.log"))
formatter = logging.Formatter('%(asctime)s %(name)-12s %(levelname)-8s %(message)s')
hdlr.setFormatter(formatter)
logger.addHandler(hdlr)
logger.setLevel(logging.WARNING)


class BuildFrame(wx.Frame):
    """Constructor"""
    def __init__(self):
        wx.Frame.__init__(self, None, title='massoc', size=(750, 700))

        ico = wx.Icon(resource_path("massoc.png"), wx.BITMAP_TYPE_PNG)
        self.SetIcon(ico)

        p = wx.Panel(self)
        self.nb = wx.Notebook(p)
        self.tab1 = IntroPanel(self.nb)
        self.tab2 = InputPanel(self.nb)
        self.tab3 = ProcessPanel(self.nb)
        self.tab4 = NetworkPanel(self.nb)
        self.tab5 = DataPanel(self.nb)

        self.nb.AddPage(self.tab1, "Start")
        self.nb.AddPage(self.tab2, "Input files")
        self.nb.AddPage(self.tab3, "Preprocessing")
        self.nb.AddPage(self.tab4, "Network inference")
        self.nb.AddPage(self.tab5, "Network database")

        self.settings = general_settings

        sizer = wx.BoxSizer()
        sizer.Add(self.nb, 1, wx.EXPAND)
        p.SetSizer(sizer)

        # listens to help messages from uncoupled tab files
        self.CreateStatusBar()
        pub.subscribe(self.change_statusbar, 'change_statusbar')
        self.Show()
        pub.subscribe(self.format_settings, 'update_settings')

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

    def change_statusbar(self, msg):
        self.SetStatusText(msg)


if __name__ == "__main__":
    app = wx.App(False)
    frame = BuildFrame()
    app.MainLoop()
