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
from wx.lib.pubsub import pub
from massoc.GUI.intro import IntroPanel
from massoc.GUI.input import InputPanel
from massoc.GUI.process import ProcessPanel
from massoc.GUI.network import NetworkPanel


class BuildFrame(wx.Frame):
    """Constructor"""
    def __init__(self):
        wx.Frame.__init__(self, None, title='massoc', size=(750, 700))

        ico = wx.Icon("massoc.png", wx.BITMAP_TYPE_PNG)
        self.SetIcon(ico)

        p = wx.Panel(self)
        nb = wx.Notebook(p)
        self.tab1 = IntroPanel(nb)
        self.tab2 = InputPanel(nb)
        self.tab3 = ProcessPanel(nb)
        self.tab4 = NetworkPanel(nb)

        nb.AddPage(self.tab1, "Start")
        nb.AddPage(self.tab2, "Input files")
        nb.AddPage(self.tab3, "Preprocessing")
        nb.AddPage(self.tab4, "Network inference")

        self.settings = {'biom_file': None, 'otu_table': None, 'tax_table': None, 'sample_data': None,
                         'otu_meta': None, 'cluster': None, 'split': None, 'prev': None, 'fp': None,
                         'levels': None, 'tools': None, 'spiec': None, 'conet': None, 'spar': None, 'spar_pval': None,
                         'spar_boot': None, 'nclust': None, 'name': None, 'cores': None, 'rar': None, 'min': None}

        sizer = wx.BoxSizer()
        sizer.Add(nb, 1, wx.EXPAND)
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
