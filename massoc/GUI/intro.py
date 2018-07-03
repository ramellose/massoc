"""
The introduction panel contains contact information and provides an interface for HTML text files.
These HTML text files contain a FAQ and step-by-step tutorial for running massoc. 
"""

__author__ = 'Lisa Rottjers'
__email__ = 'lisa.rottjers@kuleuven.be'
__status__ = 'Development'
__license__ = 'Apache 2.0'

import os

import wx
import wx.html
import wx.lib.wxpTag

import massoc


class IntroPanel(wx.Panel):
    def __init__(self, parent):
        wx.Panel.__init__(self, parent)
        self.frame = parent
        self.leftsizer   = wx.BoxSizer(wx.HORIZONTAL)
        self.topsizer = wx.BoxSizer(wx.VERTICAL)

        btnsize = (150, -1)

        path = (os.path.dirname(massoc.__file__) + '/docs/welcome.txt')
        path = path.replace('\\', '/')
        welcome_file = open(path, 'r')
        welcome = welcome_file.read()
        welcome_file.close()

        path = (os.path.dirname(massoc.__file__) + '/docs/intro.txt')
        path = path.replace('\\', '/')
        intro_file = open(path, 'r')
        intro = intro_file.read()
        intro_file.close()

        path = (os.path.dirname(massoc.__file__) + '/docs/input.txt')
        path = path.replace('\\', '/')
        input_file = open(path, 'r')
        input = input_file.read()
        input_file.close()

        path = (os.path.dirname(massoc.__file__) + '/docs/proc.txt')
        path = path.replace('\\', '/')
        proc_file = open(path, 'r')
        proc = proc_file.read()
        proc_file.close()

        path = (os.path.dirname(massoc.__file__) + '/docs/net.txt')
        path = path.replace('\\', '/')
        net_file = open(path, 'r')
        net = net_file.read()
        net_file.close()

        path = (os.path.dirname(massoc.__file__) + '/docs/demo.txt')
        path = path.replace('\\', '/')
        demo_file = open(path, 'r')
        demo = demo_file.read()
        demo_file.close()

        from massoc.docs.img1 import img1
        from massoc.docs.img2 import img2
        from massoc.docs.img3 import img3

        mfs = wx.MemoryFSHandler()
        wx.FileSystem.AddHandler(mfs)
        mfs.AddFile("input.png", img1.GetImage(), wx.BITMAP_TYPE_PNG)
        mfs.AddFile("process.png", img2.GetImage(), wx.BITMAP_TYPE_PNG)
        mfs.AddFile("network.png", img3.GetImage(), wx.BITMAP_TYPE_PNG)

        self.htmlbox = wx.html.HtmlWindow(self, -1, size=(800,1500))
        self.htmlbox.SetPage(welcome)
        self.topsizer.Add(self.htmlbox, 1, wx.ALIGN_CENTER_HORIZONTAL | wx.EXPAND | wx.ALL, 40)

        self.ico = wx.StaticBitmap(self, -1, wx.Bitmap("massoc_large.png", wx.BITMAP_TYPE_ANY))

        self.menusizer = wx.BoxSizer(wx.VERTICAL)
        self.menutitle = wx.StaticText(self, label="Documentation")
        font1 = wx.Font(18, wx.DECORATIVE, wx.NORMAL, wx.BOLD)
        self.menutitle.SetFont(font1)
        self.welcome_btn = wx.Button(self, label="Welcome to massoc", size=btnsize)
        self.welcome_btn.Bind(wx.EVT_BUTTON, self.change_html)
        self.intro_btn = wx.Button(self, label="Introduction", size=btnsize)
        self.intro_btn.Bind(wx.EVT_BUTTON, self.change_html)
        self.input_btn = wx.Button(self, label="Input files", size=btnsize)
        self.input_btn.Bind(wx.EVT_BUTTON, self.change_html)
        self.proc_btn = wx.Button(self, label="Preprocessing", size=btnsize)
        self.proc_btn.Bind(wx.EVT_BUTTON, self.change_html)
        self.net_btn = wx.Button(self, label="Network inference", size=btnsize)
        self.net_btn.Bind(wx.EVT_BUTTON, self.change_html)
        self.demo_btn = wx.Button(self, label="Demo", size=btnsize)
        self.demo_btn.Bind(wx.EVT_BUTTON, self.change_html)

        self.htmldict = {self.welcome_btn: welcome, self.intro_btn: intro, self.input_btn: input,
                         self.proc_btn: proc, self.net_btn: net, self.demo_btn: demo}

        self.menusizer.AddSpacer(50)
        self.menusizer.Add(self.ico, 1, wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.menusizer.AddSpacer(50)
        self.menusizer.Add(self.menutitle, 1, wx.EXPAND | wx.ALIGN_CENTER_HORIZONTAL, 5)
        self.menusizer.Add(self.welcome_btn, 1, wx.EXPAND, 5)
        self.menusizer.Add(self.intro_btn, 1, wx.EXPAND, 5)
        self.menusizer.Add(self.input_btn, 1, wx.EXPAND, 5)
        self.menusizer.Add(self.proc_btn, 1, wx.EXPAND, 5)
        self.menusizer.Add(self.net_btn, 1, wx.EXPAND, 5)
        self.menusizer.Add(self.demo_btn, 1, wx.EXPAND, 5)
        self.menusizer.AddSpacer(100)

        self.leftsizer.AddSpacer(20)
        self.leftsizer.Add(self.menusizer)
        self.leftsizer.AddSpacer(50)
        self.leftsizer.Add(self.topsizer)
        self.SetSizerAndFit(self.leftsizer)

    # write documentation in HTML and embed
    def change_html(self, event):
        name = event.GetEventObject()
        self.htmlbox.SetPage(self.htmldict[name])
