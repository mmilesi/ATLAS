"""
File    : Plotter/__init__.py
Authors : KG <Kong.Guan.Tan@cern.ch>

Folder should contain individual/institution plotting scripts/framework.
"""

# Intentionally blank
import optparse

def parseTauIDInputArgs():
    parser = optparse.OptionParser(description='TauID TemplateFit job configuration.')
    parser.add_option('-o', '--outdir', default=None, help="output dir")
    parser.add_option('-p', '--ptbin', default=None, help="select pt bin")
    parser.add_option('-e', '--etabin', default=None, help="select eta bin")
    parser.add_option('-t', '--bdtbin', default=None, help="select bdt bin")
    parser.add_option('-m', '--dotoys', default=False, help="run pseudo experiments", action="store_true")
    parser.add_option('-s', '--syst', default=None, help="systematics")
    parser.add_option('-d', '--systdir', default=None, help="systematics direction")
    parser.add_option('-w', '--mcsrshape', default=False, help="[For template fit robusness test] Reweight data driven W according to MC SR/WCR", action="store_true")
    parser.add_option('-b', '--batch', default=False, help="select batch mode", action="store_true")
    (options, args) = parser.parse_args()
    return options

# Create a class to allow inserting a newline in argparse help text:
#
# To use it, just import it from *this* module, and pass it as optional argument to ArgumentParser:
#
# parser = ArgumentParser(description='test', formatter_class=SmartFormatter)
#
# and specify a special intro for the options that should be handled "raw" (I use "R|rest of help"):
#
# Source : http://stackoverflow.com/questions/3853722/python-argparse-how-to-insert-newline-in-the-help-text

import argparse

class SmartFormatter(argparse.HelpFormatter):

    def _split_lines(self, text, width):
        if text.startswith('R|'):
            return text[2:].splitlines()
        # this is the RawTextHelpFormatter._split_lines
        return argparse.HelpFormatter._split_lines(self, text, width)
