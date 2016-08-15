try:
    import colorama
    colorama.init()
    reset  = colorama.Fore.RESET
    green  = colorama.Fore.GREEN
    yellow = colorama.Fore.YELLOW
    blue   = colorama.Fore.BLUE
    red    = colorama.Fore.RED
    colorama.deinit()
except ImportError as e:
    reset  = "\033[0m"
    green  = "\033[32m"
    yellow = "\033[33m"
    blue   = "\033[34m"
    red    = "\033[31m"

from sys import argv
from os import path

exec_str      = path.basename(argv[0])+": "
green_indent  = green+exec_str+reset
blue_indent   = blue+exec_str+reset
yellow_indent = yellow+exec_str+reset
red_indent    = red+exec_str+reset
indent        = len(exec_str)*" "

def error(string):
    lines = [l.rstrip() for l in string.split("\n")]
    if len(lines)==0: lines = [""]
    print red_indent+lines[0]
    for line in lines[1:]:
        print indent+line+reset

def warning(string):
    lines = [l.rstrip() for l in string.split("\n")]
    if len(lines)==0: lines = [""]
    print yellow_indent+lines[0]
    for line in lines[1:]:
        print indent+line+reset

def progress(string):
    lines = [l.rstrip() for l in string.split("\n")]
    if len(lines)==0: lines = [""]
    print green_indent+lines[0]
    for line in lines[1:]:
        print indent+line+reset
