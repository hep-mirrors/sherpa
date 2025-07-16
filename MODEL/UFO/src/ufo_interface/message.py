try:
    import colorama
    colorama.init()
    RESET = colorama.Fore.RESET
    GREEN = colorama.Fore.GREEN
    YELLOW = colorama.Fore.YELLOW
    BLUE = colorama.Fore.BLUE
    RED = colorama.Fore.RED
    colorama.deinit()
except ImportError:
    RESET = "\033[0m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    RED = "\033[31m"

from sys import argv
from os import path

EXEC_STR = path.basename(argv[0])+": "
GREEN_INDENT = GREEN+EXEC_STR+RESET
BLUE_INDENT = BLUE+EXEC_STR+RESET
YELLOW_INDENT = YELLOW+EXEC_STR+RESET
RED_INDENT = RED+EXEC_STR+RESET
INDENT = len(EXEC_STR)*" "


def color_print(string, color):
    lines = [line.rstrip() for line in string.split("\n")]
    if len(lines) == 0:
        lines = [""]
    print(color+lines[0])
    for line in lines[1:]:
        print(INDENT+line+RESET)


def error(string):
    color_print(string, RED_INDENT)


def warning(string):
    color_print(string, YELLOW_INDENT)


def progress(string):
    color_print(string, GREEN_INDENT)
