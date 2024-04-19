from importlib.machinery import SourceFileLoader
import argparse
import os
from .model_writer import ModelWriter
from .lorentz_writer import LorentzWriter
from pathlib import Path
import subprocess
from .message import error, warning, progress
import re
from . import s_color
import time


class Sherpa:
    def __init__(self):
        self.model_flags = '-g -O0 -fno-var-tracking'
        self.lorentz_flags = '-g -O2 -ffast-math'
        self.root_dir = '/home/isaacson/Documents/Projects/Sherpa/549-ufo-20/build/outputs/share/SHERPA-MC/'
        self.install_dir = self.root_dir + '/../../lib/SHERPA-MC/'


def parse_args():
    sherpa = Sherpa()
    parser = argparse.ArgumentParser(description='Generate a Sherpa model from a UFO model')
    parser.add_argument('ufo_path', type=str,
                        help='Path to the UFO model')
    parser.add_argument('--ncore', type=int, default=os.cpu_count(),
                        help='Number of cores to use')
    parser.add_argument('--root_dir', type=str,
                        default=sherpa.root_dir,
                        help='Path to Sherpa cmake config files')
    parser.add_argument('--output_dir', type=str,
                        default='.sherpa',
                        help='Path to write the Sherpa model')
    parser.add_argument('--install_dir', type=str,
                        default=sherpa.install_dir,
                        help='Path to Sherpa cmake config files')
    parser.add_argument('--model_flags', type=str,
                        default=sherpa.model_flags,
                        help='Flags to compile the model')
    parser.add_argument('--lorentz_flags', type=str,
                        default=sherpa.lorentz_flags,
                        help='Flags to compile the lorentz files')
    return parser.parse_args()


def check_color(color, particles):
    new_color = color
    pattern = r'Identity\(({}),(\d+)\)'
    replacement = r'Identity(\2,\1)'
    replacement_oct = r'IdentityG(\1,\2)'

    # Collect all anti-fundamental indices
    af_idxs = [i+1 for i, p in enumerate(particles) if p.color == -3]
    for idx in af_idxs:
        new_color = re.sub(pattern.format(idx), replacement, new_color)

    # Collect all octet indices
    oc_inds = [i+1 for i, p in enumerate(particles) if p.color == 8]
    for idx in oc_inds[::-1]:
        new_color = re.sub(pattern.format(idx), replacement_oct, new_color)

    return new_color


def check_model(model_name, ufo_mod):
    if model_name[0].isdigit():
        error(f'Model name {model_name} cannot start with a number')
        return False

    if model_name in ["SM", "HEFT", "TauPi"]:
        error(f'Model name {model_name} is reserved. Please rename UFO directory')
        return False

    return True


def prepare_output_dir(path):
    if os.path.exists(path):
        warning(f'Output directory {path} already exists.')
        # Ask if we should overwrite
        if input('Overwrite? (y/n): ').lower()[0] != 'y':
            exit(-1)
    else:
        os.makedirs(path)


def build_model(src_dir, install_dir, ncore):
    pwd = os.getcwd()
    os.chdir(src_dir)
    cmake_config_args = ['cmake', '-S', '.', '-B', 'build',
                         f'-DCMAKE_INSTALL_PREFIX={install_dir}']
    subprocess.run(cmake_config_args)
    cmake_build_args = ['cmake', '--build', 'build', '--', f'-j{ncore}']
    subprocess.run(cmake_build_args)
    cmake_install_args = ['cmake', '--install', 'build']
    subprocess.run(cmake_install_args)
    os.chdir(pwd)


def try_import(ufo_path):
    try:
        ufo_src = SourceFileLoader("ufo", ufo_path.as_posix())
        ufo_mod = ufo_src.load_module("ufo")
    except FileNotFoundError:
        error(f'UFO model not found at {ufo_path}.')
        error('Please check the path and try again.')
        return
    except ModuleNotFoundError:
        warning(f'UFO model at {ufo_path} may be written in python 2.')
        response = input('Would you like to convert to python3? (y/n): ')
        if response.lower()[0] == 'y':
            model_name = ufo_path.parts[-2] + '_py3'
            new_path = ufo_path.parent.parent / model_name
            print(f'Writing python3 model to {new_path}')
            # TODO: Update to different method for when 2to3 is removed
            subprocess.run(['2to3', '-W', ufo_path.parent.as_posix(),
                            '-o', new_path, '-n', '--no-diffs'])
            ufo_mod = try_import(new_path / '__init__.py')

    return ufo_mod


def main():
    args = parse_args()
    ufo_path = Path(args.ufo_path) / '__init__.py'
    ufo_mod = try_import(ufo_path)

    model_name = ufo_path.parts[-2]

    progress(f'Preparing output directory for {model_name}')
    out_dir = args.output_dir
    prepare_output_dir(out_dir)

    progress(f'Checking model {model_name}')
    # Check if the model is valid
    if not check_model(model_name, ufo_mod):
        return

    # Ensure color is in the correct format
    for vert in ufo_mod.all_vertices:
        for i, color in enumerate(vert.color):
            vert.color[i] = check_color(color, vert.particles)

    # Calculate time for color factor
    start = time.time()
    colors = set(sum([vert.color for vert in ufo_mod.all_vertices], []))
    colors = [color for color in colors if color != '1']
    color_files = []
    for color in colors:
        scolor = s_color.s_color(color)
        progress(f'Writing color file for color {scolor.name()}')
        scolor.write(out_dir)
        color_files.append(f'{scolor.name()}.C')
    end = time.time()
    progress(f'Color factor calculation took {end-start} seconds')

    start = time.time()
    progress(f'Writing lorentz files for {model_name}')
    lorentz_writer = LorentzWriter(out_dir)
    lorentz_writer.write_all(ufo_mod.all_lorentz, args.ncore)
    end = time.time()
    progress(f'Lorentz factor calculation took {end-start} seconds')

    opts = {
        'root_dir': args.root_dir,
        'install_dir': args.install_dir,
        'model_flags': args.model_flags,
        'lorentz_flags': args.lorentz_flags,
        'color_files': color_files,
    }
    progress(f'Writing model {model_name}')
    model_writer = ModelWriter(model_name, ufo_mod, opts)
    model_writer.write(out_dir)

    progress(f'Compiling model {model_name}')
    build_model(out_dir, args.install_dir, args.ncore)


if __name__ == '__main__':
    main()
