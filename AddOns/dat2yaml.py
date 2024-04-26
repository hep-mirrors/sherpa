#!/usr/bin/env python3

import argparse
import re

import yaml

print('WARNING : This is intended as a first step to reduce manual work in the conversion.')
print('WARNING : The conversion will be incomplete in most cases, so use with care.')

parser = argparse.ArgumentParser(description='translate .dat runcard into new .yaml file')
parser.add_argument('-f',default='Run.dat',type=str,help='run card filepath',dest='rcpath')
parser.add_argument('-o',default='Sherpa.yaml',type=str,help='outfile filepath',dest='outpath')

args = parser.parse_args()

def qcd_order(order):
    if(order=='*'):
        return 'Any'
    else:
        return order

class runcard():
    def __init__(self, rcin, rcout):
        self.tags          = {'TAGS'          : {}}
        self.hard_decays   = {'HARD_DECAYS'   : {}}
        self.particle_data = {'PARTICLE_DATA' : {}}
        self.options       = {}
        self.process       = {'PROCESSES'     : []}
        self.selectors     = {'SELECTORS'     : []}
        self.datafile      = rcin
        self.outfile       = rcout
        self.uniformed     = self.uniform_dat_file()

        self.yaml_runcard  = {'tags'          : self.tags,
                              'hard_decays'   : self.hard_decays,
                              'particle_data' : self.particle_data,
                              'options'       : self.options,
                              'processes'     : self.process,
                              'selectors'     : self.selectors}

    def uniform_dat_file(self):
        '''
        Uniform Run.dat file to a single readable format.

        This function takes care of mis-conventions:
        notably more definitions on a single line, or non space after ';' in a single line.

        Parameters
        ----------

        Returns
        -------
        string
        returns uniformed runcard
        '''
        uniformed_text = ''
        for l in open(self.datafile):
            if(l[l.find(';')+1] != ' '):
                l = l.replace(';','; ')
            endendline = l.replace('; ',';\n').strip(' ')
            uniformed_text = uniformed_text+endendline
        return uniformed_text

    def get_tags(self,line):
        '''
        extract TAGNAME:=TAGVALUE from line

        extract TAGNAME:=TAGVALUE from line.

        Parameters
        ----------
        arg2 : str
        the line of a runcard

        Returns
        -------
        void
        appends to self.tags the founded tags.
        '''

        search_tags = re.compile('\:=', re.X)
        if(search_tags.search(line)):
            l_list = line.split(':=')
            tagname, tagval = l_list[0], l_list[1].rstrip(';')
            # the new format allows only for ranges and not
            # for comma-separated fields in LJETS
            if(tagname=='LJET'):
                if(tagval.find(',')!=-1):
                    ini,fnl  = tagval.split(',')[0], tagval.split(',')[-1]
                    tagval = ini+'-'+fnl
            self.tags['TAGS'][tagname] = tagval
        return

    def get_options(self, line):
        line = line.replace('\'','')
        try:
            oname = line.rstrip(';').split(' ')[0]
            if (oname=='CKKW' or oname=='CKKW-L'):
                return
            oval  = line.rstrip(';').split(' ')[1:]
            # if more than one option (like for ME_GENERATORS) take
            # full list, otherwise only first element
            if len(oval) > 1:
                if '=' in oval:
                   oval = oval[1]
                self.options[oname] = oval
            else:
                self.options[oname] = oval[0]

        except IndexError:
            oname = line.rstrip(';').split('=')[0]
            oval  = line.rstrip(';').split('=')[1:]
            if len(oval) > 1:
                if '=' in oval:
                   oval = oval[1]
                self.options[oname] = oval
            else:
                self.options[oname] = oval[0]

    def get_process(self):
        for proclist in re.findall('processes(.*?)processes', self.uniformed, re.S):
            for procline in re.findall('Process (.*?)End process', proclist, re.S):
                process_entry = {'dummy procdef': {}}
                procline = procline.replace('\n\n','\n').replace(';','')
                for l in iter(procline.splitlines()):
                    par_info = None
                    try:
                        par_info = [r for r in re.findall('\{(.*?)\}', l, re.S)][-1];
                    except IndexError:
                        pass
                    if l.find('->')!=-1:
                        procdef = l;
                        if(par_info):
                            procdef = l.replace(par_info,'$('+par_info+ ')')
                        process_entry[procdef] = process_entry['dummy procdef']
                        del process_entry['dummy procdef']
                        continue
                    if l.find('Order')!=-1:
                        # assumes no model coupling (like HEFT)
                        orders = l.split(' ')[1].strip('\(\)').split(',')
                        process_entry[list(process_entry.keys())[0]]['Order'] = {'QCD':qcd_order(orders[0]),'EW':orders[1]};
                        continue
                    if l.find('CKKW')!=-1:
                        process_entry[list(process_entry.keys())[0]]['CKKW'] = '$(QCUT)'
                        continue

                    # above are all the common options, all other
                    # options may be applied only to a given (or a range of)
                    # multiplicity, this info is saved in par_info.
                    if(par_info):
                        if(par_info == 'LJET'):
                            par_info = '$('+par_info+ ')'
                        # the new format allows only for ranges and not
                        # for comma-separated fields
                        if(par_info.find(',')!=-1):
                            ini,fnl  = par_info.split(',')[0], par_info.split(',')[-1]
                            par_info = ini+'-'+fnl
                        entry = l.split(' ')[1]
                        if(entry in self.tags['TAGS']):
                            entry = '$('+entry+ ')'
                        try:
                            process_entry[list(process_entry.keys())[0]][par_info][l.split(' ')[0]] = entry
                        except (KeyError,AttributeError):
                            process_entry[list(process_entry.keys())[0]][par_info] = dict({l.split(' ')[0]:entry})

                        continue
                    else:
                        entry = l.split(' ')[1]
                        if(entry in self.tags['TAGS']):
                            entry = '$('+entry+ ')'
                        process_entry[list(process_entry.keys())[0]][l.split(' ')[0]] = entry
                self.process['PROCESSES'].append(process_entry)

    def get_selectors(self):
        for selline in re.findall('selector(.*?)selector', self.uniformed, re.S):
            print('WARNING : Selectors have for now to be done manually')

    def contruct_particle_data(self):
        '''
        Look inside loaded options, and construct particle data info

        Search for MASSIVE, MASS, ACTIVE, WIDTH or STABLE keywords and
        contruct the relevant particle data info

        Parameters
        ----------

        Returns
        -------
        void
        fills self.particle_data
        '''
        pd_keywords = 'MASSIVE MASS ACTIVE WIDTH STABLE'.split(' ')
        pnumbers = set([r.lstrip('[').rstrip(']') for r in re.findall('\[\d+\]', self.uniformed, re.S)]);
        if(not pnumbers):
            return

        for pn in pnumbers:
            self.particle_data['PARTICLE_DATA'][pn] = {}
        options_to_be_removed = []
        for o in self.options:
            if(o.startswith(tuple(pd_keywords))):
                onumber    = [r for r in re.findall('\[\d+\]', o, re.S)][0]
                p_property = o.replace(onumber, '').split(' ')[0][0] + o.replace(onumber, '').split(' ')[0][1:].lower()
                p_value    = self.options[o]
                onumber    = onumber.lstrip('[').rstrip(']')
                self.particle_data['PARTICLE_DATA'][onumber][p_property] = p_value
                options_to_be_removed.append(o)

        for o in options_to_be_removed:
            self.options.pop(o)

    def translate_beams(self):
        try:
            self.options['BEAMS']         = [self.options['BEAM_1'],
                                             self.options['BEAM_2']]
            self.options['BEAM_ENERGIES'] = [self.options['BEAM_ENERGY_1'],
                                             self.options['BEAM_ENERGY_2']]
            old_beams = 'BEAM_1 BEAM_2 BEAM_ENERGY_1 BEAM_ENERGY_2'.split(' ')
            for o in old_beams:
                self.options.pop(o)
        except KeyError:
            self.options['BEAMS']         = [self.options['BEAM_1'][0],
                                             self.options['BEAM_2'][0]]
            self.options['BEAM_ENERGIES'] = [self.options['BEAM_1'][1],
                                             self.options['BEAM_2'][1]]
            old_beams = 'BEAM_1 BEAM_2'.split(' ')
            for o in old_beams:
                self.options.pop(o)

    def construct_hard_decays(self):
        # this may not be exhaustive
        hdc_old_options = 'HARD_DECAYS HARD_MASS_SMEARING'.split(' ')
        hdc_new_options = 'Enabled Mass_Smearing'.split(' ')

        for i_hdc in range(len(hdc_old_options)):
            try:
                o_hdc = hdc_old_options[i_hdc]
                n_hdc = hdc_new_options[i_hdc]
                self.hard_decays['HARD_DECAYS'][n_hdc] = self.options[o_hdc]
                self.options.pop(o_hdc)
            except KeyError:
                return

    def translate(self):
        '''
        translate the loaded runcard

        Extended description of function.

        Parameters
        ----------
        arg1 : int
        Description of arg1
        arg2 : str
        Description of arg2

        Returns
        -------
        int
        Description of return value
        '''

        all_capital_regex = re.compile('^[A-Z]+(?:_[A-Z,0-9]+)*[\s,\[,=]',re.X)
        for l in iter(self.uniformed.splitlines()):
            self.get_tags(l)
            if(all_capital_regex.search(l)):
                self.get_options(l);
        # translate beams in new format
        self.translate_beams()
        # Extract process info
        self.get_process()
        # Extract selector info
        self.get_selectors()
        # Extract all other options (that are all-capital words + underscore combinations and or numbers)
        self.contruct_particle_data()
        self.construct_hard_decays()
        with open(self.outfile,'w') as out:
            for it in self.yaml_runcard:
                yaml.safe_dump(self.yaml_runcard[it],stream=out,indent=4, allow_unicode=False,default_flow_style=False)


rc = runcard(args.rcpath, args.outpath)
rc.translate()
