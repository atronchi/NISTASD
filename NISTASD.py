#!/usr/bin/python

import urllib,HTMLParser
import re
from copy import copy
from math import isnan
import pylab as pl
from mpldatacursor import datacursor


class NISTASD(object):
    '''
    Import data from the NIST Atomic Spectra Database
    i.e. aspec = NISTASD(spec='He')
    '''    
    def __init__(self, spec = 'H', lowwl=10,uppwl=10000, order=1, verbose=False, plot=True):
        self.verbose = verbose

        self.spec = spec
        self.lowwl = lowwl
        self.uppwl = uppwl
        self.order = order

        self.get_asd()
        self.parse_asd()
        if plot: self.plot()


    def get_asd(self):
        spec = self.spec
        lowwl = self.lowwl/self.order
        uppwl = self.uppwl/self.order
       
        # build the web request
        self.nist_URL = 'http://physics.nist.gov/cgi-bin/ASD/lines1.pl'
        spec_plus=spec.strip().replace(' ','+') # HTML post needs + instead of ' '
        self.post_data = ('encodedlist=XXT1XXR0q0qVqVIII' + '&' # some key to make it work?
                        + 'spectra=' + spec_plus + '&' # eg 'He' or 'He+I' or 'He+II', no spaces
                        + 'low_wl=' + str(lowwl) + '&'
                        + 'upp_wl=' + str(uppwl) + '&'
                        + 'unit=0' + '&' # wl unit 0=Angstroms, 1=nm, 2=um
                        + 'en_unit=1' + '&' # energy unit 0 cm^-1, 1 eV, 2 Rydberg
                        + 'low_wn=' + '&'
                        + 'upp_wn=' + '&'
                        + 'submit=Retrieve+Data' + '&'
                        + 'temp=' + '&'
                        + 'doppler=' + '&'
                        + 'eden=' + '&'
                        + 'iontemp=' + '&'
                        + 'java_window=3' + '&'
                        + 'java_mult=' + '&'
                        + 'tsb_value=0' + '&'
                        + 'format=1' + '&' # 0 HTML output, 1 ascii output
                        + 'remove_js=on' + '&' # cleans up output for easier parsing
                        + 'output=0' + '&' # 0 return all output, 1 return output in pages
                        + 'page_size=15' + '&'
                        + 'line_out=0' + '&' # 0 return all lines, 1 only w/trans probs, 2 only w/egy levl, 3 only w/obs wls
                        + 'order_out=0' + '&' # output ordering: 0 wavelength, 1 multiplet
                        + 'show_av=2' + '&' # show wl in Vacuum (<2000A) Air (2000-20000A) Vacuum (>20,000A)
                        + 'max_low_enrg=' + '&' # maximum lower level energy
                        + 'max_upp_enrg=' + '&' # maximum upper level energy
                        + 'min_str=' + '&' # minimum transition strength
                        + 'max_str=' + '&' # maximum transition strength
                        + 'min_accur=' + '&' # minimum line accuracy, eg AAA AA A B C
                        + 'min_intens=' + '&' # minimum relative intensity to return
                        + 'show_obs_wl=1' + '&' # show observed wavelength
                        + 'show_calc_wl=1' + '&' # show calculated (Ritz) wavelength
                        + 'A_out=0' + '&' # show $
                        + 'intens_out=on' + '&' # show relative intensity
                        + 'allowed_out=1' + '&' # show allowed transitions
                        + 'forbid_out=1' + '&' # show forbidden transitions
                        + 'conf_out=on' + '&' # show electron configuration
                        + 'term_out=on' + '&' # show terms
                        + 'enrg_out=on' + '&' # show transition energies
                        + 'J_out=on' + '&' # show J (total angular momentum)
                        + 'g_out=on' ) # show g (statistical weight?)

        # issue wget to pull the data from nist and use sed to split off the desired info
        #  -q 'quiet' suppresses wget messages
        #  -O - directs results to standard output
        self.full_URL = self.nist_URL + '?' + self.post_data # This issues as a GET instead of POST, but it works ok anyway

        self.cmd = ( 'wget -q -O - \'' + self.full_URL + '\' '
                   + '| sed -n \'/<pre*/,/<\/pre>/p\' ' # select lines between <pre> tags
                   + '| sed \'/<*pre>/d\' ' # remove <pre> lines
                   + '| iconv -f ISO-8859-1 -t ASCII' ) # convert the web encoding to something IDL can understand...
                   # '| sed \'/----*/d\'' # remove ---- lines
        #sys.spawnl(cmd)

        self.nist_read = urllib.urlopen(self.full_URL).readlines()

        # select lines between <pre> tags as the asd_lines table
        self.asd_lines = []
        found_pre = False
        for ln in self.nist_read:
            if re.search('<.*?pre>',ln) != None:
                found_pre = not found_pre
                continue
            if found_pre:
                #self.asd_lines.append(ln)
                # convert ISO-8859-1 to ASCII or UTF-8 or unicode or something...
                self.asd_lines.append( HTMLParser.HTMLParser().unescape(ln) )
        if self.asd_lines == []:
            raise Exception('NoASDlines','No ASD lines were found.')


    # parse the imported asd_lines into data arrays
    def parse_asd(self):
        asd = copy(self.asd_lines)
        isec = -1
        self.header = []
        self.lines = []

        while len(asd) > 2:
            isec += 1
            self.parse_section(asd,isec)


    def parse_section(self,asd,isec = 0):
        # first do the header
        asd.pop(0) # first line is a break...
        if self.verbose:
            print asd[0]
        hd0 = [l.strip() for l in re.split(   '\|', asd[0])]
        hd0.pop() # last value is a line break
        idx = [i.start() for i in re.finditer('\|', asd[0])] # indices for the dividers
        idx.insert(0,0)

        asd.pop(0)
        if self.verbose:
            print asd[0]
        hd1 = [l.strip() for l in re.split(   '\|', asd[0])]
        hd1.pop()
        for i in range(0,len(hd0)):
            if hd0[i].find('level') == -1:
                hd0[i] += ' ' + hd1[i].strip()

        asd.pop(0)
        if self.verbose:
            print asd[0]
        hd = []
        for i in range(0,len(hd0)):
            if hd0[i].find('level') == -1:
                a0 = asd[0][ idx[i]+1 : idx[i+1] ].strip()
                hd.append(hd0[i] + ' ' + a0)

            else:
                lvs = [l.strip() for l in asd[0][ idx[i]+1 : idx[i+1] ].split('|')]
                [hd.append(hd0[i] + ' ' + l) for l in lvs]
        hd = [h.strip() for h in hd]
        if self.verbose:
            print hd

        self.header.append(hd)

        # to identify if the first element is the Spectrum or not...
        ls = 1 if hd[0] == 'Spectrum' else 0

        # now parse associated data
        asd.pop(0) # first line is a break...
        asd.pop(0)

        nan=float('nan')
        while re.search('-'*172, asd[0]) == None:
            if self.verbose:
                print asd[0]

            l = [ l.strip() for l in re.split('\|', asd[0]) ]

            if l[0+ls] != '' or l[1+ls] != '':

                try:
                    # special parsing for some fields
                    str = l[2+ls]
                    (ri,ric) = (str,'')
                    for i in range(0,len(str)):
                        if not( str[i].isdigit() or str[i]=='(' or str[i] ==')' ):
                            (ri,ric) = (str[:i],str[i:])
                            break

                    EiEk = [re.sub('[^0-9\.]','',x) for x in l[5+ls].split('-')] if l[5+ls] != '' else ['nan', 'nan']
                    
                    gigk = l[12+ls].split('-') if l[12+ls] != '' else [nan, nan]
                    
                    # parse all fields into the dictionary
                    d = {'spec'      : l[0] if ls==1 else self.spec,
                         'wave_obs'  : float( re.sub('[^0-9\.]','',l[0+ls]) ) if l[0+ls] != '' else nan,
                         'wave_ritz' : float( re.sub('[^0-9\.]','',l[1+ls]) ) if l[1+ls] != '' else nan, # non-numerics seen: +
                         'rel_int'   : float( re.sub('[^0-9\.-]','',ri) ) if ri != '' else nan, # non-numerics seen: \(\)
                         'rel_int_com': ric,
                         'Aki'       : float( l[3+ls] ) if l[3+ls] != '' else nan,
                         'Acc'       : l[4+ls],
                         'Ei'        : float( EiEk[0] ),
                         'Ek'        : float( EiEk[1] ), # non-numerics seen: \[\]\?+xk
                         'lower_conf': l[6+ls],
                         'lower_term': l[7+ls],
                         'lower_J'   : l[8+ls],
                         'upper_conf': l[9+ls],
                         'upper_term': l[10+ls],
                         'upper_J'   : l[11+ls],
                         'gi'        : float( gigk[0] ) if gigk[0] != '' else nan,
                         'gk'        : float( gigk[1] ) if gigk[1] != '' else nan,
                         'type'      : l[14+ls],
                         'section'   : isec
                        }
                    d['wave'] = d['wave_ritz'] if isnan(d['wave_obs']) else d['wave_obs']
                    self.lines.append(d)
                except:
                    print 'NIST ASD parser error:\n',asd[0]
                    print '\n\tl=',l
                    print '\n\tEiEk=',EiEk, '\n\tgigk=',gigk

            else:
                pass # empty line

            asd.pop(0)


    def plot(self):
        specs = pl.array( list(set([ l['spec'] for l in self.lines ])) )
        specs.sort()
        self.specs = specs

        pl.figure()
        pl.hold('on')
        pl.grid('on')
        pl.jet()

        lines = [] 
        lines_spec = list(pl.zeros(len(specs)))
        for i in range(0,len(self.lines)):
            ispc = pl.find( specs == self.lines[i]['spec'] )
            self.colr = pl.cm.get_cmap()( float(ispc)/len(specs) )
            wl = self.lines[i]['wave']
            ri = float(self.lines[i]['rel_int'])
            lines.append( pl.plot( [wl, wl], [0., ri if not isnan(ri) else 0.], '.-', color=self.colr )[0] )
            lines_spec[ispc] = lines[-1]
        datacursor(lines,formatter='x={x:8.3f}\ny={y:8.3f}'.format)

        pl.rc('text',usetex=True)
        pl.xlabel('$\lambda ~ [\AA]$')
        pl.ylabel('relative intensity [arb]')
        pl.title('Spectrum for '+self.spec+' from NIST ASD')

        if len(specs) > 1:
            pl.legend( lines_spec,specs )

        pl.show()


