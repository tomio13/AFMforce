#!/usr/bin/env python
""" Scripts to import JPK AFM force data from force spectroscopy files
"""

import os
from numpy import linspace, zeros, fromstring, asarray, polyfit, polyval, abs
from AFMforce import AFMforce
from AFMforce.AFMforce import Smooth

###############################################################################
########## Actual routines:
class Pavoneforce(AFMforce):
    """ An envelop class to maintain Pavone indenter force curves. This is meant to
        simplify the data loading and handling and provide easy access to
        parameters and values, such as piezo position, force, deflection
        position, time, etc.
        It is based on the text dump from the original data by the instrument.
    """

    def read_data(self, filename):
        self.dataset = {}

        raw_data = pavone_importer(filename)
        if raw_data == {}:
            return
        # self.raw_data = raw_data

        #convert stuff...
        self.dataset['segments']=[]
        #we have only extend and retract segments:
        segmentnames = ['extend', 'contact', 'retract']

        self.dataset['position'] = {}
        self.dataset['position']['data'] = asarray([raw_data['X-position'][0],
                                                   raw_data['Y-position'][0]])
        self.dataset['position']['unit'] = raw_data['X-position'][-1]

        if 'segments' in raw_data:
            self.dataset['segments'] = []
        else:
            print('data not found!')
            return

        # transfer segments
        for i, name in enumerate(raw_data['segments']):
            self.dataset['segments'].append(i)
            self.dataset[i]= {}
            self.dataset[i]['name'] = name
            # pick this segment
            ds = raw_data['segments'][name]

            if 'k' in raw_data:
                k = raw_data['k'][0]
                self.dataset['force constant'] = k
            else:
                k = 1.0

            # this we can hardcode, because the deflection is in nm
            # and because the system uses other detection system
            # than a standard AFM
            sensor_response= 1.0
            self.dataset['sensor response']= sensor_response

            # what data are available here?
            if 'Cantilever' in ds:
                self.dataset[i]['force'] = {}

                data = ds['Cantilever']['data']
                self.dataset[i]['force']['deflection'] = data

                self.dataset[i]['force']['deflection_unit'] = ds['Cantilever']['unit']


                if ds['Cantilever']['unit'] == 'm':
                    self.dataset[i]['force']['data'] = k * data * 1E6
                elif ds['Cantilever']['unit'] == 'um':
                    self.dataset[i]['force']['data'] = k * data * 1E3

                else :
                    # it must be 'nm'
                    self.dataset[i]['force']['data'] = k * data

                self.dataset[i]['force']['unit'] = 'nN'

                self.dataset[i]['force']['raw_data'] = data/sensor_response
            # force block done
            if 'Piezo' in ds:
                # turn the data around because Pavone measures from 0 at
                # far away
                data = ds['Piezo']['data']
                self.dataset[i]['height'] = ds['Piezo']
                self.dataset[i]['height']['type'] = 'piezo height'
                self.dataset[i]['height']['data'] = data.max() - data if data.size >0 else data
            # height block done

            if 'Time' in ds:
                self.dataset[i]['time'] = ds['Time']
            # end transfering time

        #end for in segments

            #self.dataset['position']['data'] = asarray( [raw_data['X Offset']['value'], \
            #                    raw_data['Y Offset']['value']])

        if 'Date' in raw_data:
            self.dataset['date'] = f'{raw_data[r"Date"]} {raw_data["Time"]}'
            self.date = self.dataset['date']
        #import is complete
        self.filename = filename
    #end read_data
#end Brukerforce

###############################################################################
########## Pavone specific functions
###############################################################################
def list_to_header(txtlst, stoptext= None):
    """ take an importent list of text lines, and convert them
        to dict entries --> key: value pairs

        Parameters:
            txtlst     a list of text lines
            stoptext    stop if you find these text strings in the line

        Returns:
            a dict containing
            the list index where one should go on
            reading the data header and data,
            and a dict with the extracted key: value pairs
    """
    res = {}
    l = 0
    for thisline in txtlst:
        # cut off the tailing '\n':
        thisline = thisline[:-1]

        if thisline == '':
            l += 1
            continue
        # now, we have a line to be interpreted
        # most header are key\tvalue sets
        if stoptext is not None and all(i in thisline for i in stoptext):
            # print(txtlst[l])
            break
        # end if data table comes

        splitline = thisline.split('\t')
        N_split = len(splitline)

        if N_split % 2 == 0:
            for i in range(int(N_split/2)):
                key = splitline[2*i]
                val = splitline[2*i+1]

                # convert the value to number if possible
                try:
                    val = float(val)
                except ValueError:
                    pass

                # if the key contains a unit in the form (unit)
                # extract it
                if ' (' in key:
                    keylist = key.split(' (', 1)
                    key = keylist[0]
                    # the closing ) is assumed to be the last character
                    val = [val, keylist[1][:-1]]
                # end finding units
                res[key]= val

        # increment counter, we need it to be returned
        l += 1
    # end for in headers

    return (l, res)
# end of list_to_header


def pavone_importer(filename, dN= 200):
    """ Import a text output from Pavone, and return the values
        in a dict.
        At the moment the file does not give information about the
        segments of the force curve, so this function tries
        identifying them.
        Maximal force is most probably the end of the first segment.
        Maximal piezo, or a change from the plateu is tried
        to get the second (contact) and third (retract) segments
        separated.

        Parameters:
            filename    string, path to the file to be read
            dN          integer, length of data segment used
                        to calculate a trend at the force plateau
        Returns:
            a dict containing all loaded information
    """

    if not os.path.exists(filename):
        print('File not found!')
        return {}

    with open(filename, 'rt', encoding='utf8') as fp:
        txtlst = fp.readlines()
    # end pulling up the file

    linenumber, res = list_to_header(txtlst, stoptext=['Time', 'Load'])

    # the last line is with the stop strings, this is the
    # line of the data header:
    dataheader = txtlst[linenumber].split('\t')
    print('data header:', dataheader)

    # interpret the table part
    # take apart the header
    tableunits = {}
    for i, val in enumerate(dataheader):
        if ' (' in val:
            head, unit = val.split(' (', 1)
        else:
            head= val
            unit = ' '

        dataheader[i] = head
        tableunits[head] = unit[:-1]

    # pull the values
    table = [fromstring(i, sep='\t') for i in txtlst[linenumber+1:]]
    table = asarray(table)
    res['data'] = {}
    for i, val in enumerate(dataheader):
        res['data'][val] = table[:,i]

    # now, find the segments:
    # approach, contact and retract
    # ds = abs(diff(res['data']['Piezo']))
    #
    # use the smoothed function to eliminate effects of spikes
    # in the data. Some data sets have them.
    y = Smooth(res['data']['Cantilever'], 30, 10, 'Gauss')
    t = res['data']['Time']

    # approach ends at the (first) maximal cantilever deflection
    i0 = (y == y.max()).nonzero()[0][0]

    # the plateau is pushing the piezo further out or
    # has a slight decay
    # here also use smoothed data
    y = Smooth(res['data']['Piezo'][i0:], 30, 10, 'Gauss')
    t= t[i0:]
    # this is the maximal deflection, it will only
    # decrease. First slower, then fast.
    fit = polyfit(t[:dN], y[:dN], 1)
    Nrest = len(t)
    Nrest = int(Nrest/4) if Nrest > 10*dN else Nrest
    # to find the end of the lpateu
    # we do not need the whole curve, just a part
    t = t[:Nrest]
    dy = abs(y[:Nrest]- polyval(fit, t))
    # if we have a decaying plateu
    if fit[0] < 0:
        i1 = (dy == dy.max()).nonzero()[0][0]
    else:
        #no  decay, maybe increasing even
        i1 = (dy < -3*dy[:dN].std()).nonzero()[0]

        if len(i1) < 1:
            i1 = Nrest
        else:
            i1 = i1[0]

    # i0 is pulled off, so we have to add it back to
    # describe the original data
    indx = asarray([i0, i0+i1,table.shape[0] ])

    res['segments'] = {}
    # first the segments

    segmentnames = ['approach', 'contact', 'retract']
    i0 = 0
    dZ = []
    # print(indx)
    for i, i1 in enumerate(indx):
        # print('setting', i, segmentnames[i])
        name = segmentnames[i]
        res['segments'][name] = {}

        for head in dataheader:
            res['segments'][name][head] = {}
            res['segments'][name][head]['data'] = res['data'][head][i0:i1]
            res['segments'][name][head]['unit'] = tableunits[head]

        dZ.append(abs(res['data']['Cantilever'][i1-1] - res['data']['Cantilever'][i0]))

        # update the start index...
        i0= i1 + 1 if i1 < table.shape[0] else table.shape[0]
    # done adding segments
    print('segment force ranges:', dZ)

    # occasionally the contact part is stored as last, instead of middle
    # int this case try swapping them
    if len(dZ) > 2 and dZ[2] < dZ[1]/10.0 :
        print('swapping segments')
        #tmp = res['segments']['contact'].copy()
        #res['segments']['contact'] = res['segments']['retract'].copy()
        #res['segments']['retract'] = tmp
        tmp = res['segments']['contact']
        res['segments']['contact'] = res['segments']['retract']
        res['segments']['retract'] = tmp

    return res
# end of the importer
