#!/usr/bin/env python
""" Scripts to import JPK AFM force data from force spectroscopy files
"""

import os
from numpy import linspace, zeros, fromstring, asarray
from AFMforce import AFMforce

###############################################################################
########## Actual routines:
class Brukerforce(AFMforce):
    """ An envelop class to maintain Bruker force curves. This is meant to
        simplify the data loading and handling and provide easy access to
        parameters and values, such as piezo position, force, deflection
        position, time, etc.
    """

    def read_data(self, filename):
        self.dataset = {}

        raw_data = bruker_importer(filename)
        if raw_data == {}:
            return
        #convert stuff...
        self.dataset['segments']=[]
        #we have only extend and retract segments:
        segmentnames = ['extend', 'retract']

        for i in range(2):
            if segmentnames[i] in raw_data:
                self.dataset[i] = {}
                self.dataset['segments'].append(i)
                if segmentnames[i] == 'extend':
                    self.dataset[i]['name'] = 'approach'
                else:
                    self.dataset[i]['name'] = segmentnames[i]

                cs = raw_data[segmentnames[i]]
                if 'height' in cs :
                    cs_2 = cs['height']
                    cn = self.dataset[i]
                    cn['height'] = {}
                    cn = cn['height']

                    if 'Distance data' in cs_2:
                        cn['data'] = cs_2['Distance data']
                        cn['unit'] = cs_2['Distance unit']
                        cn['type'] = cs_2['type']
                        #Bruker stores the piezo distance the other way:
                        # 0 is full retracted, max is closest to the surface
                        cn['data'] = cn['data'].max() - cn['data']

                if 'deflection' in cs:
                    if 'Spring Constant' in raw_data:
                        self.dataset['force constant'] = \
                                float( raw_data['Spring Constant']['value'] )
                    if 'Sens. DeflSens' in raw_data and \
                            raw_data['Sens. DeflSens'] != 1.0 :
                        self.dataset['sensor response'] = \
                                float( raw_data['Sens. DeflSens']['value'] )

                    #this is the force part:
                    cs_2 = cs['deflection']
                    cn = self.dataset[i]
                    cn['force'] = {}
                    cn = cn['force']
                    if 'raw' in cs_2:
                        cn['raw_data'] = cs_2['raw']
                        cn['raw_unit'] = cs_2['raw_unit']
                    if 'Distance data' in cs_2:
                        cn['deflection'] = cs_2['Distance data']
                        cn['deflection_unit'] = cs_2['Distance unit']
                    elif 'sensor response' in self.dataset and 'raw' in cn:
                        cn['deflection'] = cn['raw']* self.dataset['sensor response']
                        cn['deflection_unit'] = \
                                raw_data['Sens. DeflSens']['unit'].split('/',1)[0]

                    if 'data' in cs_2:
                        cn['data'] = cs_2['data']
                        cn['unit'] = cs_2['unit']
                    elif 'force constant' in self.dataset and 'deflection' in cn:
                        cn['data'] = cn['deflection']*self.dataset['force constant']
                        cn['unit'] = 'nN' if 'nm' in cn['deflection_unit'] else 'N'

                    if 'unit' in cn:
                        if cn['unit'] == 'pN' and \
                            'height' in self.dataset[i] and\
                            'unit' in self.dataset[i]['height']:
                            if self.dataset[i]['height']['unit'] == 'nm':
                                cn['data'] = cn['data'] / 1000
                                cn['unit'] = 'nN'

                #end collecting defleciton data
                if 'time' in cs:
                    cs_2 = cs['time']
                    cn = self.dataset[i]
                    cn['time'] = {}
                    cn = cn['time']
                    if 'data' in cs_2:
                        cn['data'] = cs_2['data']
                    if 'unit' in cs_2:
                        cn['unit'] = cs_2['unit']

            #end for in segments
        if 'X Offset' in raw_data and 'Y Offset' in raw_data and\
            'Stage X' in raw_data and 'Stage Y' in raw_data:
                # stage values are in micrometers
                # piezo should be in nm
            self.dataset['position'] = {}
            print(raw_data['Stage X'], raw_data['Stage Y'])
            print(raw_data['X Offset'], raw_data['Y Offset'])

            if 'unit' in raw_data['X Offset']:
                self.dataset['position']['unit'] = raw_data['X Offset']['unit']
            else :
                self.dataset['position']['unit'] = 'nm'
            if self.dataset['position']['unit'] == 'nm':
                X = raw_data['Stage X']['value'] * 1000 + raw_data['X Offset']['value']
                Y = raw_data['Stage Y']['value'] * 1000 + raw_data['Y Offset']['value']
            else:
                X = raw_data['Stage X']['value'] + raw_data['X Offset']['value']
                Y = raw_data['Stage Y']['value'] + raw_data['Y Offset']['value']

            self.dataset['position']['data'] = asarray([X, Y])
            #self.dataset['position']['data'] = asarray( [raw_data['X Offset']['value'], \
            #                    raw_data['Y Offset']['value']])

        if 'Date' in raw_data:
            self.dataset['date'] = raw_data['Date']
            self.date = raw_data['Date']['value']
        #import is complete
        self.filename = filename
    #end read_data
#end Brukerforce

###############################################################################
########## Bruker specific functions
###############################################################################
def bruker_importer(filename):
    """ to import afm exported data"""

    if os.path.isfile(filename):
        fp = open(filename, encoding="ascii", errors="ignore")
    else:
        return {}

    txtbuff = fp.readlines()
    fp.close()

    header = {}

    N = len(txtbuff)
    print("Buffer is %d long" %N)

    i = 0
    t =""
    for i in range(N):
        t = txtbuff[i]
        if t[1:3] == "\\?" or t[1:3] == "?*":
            print("found header")
            break
    #end for
    i0 = i+1

    print("start at: %d" %i0);
    for i in range(i0, N):
        t = txtbuff[i]
#        print(t)

        subkey = ""
        defaultkey = ""
        default = ""
        defaultunit = ""
        unit = ""

        if "file list end" in t.lower():
            print("Header ended")
            break;

        elif t[1:3] == "\\?":
            print("found \\?")
            continue

        #header looks like "\....."
        elif t[:2] == '"\\' and t[-2] == '"':
            #readlines keeps the end of line \n character:
            tt = t[2:-2]
#            print("interpreting: ", tt)

            #if the line is "\@..., we have a special type of header
            if tt[0] == '@':
                #first strip it, and see if we have a subkey or not
                #that is "\@text:text:..."
                kv = tt[1:].split(':',2)
                if len(kv) == 2:
                    k, v = kv
                    subkey = ""
                else:
                    k, subkey, v = kv
                #end if
                vlist = v.split(None,1)
                #vlist[0] is 'S' or 'V'
                #vlist[1] is [...] or (...) or value
                #example:
                # "\@Sens. Log(Resistance): V 1.000000 "
                # "\@Sens. LogSneddonModulusSens: V 1.000000 log(Pa)/log(Arb)"
                # "\@Z center: V [Sens. Zsens] (0.006713867 V/LSB) 55.00000 V"
                # "\@2:DeflectionLimitLockIn3LSADC1: V (0.0000000057 V/LSB) 24.57600 V"
                #and we have to handle all these
                if len(vlist) == 1:
                    # v = vlist[0].strip() if v != [] else ""
                    v = vlist[0].strip() if v else ""
                else:
                    #o.k., we have more than a value here
                    vv= vlist[0].strip()

                    if vv[0] == 'S' or vv[0] == 'V':
                        vrest = vlist[1].strip()
                        vtype = vv[0]

                        #print(vrest)
                        if vrest[0] == '[':
                            v = vrest[1:].split("]",1)
                            detfaultkey = v[0]
                            vrest = ""

                            if len(v) > 1:
                                vrest = v[1].strip()

                        #end more than []...
                        if vrest[0] == '(':
                            v = vrest[1:].split(')',1)
                            vv = v[0].strip()
                            if vv != "":
                                default, defaultunit = vv.split(None,1)
                            vrest = ""

                            if len(v) > 1:
                                vrest = v[1].strip()
                        #handled default and value

                        if vrest != "":
                            #print("again",vrest)
                            v = vrest.split(None,1)

                            if len(v) == 1:
                                unit = ""
                                v = v[0].strip()
                            else:
                                unit = v[1]
                                v = v[0].strip()
                            #print(v, "unit:", unit)

                    else:
                        v = v[0].strip()


            #now, what if we have a usual line:
            else:
                kv = tt.split(':',1)

                #print(kv)
                if len(kv) < 2:
                    k = kv[0]
                    v = ""
                else:
                    k,resv = kv

                v = resv.strip().split(None, 1)
                #end if value is empty

                #number and unit, or several words:
                if len(v) == 2 :
                    unit = v[1]

                    #we may have a space in strings, but they are
                    #not value unit pairs, just strings:
                    #if v[0].isnumeric(): -> fails on negative numbers!
                    try:
                        v = float(v[0])
                    except:
                        #o.k., it was a string
                        v = resv.strip()
                        unit = ""

                elif len(v) == 1:
                    v = v[0]
                else:
                    v = ""

                #sometimes it is still a number...
                try:
                    v = float(v)
                except:
                    pass;
            #end if \@ or just \...

            #now fill in the values and subvalues under key: k

            if subkey != "":
                if k not in header:
                    header[k] = {}

                #we have to generate the lower dir structure:
                header[k][subkey] = {}

                header[k][subkey]['value'] = v

                if default != "":
                    header[k][subkey]['default'] = default.strip()
                if defaultunit != "":
                    header[k][subkey]['default unit'] = defaultunit.strip()
                if defaultkey != "":
                    header[k][subkey]['reference'] = defaultkey.strip()
                if unit != "":
                    header[k][subkey]['unit'] = unit.strip()
            elif k != "":
                if k not in header:
                    header[k] = {}

                header[k]['value'] = v
                if unit != "":
                    header[k]['unit'] = unit.strip()
            #end if how are the keys

#        else:
#            print("can not interpret %s" %t)
    #end for
    #now comes the header of the data table:
    i = i+1
    #there may be some white space lines...
    while (txtbuff[i].strip() == ""):
        #print(i)
        i = i+1
    #end killing empty lines

    #print("data start at:", i, txtbuff[i])
    datahead = txtbuff[i].split('\t')[:-1]
    i0 = i+1
    data = []

    for i in range(i0, N):
        #the last field is from \n, so drop it!
        dline = txtbuff[i].split('\t')[:-1]
        dlN = len(dline)
        #now we have to pull the data in:
        if i == i0:
            data = [[float(d)] for d in dline]
        elif len(data) == dlN:
            for j in range(len(data)):
                if dline[j] != '' and not dline[j].isspace():
                    data[j].append(float(dline[j]))
        else:
            print('Corrupted line found at', i, '/', N,':',
                    len(data), '/', dlN)
    #end for in data
    for i in range(len(data)):
        data[i] = asarray(data[i])

    #now process the header of the data table, it contains information in the
    #channel_unit_segment format, where segment is Extend or Retract: Ex, Rt
    #This I can assign to segment 0 and 1, or leave it as extend and retract

    segmentlist = ['Rt', 'Ex']
    channelconvert = {'Defl': 'deflection', 'Calc_Ramp': 'height', \
                        'Height_Sensor': 'height'}
    header['segments'] = []

    for i in  range(len(datahead)):
        currhead = datahead[i]
        #print("analyzing header:", currhead)
        currdata = data[i]
        #if the calibration was not provided, the data is not available:
        if "Not Available" in currhead:
            continue

        if '_' in currhead:
            rawchannel, unit, segment = currhead.rsplit('_',2)
        else:
            rawchannel, unit, segment = currhead.rsplit(' ',2)

        if segment not in segmentlist and unit in segmentlist:
            tmp = segment
            segment = unit
            unit = tmp
        #end if unit and segment are mixed up

        if segment == 'Rt':
            segment = 'retract'
            segment_name = 'retract'

        elif segment == 'Ex':
            segment = 'extend'
            segment_name = 'extend'
        else:
            print("unknown segment:", segment)
        #end  sorting out segment
        if segment not in header['segments']:
            header['segments'].append( segment )

        if not segment in header:
            header[segment] = {}
            header[segment]['name'] = segment_name

        if rawchannel in channelconvert:
            channel = channelconvert[rawchannel]
        else:
            channel = rawchannel.lower()
        #cleaned up channel list

        if not channel in header[segment]:
            header[segment][channel] = {}


        #deflection is tricky: we can have
        #deflection in distance, deflection in force and in raw.
        if channel in ['deflection', 'height'] :
            if 'sensor' in rawchannel.lower():
                header[segment][channel]['type'] = 'measured'
            else:
                header[segment][channel]['type'] = 'piezo'
            #print('current unit:', unit)

            if unit == "V":
                #print('filling in raw')

                header[segment][channel]['raw'] = currdata
                header[segment][channel]['raw_unit'] = unit

            elif unit in ["m", "nm", "um"]:
                #print("filling in distance nm")

                header[segment][channel]['Distance data'] = currdata
                header[segment][channel]['Distance unit'] = unit

            elif unit == 'Lsb':
                #print("filling in Lsb")

                header[segment][channel]['LSB data'] = currdata
                header[segment][channel]['LSB unit'] = unit

            else:
                #print("filling in force")

                header[segment][channel]['data'] = currdata
                header[segment][channel]['unit'] = unit
        else:
            header[segment][channel]['data'] = currdata
            header[segment][channel]['unit'] = unit


    #header['data'] = data
    header['segments'].sort()

    return header
#end bruker_importer


