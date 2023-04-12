#!/usr/bin/env python
""" Scripts to import JPK AFM force data from force spectroscopy files
"""

import os
from zipfile import ZipFile
from numpy import linspace, zeros, frombuffer, asarray
from AFMforce import AFMforce

###############################################################################
########## Actual routines:
class JPKforce(AFMforce):
    """ An envelop class to maintain JPK force curves. This is meant to
        simplify the data loading and handling and provide easy access to
        parameters and values, such as piezo position, force, deflection
        position, time, etc.
    """

    def read_data(self, filename):
        self.dataset = {}

        raw_data = ReadForceData(filename)
        if raw_data != {}:
            cal_data = CalibrateData(raw_data)
        else:
            return
        #convert stuff...
        self.dataset['segments'] = cal_data['segments']
        #position:
        if "force-scan-series" in raw_data:
            a = raw_data['force-scan-series']
            if 'header' in a and 'position' in a['header']:
                a = a['header']['position']
                self.dataset["position"]= {}
                self.dataset['position']["data"] = asarray( [a['x'],a['y']] )
                self.dataset['position']['unit'] = "m"

        for seg in self.dataset['segments']:
            self.dataset[seg] = {}
            currseg = self.dataset[seg]

            c = cal_data[seg]

            if "name" in c:
                currseg['name'] = c['name']
            #segment name?

            if 'vDeflection' in c:
                cs = c['vDeflection']
                #current segment in the dataset needs force:
                currseg['force'] = {}

                nc = currseg['force']

                #if there is no raw, then data has the raw data,
                #and no other things are defined:
                if 'raw' in cs:
                    nc['raw_data'] = cs['raw']
                    nc['raw_unit'] = "V"
                    #it is possible that the sensor response is calibrated,
                    #but the spring constnt is not
                    if "Distance data" in cs:
                        #this case, we have all defined:
                        nc['deflection'] = cs['Distance data']
                        nc['deflection_unit'] = cs['Distance unit']
                        nc['data'] = cs['data']
                        nc['unit'] = cs['unit']
                        if 'force constant' in cs:
                            self.dataset['force constant'] = cs['force constant']
                        if 'Sensor response' in cs:
                            self.dataset['sensor response'] = cs['Sensor response']

                    else:
                        #we do not have force:
                        nc['deflection'] = cs['data']
                        nc['deflection_unit'] = cs['unit']
                        if 'force constant' in cs:
                            self.dataset['force constant'] = cs['force constant']
                        if 'Sensor response' in cs:
                            seft.dataset['sensor response'] = cs['Sensor response']
                else:
                    #nothing is defined but raw data:
                    nc['raw_data'] = cs['data']
                    nc['raw_unit'] = cs['unit']
            #end taking over vDeflection
            for h in ['strainGaugeHeight','capacitiveSensorHeight',\
                                'measuredHeight', 'height']:
                if h in c:
                    cs = c[h]

                    #current segment in dataset needs height:
                    currseg['height'] = {}
                    nc = currseg['height']
                    nc['data'] = cs['data']
                    nc['unit'] = cs['unit']
                    if h in ['capacitiveSensorHeight','measuredHeight']:
                        nc['type'] = "measured"
                    else:
                        nc['type'] = "piezo"

                    break
            #end this for here

            if 'time' in c:
                #current segment in the dataset needs time:
                currseg['time'] = {}
                currseg['time']['data'] = c['time']['data']
            else:
                cs = raw_data['segment-headers'][seg]
                if 'force-segment-header' in cs:
                    cs = cs['force-segment-header']
                    if 'duration' in cs:
                        t = cs['duration']
                        N = int( cs['num-points']) if 'num-points' in cs else len(self.dataset['height']['data'])
                        currseg['time'] = {}

                        if t > 0:
                            currseg['time']['data'] = linspace(0, t, N)
                        else:
                            currseg['time']['data'] = zeros(N)
            #end if time
        #import is complete
        self.filename = filename
        if 'date' in raw_data:
            self.date = raw_data['date']
    #end read_data

#end JPKforce

###############################################################################
########## JPK specific functions
###############################################################################

def header_to_dict(txt):
    """ Convert the header text structures to a dict structure
        txt: the text content
        """
    #why do I get here a bytes object, I can not tell... but clean the things then:
    if type(txt).__name__ == 'bytes':
        txt = txt.decode('utf8')

    if '\n' in txt:
        txtlst = txt.split('\n')
    else:
        txtlst = txt

    res = {}
    for t in txtlst:
        t = t.strip()
        if len(t) < 1:
            continue;

        if t[0:2] == '###':
            #remark line
            if 'comment' in res:
                if type(res['comment']).__name__ == 'list':
                    res['comment'].append(t)
                else:
                    res['comment'] = [res['comment'],t]
            else:
                res['comment'] = t
            #print("%s" %t)
        elif t[0] == '#':
            res['date'] = t[1:]

        #these description fileas are key = value pairs if not comments
        elif '=' in t:
            key, tval = t.split('=',1)
            keylist = key.split('.')
            #clean off accidental spaces (should not be any)
            tval=tval.strip()

            try:
                #tval.isdigit does not help for engineering format
                #like 5E-7...
                #simple and dirty:
                val = float(tval)
            except:
                #not a number, it may be bool text:
                if tval.lower() == "true":
                    val = True
                elif tval.lower() == "false":
                    val = False
                #then consider it as a text:
                else:
                    val = tval
            #now we have val set
            #print("value:", val)
            #construct res with embedded dict set according to the
            #names:
            c = res
            for k in keylist[:-1]:
                if not k in c:
                    c[k] = {}
                c = c[k]
            #end for
            c[keylist[-1]] = val

        else:
            print("Invalid line: %s" %t)
    #end for
    return res
#end of text_to_dict


def ReadForceData(filename):
    """ Read the descriptors of the JPK force measurement files:
        The file contains one general desctiption as: header.properties
        and special ones for each segment as: segment-header.properties

        We convert these headers to a dict structure, containing
        all parameters from the general header and the various
        segments. Check for the keys 'segment-0'...
        Each segment contains all the segment headers. A 'channels'
        key contains a list of recorded channels.
        Each channel name is a key then, containing further data
        of that channel.

        in the newest system there is a shared folder containing
        information about the segments...

        The function constructs a dict system from the headers,
        containing all kinds of information.
        Important:
        the old system contained all calibration related information
        in the ['segment-headers'][segment]...
        while the new contains the calibration in the 'shared-headers'.
        The key for the shared-headers is in the lcd-info.

        returns
        a dict with all headers and descriptors.
    """

    if not os.path.isfile(filename):
        print("File not found error")
        return {}
    else:
        fp = ZipFile(filename)

    #in the zip file, one can reach the individual files in the
    #filelist, so we collect them to a list:
    #filenamelist = map(lambda x: x.filename, fp.filelist)
    # in python 3, use a for...
    filenamelist = [x.filename for x in fp.filelist]

    #header.properties describes the experiment:
    if not "header.properties" in filenamelist:
        print("Unknown file format!")
        print("File structure contains:")
        for i in filenamelist:
            print("%s" %i)
        return {}
    #end if header exists

    #header.properties is a text file describing the general
    #header, and available in the root of the file
    #the file is a zipped folder containing:
    # header.properties
    # segments/
    #the main result is a dict, starting with this header:
    res = header_to_dict(fp.read("header.properties"))
    #so far this was general, but the details are different from now:
    #newer (ver 2.0) versions have a shared header as well:
    #tmp = filter(lambda x: "shared-data" in x, filenamelist)
    tmp = [x for x in filenamelist if "shared-data" in x]

    if tmp != []:
        res['shared-header'] = header_to_dict(fp.read(tmp[0]))
    #end if we have shared-data

    #figure out the segment number:
    #in the old one, we may have got 4 for segments even when 2 were
    #recorded, so look for the files;
    #and ignore the new version for the moment
    #in python3 this is an iterator object:
    #tmp = filter(lambda x: "segment-header.properties" in x, filenamelist)
    tmp = [x for x in filenamelist if "segment-header.properties" in x]

    N = len(tmp)
    res['N-segment'] = N
    res['segment-headers'] = {}
    res['segments'] = {}

    #for i in range(len(res['segment-headers'])):
    for i in range(N):
        #old and new versions have this:
        #Pull up the segment headers, and store them both temporary and permanent:
        actual_header = res['segment-headers'][i]= header_to_dict(fp.read(tmp[i]))
        chlist = actual_header['channels']['list'].split()
        channels = {}
        #point to the actual position between segments:
        if 'channel' in actual_header:
            actual_channels = actual_header['channel']
        else:
            #not found, it may not even be recorded (old files)
            #so no trouble, just skip it, or do some debugging...
            # print("actual segment", i, "has no key: channel")
            # print("keys are:")
            # print(actual_header.keys())
            continue

        for j in chlist:
            #we need a data type:
            #if it is new, we have lcd-info to deal with:
            if res['file-format-version'] >= 2.0:
                #new type, we have the channel info in the shared, 1x
                #for all segments. For this we need to find out the "lcd":
                lcd = '%d' %(int(actual_channels[j]['lcd-info']['*']))
                ftype = res['shared-header']['lcd-info'][lcd]['type']
            else:
                #old way:
                ftype = actual_channels[j]['data']['type']
            #end if version
            #where is the file: reconstruct the file name
            datafile = "segments/%d/%s" %(i, \
                                actual_channels[j]['data']['file']['name'])

            #now, what we do depends on the ftype in the headers:
            if ftype == "float-data":
                datakey = ">f4"
            elif ftype == "integer-data":
                #datakey = ">i32" changed in numpy 1.10
                datakey = ">i4"
            elif ftype == "short":
                datakey = ">h"
            else:
                #not probable, but who knows?!
                raise ValueError("Unknown data type: %s" %ftype)
            #end if ftype

            #get the file:
            buff = fp.read(datafile)
            # data = fromstring(buff, dtype=datakey)
            data = frombuffer(buff, dtype=datakey)
            channels[j] = data
        #end for chlist
        res['segments'][i] = channels
    #end for segment_headers

    return res
#end of ReadForceData

def CalibrateData(jpkforcedata):
    """
        Convert JPK data to a calibrated data set based on the
        saved parameters. It does not contain any analysis except
        extracting parameters from the headers and doing the conversions
        described there.
        The result is equivalent data to the .out file of the JPK software.

        Parameters:
        a data set read with the ReadForceData

        Return:
        a dict containing the calibrated data

        - each channel is converted in each segment
        - add the end unit to each segment
        - each sensor parameter has various values:
            'data':                 the converted data
            'unit':                 unit of the converted data
            'raw':                  raw data for possible reprocessing
            'raw unit'              units of the raw data (V for Volt)
            For deflection data these are also saved:
            'Force constant':       spring constant if available
            'Force constant unit':  the unit of the spring constant
            'Distance unit':        for deflections
            'Distance':             the distance data after sensor/distance 
                                    conversion
    """

    if type(jpkforcedata).__name__ != "dict":
        raise ValueError("Expected a dict structure, see ReadForceData")
    #end if
    if 'segments' not in jpkforcedata:
        print("No measured data is found")
        return {}

    res = {}

    for seg in jpkforcedata['segments'].keys():
        #print("processing segment: %s" %segment)
        segment_headers = jpkforcedata['segment-headers']
        segment = jpkforcedata['segments'][seg]
        #the information structure for the segment/channels:
        #for the newest system, it is still not enough, we need
        #the shared headers as well
        actual_channels = segment_headers[seg]['channel']
        force_segment_header = segment_headers[seg]['force-segment-header']

        resseg = {}

        #the keys of each segment are the channels:
        for dsetkey in segment.keys():
            #print("processing segment: %s" %dsetkey)
            x = segment[dsetkey]
            #we put the results into a subdict:
            resseg[dsetkey] = {}

            if len(x) < 1:
                print("Warning: Empty data set in %s" %dsetkey)
                continue
            #end if

            #we need to find the calibration info. This will be a subdict
            #held in the convset variable.
            #if it is new, we have lcd-info to deal with:
            #old formats: 0.12, 0.5
            #new format: 2.0 from Sept 2013 on
            if jpkforcedata['file-format-version'] >= 2.0:
                #new type, we have the channel info in the shared, 1x
                #for all segments. For this we need to find out the "lcd":
                lcd = '%d' %(int(actual_channels[dsetkey]['lcd-info']['*']))
                convset = jpkforcedata['shared-header']\
                        ['lcd-info'][lcd]
                encoderset = convset
            else:
                #old way:
                convset = actual_channels[dsetkey]
                encoderset = convset['data']

            if 'encoder' in encoderset :
                encoderset = encoderset['encoder']['scaling']
                a = encoderset['multiplier']
                b = encoderset['offset']
                unit = encoderset['unit']['unit']

                if encoderset['type'] == 'linear' and\
                    encoderset['style'] == 'offsetmultiplier':
                    x = a*x + b
                else:
                    #piezoelectric actuators may have polynomial
                    #conversions,... so this is an open design
                    print("Unknown encoder scaling")
            #end if

            #backup the raw data first after encoder corrections:
            resseg[dsetkey]['raw'] = x.copy()
            #this is not the complete location, but we had to check for
            #the encoder conversion, which is not calibration, but still
            #a scaling
            #now go for the conversion:
            convset = convset['conversion-set']

            #calibrated? How?
            convkeys = convset['conversions']['list'].split()


            #default calibration is:
            for calit in convkeys:
                calib = convset['conversion'][calit]

                #if defined is false, we can not do anything
                if calib['defined'] == False:
                    #print("%d %s %s conversion is not defined" \
                    #%(seg, dsetkey, calit))
                    continue
                #end if defined

                #zoom down:
                #this should be redundant check:
                if  'scaling' not in calib:
                    print("no scaling defined")
                    continue
                #end if
                sc = calib['scaling']

                if 'unit' in sc and 'unit' in sc['unit']:
                    unit = sc['unit']['unit']
                else:
                    unit = 'unknown'
                #end unit

                #we have some more description:
                #calibration type: linear
                #calibration offset and multiplier
                ctype = sc['type']
                if ctype != "linear":
                    raise ValueError("Unknown scaling: %s in %s of %s" \
                                    %(ctype,calit, dsetkey))
                #end if

                #now get the offset and the multiplier

                #multiplier:
                a = sc['multiplier'] if 'multiplier' in sc else 1.0

                ktext = "scaling.offset"
                b = sc['offset'] if 'offset' in sc else 0.0

                scstyle = sc['style']
                if ctype  == "linear" and \
                        scstyle == 'offsetmultiplier':
                    #print("%s is done" %calit)
                    x = a*x + b
                else:
                    print("%s is invalid" %calit)
                #end if

                #record the force constant:
                #print("calit is", calit)
                #print('keys', resseg[dsetkey].keys())
                if calit == 'distance':
                    resseg[dsetkey]['Distance data']= x
                    resseg[dsetkey]['Distance unit'] = unit
                    resseg[dsetkey]['Sensor response'] = a
                    resseg[dsetkey]['Sensor response unit'] = "%s/V" %\
                            unit

                #if we have a calibration constant force, it may be
                #invalid. If the cantilever is not calibrated, we have
                #only a raw Voltage as output
                #Such still needs conversion, but there is no force constant!
                elif calit == 'force' and 'Distance unit' in resseg[dsetkey]:
                    #print('a:', a, 'b:',b)
                    #print('keys', resseg[dsetkey].keys())
                    resseg[dsetkey]['force constant'] = a
                    resseg[dsetkey]['force constant unit'] = "%s/%s" %\
                            (unit, resseg[dsetkey]['Distance unit'])
                #end if distance or force

            #end for calibration
            resseg[dsetkey]['data'] = x
            resseg[dsetkey]['unit'] = unit

        #end for dset in segment

        #what kind of segment is this?
        #standards are extend and retract, but we may have
        #various pause versions as well
        #let us handle the first two:
        if 'name' in force_segment_header:
            namedict = force_segment_header['name']
            if 'base-object-name' in namedict:
                namedict = namedict['base-object-name']

            if 'name' in namedict:
                name = namedict['name']
        else:
            name = 'unknown'
        #end if
        if 'extend' in name:
            name = "approach"
        elif 'retract' in name:
            name = "retract"
        #end if

        resseg['name'] = name
        res[seg] = resseg
    #end for segments
    res['segments'] = [x for x in res.keys()]
    return res
#end CalibrateData

