import xlrd
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from functools import partial
amount_after_acrit = 15/60

## Parsing routines ##
def opensheet(excel_file_name):
    wb = xlrd.open_workbook(excel_file_name)
    return wb.sheets()[0]

def nextpatient(sheet, start_line):
    name = sheet.cell(start_line, 0).value
    procedure = sheet.cell(start_line, 1).value
    channel = sheet.cell(start_line, 2).value
    sampletype = sheet.cell(start_line, 3).value
    description = sheet.cell(start_line, 4).value
    date = sheet.cell(start_line, 5).value
    i = start_line + 1
    d = []
    min_threshold = 15/60.0
    max_threshold = 60
    t = 0
    while i < sheet.nrows and sheet.cell(i, 0).value == '':
        t = sheet.cell(i,7).value/60.0
        #v = max([0.125, sheet.cell(i,8).value/2.0])
        d.append((t, sheet.cell(i,8).value / 2))
        i = i + 1
    return (i, {'fullname': name, 'data': d, 'procedure': procedure,
                'channel': channel, 'sampletype': sampletype, 
                'description': description, 'date': date})

def tracingsheetgenerator(sheet, f = lambda x:x, filter_fn = lambda x:True):
    last_row = sheet.nrows
    row = 1
    while row < (last_row - 1):
        (row, p) = nextpatient(sheet, row)
        if filter_fn(p):
            yield f(p)

def coagulate(meta, data):
    return {'fullname': meta[0],
            'procedure': meta[1],
            'channel': meta[2],
            'sampletype': meta[3],
            'description': meta[4],
            'date': meta[5].strip(),
            'data': data}

def tracingtextgenerator(filename ,f = lambda x:x, filter_fn = lambda x:True):
    file = open(filename)
    file.readline()
    meta = file.readline().split('\t')
    data = []
    for line in file:
        if line.strip() != '':
            l = line.split('\t')
            if l[0] != '':
                curr = coagulate(meta, data)
                if filter_fn(curr):
                    yield f(curr)
                meta = l
                data = []
            else:
                data.append((float(l[7]) / 60, float(l[8]) / 2))
    #except:
    #foo = line.split('\t')
     #   print(foo[0], foo[1], foo[2], foo[3])
    #print("whoa: " + line.split('\t'))
## -- ##

## 
def deltaamplitude(data):
    min_threshold = 0
    max_threshold = 60
    a_threshold = -1000
    f_threshold = lambda t,a,da: ((da > 0) and (a > a_threshold) and (t > min_threshold) and (t < max_threshold))
    f_dA = lambda a1,a2: (a2[0], 0.5*(a2[1]+a1[1]), a2[1]-a1[1])
    return [f_dA(d1, d2) for (d1, d2) in zip(data[:-1], data[1:]) if f_threshold(f_dA(d1, d2)[0], f_dA(d1, d2)[1], f_dA(d1, d2)[2])]

def maxdeltaamplitude(p):
    dA = deltaamplitude(p['data'])
    m = None
    if (len(dA) > 0):
        m = max(dA, key=lambda x: x[2])
    #return (m[0][0], m[0][1], m[1][1])
    return m

def maxamplitude(p):
    min_threshold = 15/60.0
    max_threshold = 60
    f_threshold = lambda t: t > min_threshold and t < max_threshold
    return max([(t,d) for (t,d) in p['data'] if f_threshold(t)], key=lambda x: x[1])


def writepatientinfo(writer_file, p):
    m = maxdeltaamplitude(p)
    s = ', '.join([str(int(p['fullname'])), str(p['date']), str(p['sampletype']), str(p['description']), str(m[2]), str(m[1])])
    writer_file.write(s + '\n')

def patientwriter(writer_file):
    return lambda p: writepatientinfo(writer_file, p)

def splittime(p):
    c = 0
    tsplit = 0
    asplit = 0
    split_thresh = 0.2
    for (t, a) in p['data']:
        if a > split_thresh:
            c = c+1
            if c == 1:
                tsplit = t
                asplit = a
            if c > 5:
                return (tsplit, asplit)
        else:
            c = 0
            a = 0
    return p['data'][-1]

def amplitude5and10(p):
    tsplit, asplit = splittime(p)
    try:
        a5 = next(filter(lambda d: d[0] >= tsplit + 5, p['data']))
        a10 = next(filter(lambda d: d[0] >= tsplit + 10, p['data']))
    except:
        a5 = p['data'][-1]
        a10 = p['data'][-1]
    return (tsplit, asplit, a5[1], a10[1])

def salientdetails(p):
    m = maxdeltaamplitude(p)
    maxA = maxamplitude(p)
    area, ds = dAintegration(p)
    (tsplit, asplit, a5, a10) = amplitude5and10(p)
    a,b,minval = fit_best(p)
    out = {'fullname':p['fullname'], 
           'date': p['date'],
           'sampletype': p['sampletype'],
           'description': p['description'],
           'Amax': maxA[1],
           'tmax': maxA[0],
           'tsplit': tsplit,
           'asplit': asplit,
           'a5': a5,
           'a10': a10,
           'dASum': area,
           'aModelParam': a,
           'bModelParam': b}
    if m != None:
        out['dAmax'] = m[2]
        out['Acrit'] = m[1]
        out['tcrit'] = m[0]
    else:
        print("No dA values. Problem with this patient: %i %s %s %s"%(out['fullname'], makeexceldatestring(out['date']), out['sampletype'], out['description']))
    return out

def simpledetails(p):
    m = maxdeltaamplitude(p)
    maxA = maxamplitude(p)
    (tsplit, asplit, a5, a10) = amplitude5and10(p)
    out = {'fullname':p['fullname'], 
           'date': p['date'],
           'sampletype': p['sampletype'],
           'description': p['description'],
           'Amax': maxA[1],
           'tmax': maxA[0],
           'tsplit': tsplit,
           'asplit': asplit,
           'a5': a5,
           'a10': a10}
    if m != None:
        out['dAmax'] = m[2]
        out['Acrit'] = m[1]
        out['tcrit'] = m[0]
    else:
        print("No dA values. Problem with this patient: %s %s %s %s"%(str(out['fullname']), out['date'], out['sampletype'], out['description']))
    return out

## Plotting routines ##
def plottimeseriestoaxes(p, ax):
    m = maxdeltaamplitude(p)
    maxA = maxamplitude(p)
    tsplit, asplit, a5, a10 = amplitude5and10(p)
    #ax.set_title('Patient ' + str(p['fullname']))
    ax.plot([i[0] for i in p['data']], [i[1] for i in p['data']])
    ax.scatter(maxA[0], maxA[1], c='g')
    #d = data_after_split(p)
    d=p['data']
    dAs = deltaamplitude(d)
    tdAmax, Acrit, dAmax = maxdeltaamplitude(p)
    amps = [a for t,a,da in dAs if t <= tdAmax+amount_after_acrit]
    ts = [t for t,a,da in dAs if t <= tdAmax + amount_after_acrit]
    ##ax.plot(ts, amps, c='r', linewidth=3)
    if m != None:
        ax.scatter(m[0], m[1], c='r')
        s = "Acrit (dA = %.2f, A = %.2f)" % (m[2], m[1])
        ax.annotate(s, xy=(m[0],m[1]), xytext=(m[0]+1,m[1]), fontsize=8)
    s = "Amax (A = %.2f)" % (maxA[1])
    ax.annotate(s, xy=maxA, xytext=(maxA[0]+1, maxA[1]), fontsize=8)
    ax.set_xlabel("Time (min)")
    ax.set_ylabel("Amplitude (mm)")
    ax.scatter(tsplit, asplit, c='b')
    ax.annotate("Split (t = %2f)" % (tsplit), xy=(tsplit, asplit), fontsize=8)
    ax.scatter(tsplit+5, a5, c='b')
    ax.annotate("A5 (t = %.2f, A = %.2f)" % (tsplit+5, a5), xy=(tsplit+5, a5), fontsize = 8)
    ax.scatter(tsplit+10, a10, c='b')
    ax.annotate("A10 (t = %.2f, A = %.2f)" % (tsplit+10, a10), xy=(tsplit+10, a10), fontsize=8)

def plotbestfit(p, ax, pastAcrit=numpy.Inf, col='r', text_model=False):
    #a,b,minval = fit_best(p)
    if pastAcrit == numpy.Inf:
        a,b,minval = fit_best(p)
    else:
        a,b,minval = fit_best_uptoAcrit(p, pastAcrit)
    dAs = deltaamplitude(p['data'])
    #ax.plot([a for t,a,da in dAs], [da for t,a,da in dAs])
    maxA = max([a for t,a,da in dAs])
    avals = numpy.linspace(0.01, maxA, 100)
    ax.plot(avals, solve_dA(a,b)(avals), c=col, linewidth=2)
    if text_model:
        ax.annotate("Model: a=%.3f, b=%.3f" % (a, b), xy=(maxA/5, 0))
    
def plotdavatoaxes(p, ax):
    dA = deltaamplitude(p['data'])
    if (dA == None):
        return
    tdAmax, Acrit, dAmax = maxdeltaamplitude(p)
    ax.scatter([a for (t,a,da) in dA], [da for (t,a,da) in dA])
    ax.scatter([a for (t,a,da) in dA if t <= tdAmax+amount_after_acrit], [da for (t,a,da) in dA if t <= tdAmax+amount_after_acrit], c='r')
    area, d = dAintegration(p)
    plotbestfit(p, ax, numpy.Inf, 'g', False)
    plotbestfit(p, ax, 0, 'b', False)
    plotbestfit(p, ax, amount_after_acrit, 'r', True)
    #plotbestfit(p, ax, amount_after_acrit, 'r', True)
    ax.legend(['Fit up to Amax', 'Fit up to Acrit', 'Fit up to Acrit + 15 sec'])
    #ax.fill_between([a for (t,a,da) in d], [da for (t,a,da) in d], 0, facecolor='r', alpha=0.7)
    #foo=(d[-1][1]/2, 0)
    #ax.annotate("Area from split to crit: %.2f"%(area), xy=foo)
#     min_threshold = 15/60.0
#     max_threshold = 60
#     f_threshold = lambda t: t > min_threshold and t < max_threshold
#     A = [a for (t,a) in p['data'][1:] if f_threshold(t)]
#     deltaA = [da for (t,da) in dA if f_threshold(t)]
# #    A = [i[1] for i in p['data'][1:]]
# #    deltaA = [i[1] for i in dA]
    # ax.scatter([a for (a, da) in zip(A, deltaA) if da > 0],
    #            [da for (a, da) in zip(A, deltaA) if da > 0])
    ax.set_xlabel("Amplitude (mm)")
    ax.set_ylabel("dA (mm)")
    #ax.set_title('Patient ' + str(p['fullname']))

def plottimeseriesdata(p):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plottimeseriestoaxes(p, ax)
    plt.show()

def plottimedomain(p):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    try:
        plottimeseriestoaxes(p, ax)
    except:
        fig = None
    return (fig, p['fullname'] + "; " + p['sampletype'] + "; " + p['description'])
    
def plotamplitudedomaindata(p):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plotdavatoaxes(p, ax)
    plt.show()

def superplot(p):
    fig, (a1, a2) = plt.subplots(2, 1)
    title = "Patient ID: %s, date: %s, sample type: %s, description: %s" %(str(p['fullname']), str(p['date']), p['sampletype'], p['description'])
    try:
        plottimeseriestoaxes(p, a1)
        plotdavatoaxes(p, a2)
        return fig, title
    except:
        return None, title
        
def superplot_tk(p):
    fig, (a1, a2) = plt.subplots(2, 1)
    plottimeseriestoaxes(p, a1)
    plotdavatoaxes(p, a2)
    plt.show()

def findallpatiententries(sheet, patientid, f=simpledetails, filter_fn = lambda x:True):
    f_patient = lambda p: p['fullname'] == patientid
    #or int(p['fullname']) == patientid
    g = patientgenerator(sheet, f, lambda p: f_patient(p) and filter_fn(p))
    return [m for m in g]

## PDF generation routines ##
def makePDF(fig_generator, outfilename):
    pp = PdfPages(outfilename)
    for f, title in fig_generator:
        if not f == None:
            f.suptitle(title)
            pp.savefig(f)
        else:
            print("problem with figure. " + title)
            
    pp.close()
    plt.close('all')

def makePDFpage(pdfcontext, patientid, sheet, filter_fn = lambda x:True):
    patients = findallpatiententries(sheet, patientid, lambda p:p, filter_fn)
    patients = sorted(patients, key=lambda x: x['date'])
    for p in patients:
        fig, (ax1, ax2)  = plt.subplots(2,1)
        title = "Patient ID: %s, date: %s, sample type: %s, description: %s" %(str(p['fullname']), makeexceldatestring(p['date']), p['sampletype'], p['description'])
        fig.suptitle(title)
        plottimeseriestoaxes(p, ax1)
        plotdavatoaxes(p, ax2)
        ppdfcontext.savefig(fig)
        
def makefullPDF(pdffilename, sheet, filter_fn = lambda x:True):
    g = patientgenerator(sheet, simpledetails, filter_fn)
    ids = set([i['fullname'] for i in g if i != None])
    pp = PdfPages(pdffilename)
    for i in ids:
        makePDFpage(pp, i, sheet, filter_fn)
    pp.close()
    plt.close('all')

## -- ##
    
def groupbypatient(patients):
    d = {}
    for p in patients:
        if p['fullname'] in d:
            d[p['fullname']].append(p)
        else:
            d[p['fullname']] = [p]
    for p in d:
        d[p] = sorted(d[p], key=lambda x:x['date'])
    return d

def makeexceldatestring(datethingie):
    (y, m, d, h, minute, s)=xlrd.xldate.xldate_as_tuple(datethingie, 0)
    s = "%d/%d/%d %d:%02d:%02d" % (m, d, y, h, minute, s)
    #s = str(m) + "/" + str(d) + "/" + str(y) + " " + str(h) + ":" + str(minute)
    return s

def writeoutgroupedpatients(grouped_patients, outfilename):
    names = grouped_patients.keys()
    names = sorted(names)
    outfile = open(outfilename, 'w')
    outfile.write("patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), A5 (mm), A10 (mm), dA integral from split to crit, fit model param a, fit model param b\n")
    for n in names:
        if isinstance(n, (int, float, complex)):
            n = int(n)
        for e in grouped_patients[n]:
            s = ', '.join([str(n), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']), str(e['a5']), str(e['a10']), str(e['dASum']), str(e['aModelParam']), str(e['bModelParam'])])
            outfile.write(s + "\n")

    outfile.close()

def writeoutgroupedpatientsforsampletype(grouped_patients, outfilename, sampletype):
    names = grouped_patients.keys()
    names = sorted(names)
    
    outfile = open(outfilename, 'w')
    outfile.write("patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), A5 (mm), A10 (mm), dA integral from split to crit\n")

    for n in names:
        if isinstance(n, (int, float, complex)):
            n = int(n)
        for e in grouped_patients[n]:
            if e['sampletype'] == sampletype:
                try:
                    s = ', '.join([str(n), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']), str(e['a5']), str(e['a10']), str(e['dASum'])])
                    outfile.write(s + "\n")
                except:
                    print("Problem with patient: %s, %s, %s, %s" % (str(e['fullname']), makeexceldatestring(e['date']), e['sampletype'], e['description']))
    outfile.close()


def groupbysampletype(patients):
    d = {}
    for p in patients:
        if p['sampletype'] in d:
            d[p['sampletype']].append(p)
        else:
            d[p['sampletype']] = [p]
    return d


def timeintegral(p, fun = lambda x: True):
    data = filter(fun, p['data'])
    t0,a0 = next(data)
    acc = 0
    for t1,a1 in data:
        tdelta = t1 - t0
        a = 0.5 * (a0 + a1)
        acc = acc + a * tdelta
        t0,a0 = t1,a1
    return acc

def ly30(p):
    t, amax = maxamplitude(p)
    f = lambda x: x[0] >= t and x[0] <= t + 30
    actual = timeintegral(p, f)
    extrapolated = amax * 30
    return (extrapolated - actual) / extrapolated

def indexanddelta(inhib, other):
    a = other
    b = inhib
    tasplit,aasplit = splittime(a)
    tbsplit,absplit = splittime(b)
    ta = [x-tasplit for x,y in a['data'] if x >= tasplit]
    tb = [x-tbsplit for x,y in b['data'] if x >= tbsplit]
    ya = [y for x,y in a['data'] if x >= tasplit]
    yb = [y for x,y in b['data'] if x >= tbsplit]
    tamax,aamax = maxamplitude(a)
    tbmax,abmax = maxamplitude(b)
    adelt = deltaamplitude(a['data'])
    foo = [(x,y,z) for x,y,z in adelt if x>=tbmax][:4]
    t = tbmax - tbsplit
    yval = [y for x,y in a['data'] if x >= tbmax][0]
    da = (foo[-1][1] - foo[0][1])/(foo[-1][0]-foo[0][0])
    return da

def findmatchie(p, searchdatabase, findfn):
    try:
        gen = (i for i in searchdatabase() if findfn(p, i))
        match = next(gen)
    except:
        print("Could not find match for " + p['fullname'])
        match = None
    return (p, match)

def findmatchie2(p, sd, findfn, sd2, findfn2):
    try:
        gen = (i for i in sd() if findfn(p, i))
        match = next(gen)
    except:
        print("Could not find match for " + p['fullname'])
        match = None
    try:
        gen = (i for i in sd2() if findfn2(p, i))
        match2 = next(gen)
    except:
        print("Could not find match for " + p['fullname'])
        match2 = None
    return (p, match, match2)

def findmatch(generatormanufacturer, findfn):
    return lambda p: findmatchie(p, generatormanufacturer, findfn)

def findmatch2(gm, findfn, gm2, findfn2):
    return lambda p: findmatchie2(p, gm, findfn, gm2, findfn2)

def makehealthymatchedstring(patients):
    metastring = ", ".join(["ID", "CRT Description", "CRT Amax (mm)", "CRT t at Amax (min)", "CRT LY30 (mm^2)", "CFF Description", "CFF Amax (mm)", "CFF t at Amax (min)", "CFF LY30 (mm^2)", "CFTX Description", "CFTX Amax (mm)", "CFTX t at Amax (min)", "CFTX LY30 (mm^2)", "CFTX dA @ CFF Amax (mm)"])
    ck = patients[0]
    ckd = simpledetails(ck)
    s = ckd['fullname'].replace(","," ") + ', ' + ckd['description'].replace(","," ") + ', '
    s = s + str(ckd['Amax']) + ', ' + str(ckd['tmax']) + ', '
    s = s + str(ly30(ck)) + ', '
    cff = patients[1]
    if cff is not None:
        cffd = simpledetails(cff)
        s = s + cffd['description'].replace(",", " ") + ", "
        s = s + str(cffd['Amax']) + ', ' + str(cffd['tmax']) + ', '
        s = s + str(ly30(cff)) + ', '
    else:
        s = s + 'N/A, N/A, N/A, N/A, '
    inhib = patients[2]
    if inhib is not None:
        inhibd = simpledetails(inhib)
        s = s + inhibd['description'].replace(",", " ") + ", "
        s = s + str(inhibd['Amax']) + ', ' + str(inhibd['tmax']) + ', '
        s = s + str(ly30(inhib)) + ", "
    else:
        s = s + 'N/A, N/A, N/A, N/A, '
    if inhib is not None and cff is not None:
        try:
            s = s + str(indexanddelta(cff, inhib))
        except:
            s = s + 'N/A'
            print("Problem indexing and finding dA.")
        
    return (metastring, s)

def healthygetmatched():
    findfn = lambda a,b:a['fullname'].lower()==b['fullname'].lower() and a['description'][0:2].lower()==b['description'][0:2].lower()
    gen_inhib = lambda : tracingtextgenerator("./data/inhib.crd")
    gen_cff = lambda : tracingtextgenerator("./data/cff.crd")
    ##g = tracingtextgenerator("./data/ck.crd", findmatch2(gen_cff, findfn, gen_inhib, findfn))
    fn = lambda p: makehealthymatchedstring(findmatch2(gen_cff, findfn, gen_inhib, findfn)(p))
    g = tracingtextgenerator("./data/ck.crd", fn)
    write_function_csv("whoaglop.csv", g)

def tapgetmatched():
    cff_findfn = lambda a,b:a['fullname'].lower()==b['fullname'].lower() and b['sampletype'][:3]=="CFF" and b['description'].lower().find("txa")<0
    inhib_findfn = lambda a,b:a['fullname'].lower()==b['fullname'].lower() and b['sampletype'][:3]=="CFF" and b['description'].lower().find("txa")>=0
    gen_inhib = lambda : tracingtextgenerator("./data/Combat/tap.crd")
    gen_cff = lambda : tracingtextgenerator("./data/Combat/tap.crd")
    ##g = tracingtextgenerator("./data/ck.crd", findmatch2(gen_cff, findfn, gen_inhib, findfn))
    fn = lambda p: makehealthymatchedstring(findmatch2(gen_cff, cff_findfn, gen_inhib, inhib_findfn)(p))
    g = tracingtextgenerator("./data/Combat/tap.crd", fn, lambda p: p['sampletype'][:3]=="CRT")
    write_function_csv("whoa.csv", g)
    #return next(g)

def combatgetmatched():
    cff_findfn = lambda a,b:a['fullname'][:5].lower()==b['fullname'][:5].lower() and b['sampletype'][:3]=="CFF" and b['description'].find(a['description'])>=0 and b['description'].lower().find("txa")<0
    inhib_findfn = lambda a,b:a['fullname'][:5].lower()==b['fullname'][:5].lower() and b['sampletype'][:3]=="CFF" and b['description'].find(a['description'])>=0 and b['description'].lower().find("txa")>=0

    gen_inhib = lambda : tracingtextgenerator("./data/Combat/combat.crd")
    gen_cff = lambda : tracingtextgenerator("./data/Combat/combat.crd")
    ##g = tracingtextgenerator("./data/ck.crd", findmatch2(gen_cff, findfn, gen_inhib, findfn))
    fn = lambda p: makehealthymatchedstring(findmatch2(gen_cff, cff_findfn, gen_inhib, inhib_findfn)(p))
    g = tracingtextgenerator("./data/Combat/combat.crd", fn, lambda p: p['sampletype'][:3]=="CRT")
    write_function_csv("whoa.csv", g)
    #return next(g)
    
def dAintegration(p, a_right = 100):
    tsplit, asplit, a5, a10 = amplitude5and10(p)
    #tdAmax, Acrit, dAmax = maxdeltaamplitude(p)
    ds = deltaamplitude(p['data'])
    f_threshold = lambda d: d[0] >= tsplit and d[1] < a_right#and d[0] <= tdAmax
    df = [i for i in filter(f_threshold, ds)]
    s = sorted(df, key=lambda x: x[1])
    if len(s) < 2:
        print("danger: ", p['sampletype'])
    cum = 0
    for i in range(1, len(s)):
        pt, pA, pdA = s[i-1]
        ct, cA, cdA = s[i]
        dx = cA - pA
        area = dx * 0.5 * (pdA + cdA)
        cum = cum + area
    return (cum, s)

def solve_dA(a, b):
    f = lambda x: [(a * xi**2) / (numpy.exp(b * xi)-1) for xi in x]
    return f

def squared_error(amplitudes, dAs, model, a_at_max_dA = numpy.Inf):
    cum = 0
    mvals = model(amplitudes)
    amax = max(amplitudes)
    #adAmax, daMax = max(zip(amplitudes, dAs), key=lambda x: x[1])
    n = 0
    for a, da, mv in zip(amplitudes, dAs, mvals):
        #if a > 2*amax / 3:
        #    cum = cum + 0.5*(mv - da)**2
        #else:
        #    cum = cum + (mv - da)**2
        if a <= a_at_max_dA:
            n = n + 1
            cum = cum + (mv - da)**2
    return cum/n
    
def corr_coef(amplitudes, dAs, model):
    m = numpy.mean(dAs)
    sst = 0
    sse = 0
    for da,damodel in zip(dAs, model(amplitudes)):
        sst = sst + (da - m)**2
        sse = sse + (da - damodel)**2
    return 1 - (sse / sst)

def model_gradient(a, b, modelerror, amplitudes, dAs, Acrit=numpy.Inf):
    step = 0.00005
    modela = solve_dA(a + step, b)
    erra = squared_error(amplitudes, dAs, modela, Acrit)
    modelb = solve_dA(a, b + step)
    errb = squared_error(amplitudes, dAs, modelb, Acrit)
    #print((erra - modelerror) / step, (errb - modelerror) / step)
    return ((erra - modelerror) / step, (errb - modelerror) / step)

def find_best_init_val(amplitudes, dAs, a_at_damax = numpy.Inf):
    avals = list(numpy.linspace(0.001, 0.5, 50)) + [0.75, 1.0]
    bvals = list(numpy.linspace(0.001, 0.5, 30)) + [0.75, 1.0]
    amin, bmin = (0,0)
    valmin = numpy.Inf
    for a in avals:
        for b in bvals:
            model = solve_dA(a, b)
            e = squared_error(amplitudes, dAs, model, a_at_damax)
            if e < valmin:
                valmin = e
                amin, bmin = a, b
    return (amin, bmin, e)

def fit_best_dAvA_curve(amplitudes, dAs, a_at_damax=numpy.Inf):
    aold, bold, e = find_best_init_val(amplitudes, dAs, a_at_damax)
    print(aold, bold, e)
    modelold = solve_dA(aold, bold)
    errold = squared_error(amplitudes, dAs, modelold, a_at_damax)
    precision = 0.002
    da, db = model_gradient(aold, bold, errold, amplitudes, dAs, a_at_damax)
    delta = numpy.sqrt(da**2 + db**2)
    i = 0
    step = 0.0005
    while delta > precision and i < 12000:
        aold = aold - step * da
        bold = bold - step * db
        modelold = solve_dA(aold, bold)
        errold = squared_error(amplitudes, dAs, modelold, a_at_damax)
        da, db = model_gradient(aold, bold, errold, amplitudes, dAs, a_at_damax)
        delta = numpy.sqrt(da**2 + db**2)
        i = i + 1
        
    print(aold, bold, errold, i)
    return aold, bold, errold
    
def data_after_split(p):
    tsplit, asplit, a5, a10 = amplitude5and10(p)
    return [(t, a) for t,a in p['data'] if t >= tsplit]

def fit_best(p):
    d = p['data']
    dAs = deltaamplitude(d)
    return fit_best_dAvA_curve([a for t,a,da in dAs], [da for t,a,da in dAs])

def fit_best_uptoAcrit(p, amount_after_Acrit = 0):
    d = p['data']
    dAs = deltaamplitude(d)
    tdAmax, Acrit, dAmax = maxdeltaamplitude(p)
    amps = [a for t,a,da in dAs if t <= tdAmax + amount_after_Acrit]
    davals = [da for t,a,da in dAs if t <= tdAmax + amount_after_Acrit]
    return fit_best_dAvA_curve(amps, davals)

def plot_best_fit(p):
    a,b,minval = fit_best(p)
    dAs = deltaamplitude(p['data'])
    plt.plot([a for t,a,da in dAs], [da for t,a,da in dAs])
    avals = numpy.linspace(1, 60, 100)
    plt.plot(avals, solve_dA(a,b)(avals))
    plt.show()

def calcExpectedArea(a, b):
    return (2*a*1.20206)/b**3

def numericIntegrationOfModel(model, x1, x2):
    n = 1000
    delta = (x2 - x1) / n
    x = numpy.linspace(x1, x2, n)
    cum = 0
    y = model(x)
    for i in range(1, n):
        y1 = y[i-1]
        y2 = y[i]
        cum += delta * 0.5 * (y1 + y2)
    return cum

def grabAreaStuff(p):
    a,b,err = fit_best_uptoAcrit(p, 0)
    #d = data_after_split(p)
    d=p['data']
    dAs = deltaamplitude(d)
    r2 = corr_coef([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a,b))
    return ((a,b), dAintegration(p)[0], calcExpectedArea(a,b), numericIntegrationOfModel(solve_dA(a,b), 0.000001, 100), r2)
    
# "patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, numerical integration to 100mm, analytical integral of model"

def area_details(p):
    e = simpledetails(p)
    a1, b1, e1 = fit_best(p)
    a2, b2, e2 = fit_best_uptoAcrit(p)
    a3, b3, e3 = fit_best_uptoAcrit(p, amount_after_acrit)
    #d = data_after_split(p)
    d=p['data']
    dAs = deltaamplitude(d)
    #r21 = corr_coef([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a1, b1))
    rmse1 = numpy.sqrt(squared_error([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a1, b1), numpy.Inf))
    rmse2 = numpy.sqrt(squared_error([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a2, b2), numpy.Inf))
    rmse3 = numpy.sqrt(squared_error([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a3, b3), numpy.Inf))
    #r22 = corr_coef([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a2, b2))
    area = dAintegration(p)[0]
    n1 = numericIntegrationOfModel(solve_dA(a1, b1), 0.000001, 100)
    n2 = numericIntegrationOfModel(solve_dA(a2, b2), 0.000001, 100)
    n3 = numericIntegrationOfModel(solve_dA(a3, b3), 0.000001, 100)
    c1 = calcExpectedArea(a1, b1)
    c2 = calcExpectedArea(a2, b2)
    c3 = calcExpectedArea(a3, b3)
    return (area, (a1, b1, rmse1, n1, c1), (a2, b2, rmse2, n2, c2), (a3, b3, rmse3, n3, c3))
    #s = ', '.join([str(n), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']))])

def stringify_area_stuffs(p):
    area, full, partial, partial2 = area_details(p)
    e = simpledetails(p)
    name = e['fullname']
    if isinstance(name, (int, float, complex)):
        name = int(name)
        
    s = ', '.join([str(name), str(e['date']), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']), str(area), str(full[0]), str(full[1]), str(full[2]), str(full[3]), str(full[4]), str(partial[0]), str(partial[1]), str(partial[2]), str(partial[3]), str(partial[4]), str(partial2[0]), str(partial2[1]), str(partial2[2]), str(partial2[3]), str(partial2[4])])
    meta = "patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, RMSE, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, RMSE, numerical integration to 100mm, analytical integral of model, a fit to Acrit+15sec, b fit to Acrit+15sec, RMSE, numerical integration to 100mm, analytical integral of model"
    return meta, s
    
def stringify_area_stuffs_excel(p):
    area, full, partial, partial2 = area_details(p)
    e = simpledetails(p)
    name = e['fullname']
    if isinstance(name, (int, float, complex)):
        name = int(name)
        
    s = ', '.join([str(name), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']), str(area), str(full[0]), str(full[1]), str(full[2]), str(full[3]), str(full[4]), str(partial[0]), str(partial[1]), str(partial[2]), str(partial[3]), str(partial[4]), str(partial2[0]), str(partial2[1]), str(partial2[2]), str(partial2[3]), str(partial2[4])])
    meta = "patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, RMSE, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, RMSE, numerical integration to 100mm, analytical integral of model, a fit to Acrit+15sec, b fit to Acrit+15sec, RMSE, numerical integration to 100mm, analytical integral of model"
    return meta, s

def write_function_csv(outfilename, patient_generator):
    outfile = open(outfilename, 'w')
    d = next(patient_generator)
    outfile.write(d[0] + "\n")
    fn = lambda x: outfile.write(x[1] + "\n")
    fn(d)
    while True:
        try:
            d = next(patient_generator)
            fn(d)
        except:
            print("Problem")
    #[fn(l) for l in patient_generator]
    outfile.close()
    
def write_function_to_spreadsheet(filename, sheet, patient_fn, sampletype = "all"):
    outfile = open(filename, 'w')
    if sampletype == "all":
        g = patientgenerator(sheet, patient_fn)
    else:
        g = patientgenerator(sheet, patient_fn, lambda p: p['sampletype'] == sampletype)
        
    d = next(g)
    outfile.write(d[0] + "\n")
    fn = lambda x: outfile.write(x[1] + "\n")
    fn(d)
    [fn(l) for l in g]
    
    outfile.close()

def profile_goodness(amounts, p):
    area, d = dAintegration(p, 10000)
    delts = []
    for am in amounts:
        a,b,e = fit_best_uptoAcrit(p, am)
        delts.append((a, b, area, calcExpectedArea(a, b)))
    return delts

def profile_goodness_fit(amounts):
    return partial(profile_goodness, amounts)



def foobar():
    # g1 = groupbypatient(tracingtextgenerator("./data/inhib.crd"))
    # g2 = groupbypatient(tracingtextgenerator("./data/cff.crd"))
    # for i in g1:
    #     try:
    #         a = g1[i]
    #         b = g2[i]
    #         c = {}
    #         d = {}
    #         for j in a:
    #             c[j['description'][:2]] = a[j]
    #         for j in b:
    #             d[j['description'][:2]] = b[j]
    #         print(str(len(c)))
    #         for x in c:
    #             print(d[x]['description'])
    #     except:
    #         print("danger!: " + str(i))
    pp = PdfPages("glop.pdf")
    g1 = tracingtextgenerator("./data/inhib.crd")
    i = 0
    for a in g1:
        #if i == 5:
        #pp.close()
        #return
        i = i + 1
        g2 = tracingtextgenerator("./data/cff.crd")
        for b in g2:
            if a['fullname']==b['fullname'] and a['description'][:2]==b['description'][:2]:
                try:
                    tasplit,aasplit = splittime(a)
                    tbsplit,absplit = splittime(b)
                    ta = [x-tasplit for x,y in a['data'] if x >= tasplit]
                    tb = [x-tbsplit for x,y in b['data'] if x >= tbsplit]
                    ya = [y for x,y in a['data'] if x >= tasplit]
                    yb = [y for x,y in b['data'] if x >= tbsplit]
                    tamax,aamax = maxamplitude(a)
                    tbmax,abmax = maxamplitude(b)
                    adelt = deltaamplitude(a['data'])
                    #plt.plot((t, t), (abmax, yval[0]), 'k')
                    foo = [(x,y,z) for x,y,z in adelt if x>=tbmax][:4]
                    t = tbmax - tbsplit
                    yval = [y for x,y in a['data'] if x >= tbmax][0]
                    da = (foo[-1][1] - foo[0][1])/(foo[-1][0]-foo[0][0])
                    fig = plt.figure()
                    ax = fig.add_subplot(111)
                    fig.suptitle(a['fullname'] + ", " + b['description'])
                    ax.plot(ta, ya, 'g')
                    ax.plot(tb, yb, 'b')
                    ax.legend(( "Inhib", "CFF"), loc=4)
                    ax.set_xlabel('Time (min)')
                    ax.set_ylabel('Amplitude (mm)')
                    ax.scatter(t, yval)
                    ax.annotate("dA/dt=%.2fmm/min"%(da), xy=(t, yval-1))
                    #plt.show()
                    pp.savefig(fig)
                except:
                    print("Bah: " + a['fullname'])
                    
                break
    pp.close()
    plt.close('all')
            #ax.annotate("Area from split to crit: %.2f"%(area), xy=foo)
# def foo(data, outfile, amounts):
#     meta = "patient ID, date-time, sample type, description, empirical dA integration (mm^2), offset from Acrit (min), model a, model b, model dA integration (mm^2), difference between model and empirical area (mm^2)"
#     for i in patientgenerator(sheet, simpledetails, 'CN'):
#         s = ','.join(i['fullname'], makeexceldatetime(i['date']), i['sampletype'], 
    
#     for d in data:
#         s = 
#         for am in amounts:
            
def foobar2():
    # g1 = groupbypatient(tracingtextgenerator("./data/inhib.crd"))
    # g2 = groupbypatient(tracingtextgenerator("./data/cff.crd"))
    # for i in g1:
    #     try:
    #         a = g1[i]
    #         b = g2[i]
    #         c = {}
    #         d = {}
    #         for j in a:
    #             c[j['description'][:2]] = a[j]
    #         for j in b:
    #             d[j['description'][:2]] = b[j]
    #         print(str(len(c)))
    #         for x in c:
    #             print(d[x]['description'])
    #     except:
    #         print("danger!: " + str(i))
    g1 = tracingtextgenerator("./data/inhib.crd")
    i = 0
    outfile = open("glop.csv", 'w')
    desc ="ID, description, t @ Amax (min), Amax (mm), Inhib dA/dt @ Amax (mm/min)"
    print(desc)
    outfile.write(desc + '\n')
    for a in g1:
        #if i == 5:
        #pp.close()
        #return
        i = i + 1
        g2 = tracingtextgenerator("./data/cff.crd")
        for b in g2:
            if a['fullname']==b['fullname'] and a['description'][:2]==b['description'][:2]:
                try:
                    tasplit,aasplit = splittime(a)
                    tbsplit,absplit = splittime(b)
                    ta = [x-tasplit for x,y in a['data'] if x >= tasplit]
                    tb = [x-tbsplit for x,y in b['data'] if x >= tbsplit]
                    ya = [y for x,y in a['data'] if x >= tasplit]
                    yb = [y for x,y in b['data'] if x >= tbsplit]
                    tamax,aamax = maxamplitude(a)
                    tbmax,abmax = maxamplitude(b)
                    adelt = deltaamplitude(a['data'])
                    #plt.plot((t, t), (abmax, yval[0]), 'k')
                    foo = [(x,y,z) for x,y,z in adelt if x>=tbmax][:4]
                    t = tbmax - tbsplit
                    yval = [y for x,y in a['data'] if x >= tbmax][0]
                    da = (foo[-1][1] - foo[0][1])/(foo[-1][0]-foo[0][0])
                    name = b['fullname']
                    tthing = b['description']
                    sfoo = name + ", " + tthing + ", " + str(t) + ", " + str(abmax) +  ", " + str(da)
                    print(sfoo)
                    outfile.write(sfoo + '\n')
                except:
                    print("Bah: " + a['fullname'] + ", " + a['description'])
                break
    outfile.close()
                    
    
    

    

