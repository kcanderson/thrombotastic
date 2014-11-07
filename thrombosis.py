import xlrd
import numpy
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

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
        v = max([0.125, sheet.cell(i,8).value/2.0])
        d.append((t, sheet.cell(i,8).value / 2))
        i = i + 1
    return (i, {'fullname': name, 'data': d, 'procedure': procedure,
                'channel': channel, 'sampletype': sampletype, 
                'description': description, 'date': date})

def patientgenerator(sheet, f = lambda x:x, filter_fn = lambda x:True):
    last_row = sheet.nrows
    row = 1
    while row < (last_row - 1):
        (row, p) = nextpatient(sheet, row)
        if filter_fn(p):
            yield f(p)

def deltaamplitude(data):
    min_threshold = 0#15 / 60.0 
    max_threshold = 60
    a_threshold = -1000
    f_threshold = lambda t,a,da: ((da > 0) and (a > a_threshold) and (t > min_threshold) and (t < max_threshold))
    f_dA = lambda a1,a2: (a2[0], 0.5*(a2[1]+a1[1]), a2[1]-a1[1])
    return [f_dA(d1, d2)  for (d1, d2) in zip(data[:-1], data[1:]) if f_threshold(f_dA(d1, d2)[0], f_dA(d1, d2)[1], f_dA(d1, d2)[2])]

def maxdeltaamplitude(p):
    dA = deltaamplitude(p['data'])
    m = None
    if (len(dA) > 0):
        m = max(dA, key=lambda x: x[2])
    #return (m[0][0], m[0][1], m[1][1])
    return m

def writepatientinfo(writer_file, p):
    m = maxdeltaamplitude(p)
    s = ', '.join([str(int(p['fullname'])), str(p['date']), str(p['sampletype']), str(p['description']), str(m[2]), str(m[1])])
    writer_file.write(s + '\n')

def patientwriter(writer_file):
    return lambda p: writepatientinfo(writer_file, p)

def maxamplitude(p):
    min_threshold = 15/60.0
    max_threshold = 60
    f_threshold = lambda t: t > min_threshold and t < max_threshold
    return max([(t,d) for (t,d) in p['data'] if f_threshold(t)], key=lambda x: x[1])

amount_after_acrit = 10

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
    amps = [a for t,a,da in dAs if a <= Acrit+amount_after_acrit]
    ts = [t for t,a,da in dAs if a <= Acrit + amount_after_acrit]
    #ax.plot(ts, amps, c='r', linewidth=3)
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
    ax.scatter([a for (t,a,da) in dA], [da for (t,a,da) in dA])
    area, d = dAintegration(p)
    plotbestfit(p, ax, numpy.Inf, 'g', False)
    plotbestfit(p, ax, 0, 'r', True)
    #plotbestfit(p, ax, amount_after_acrit, 'r', True)
    ax.legend(['Fit up to Amax', 'Fit up to Acrit'])
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

def superplot_tk(p):
    fig, (a1, a2) = plt.subplots(2, 1)
    plottimeseriestoaxes(p, a1)
    plotdavatoaxes(p, a2)
    plt.show()

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
        print("No dA values. Problem with this patient: %s %s %s %s"%(str(out['fullname']), makeexceldatestring(out['date']), out['sampletype'], out['description']))
    return out

def findallpatiententries(sheet, patientid, f=simpledetails, filter_fn = lambda x:True):
    f_patient = lambda p: p['fullname'] == patientid
    #or int(p['fullname']) == patientid
    g = patientgenerator(sheet, f, lambda p: f_patient(p) and filter_fn(p))
    return [m for m in g]

def makePDFpage(pdfcontext, patientid,  sheet, filter_fn = lambda x:True):
    patients = findallpatiententries(sheet, patientid, lambda p:p, filter_fn)
    patients = sorted(patients, key=lambda x: x['date'])
    for p in patients:
        fig, (ax1, ax2)  = plt.subplots(2,1)
        title = "Patient ID: %s, date: %s, sample type: %s, description: %s" %(str(p['fullname']), makeexceldatestring(p['date']), p['sampletype'], p['description'])
        fig.suptitle(title)
        plottimeseriestoaxes(p, ax1)
        plotdavatoaxes(p, ax2)
        pdfcontext.savefig(fig)
        
def makefullPDF(pdffilename, sheet, filter_fn = lambda x:True):
    g = patientgenerator(sheet, simpledetails, filter_fn)
    ids = set([i['fullname'] for i in g if i != None])
    pp = PdfPages(pdffilename)
    for i in ids:
        makePDFpage(pp, i, sheet, filter_fn)
    pp.close()
    plt.close('all')
        
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

def plotgroups(grouped_patients):
    colors = ['k', 'b', 'g', 'r', 'm', 'y']
    legend_values = []
    for (k, col) in zip(grouped_patients.keys(), colors):
        x = [i['Acrit'] for i in grouped_patients[k]]
        y = [i['dAmax'] for i in grouped_patients[k]]
        plt.scatter(x, y, c=col)#, alpha=0.75)
        legend_values.append(k)
    plt.xlabel('Amplitude')
    plt.ylabel('Acrit')
    plt.legend(legend_values)
    plt.title("Acrit versus amplitude grouped by sample type")
    plt.show()

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
    #d = data_after_split(p)
    d = p['data']
    dAs = deltaamplitude(d)
    return fit_best_dAvA_curve([a for t,a,da in dAs], [da for t,a,da in dAs])

def fit_best_uptoAcrit(p, amount_after_Acrit = 0):
    #d = data_after_split(p)
    d = p['data']
    dAs = deltaamplitude(d)
    tdAmax, Acrit, dAmax = maxdeltaamplitude(p)
    amps = [a for t,a,da in dAs if a <= Acrit + amount_after_Acrit]
    davals = [da for t,a,da in dAs if a <= Acrit + amount_after_Acrit]
    return fit_best_dAvA_curve(amps, davals)

def plot_best_fit(p):
    a,b,minval = fit_best(p)
    dAs = deltaamplitude(p['data'])
    plt.plot([a for t,a,da in dAs], [da for t,a,da in dAs])
    avals = numpy.linspace(1, 60, 100)
    plt.plot(avals, solve_dA(a,b)(avals))
    plt.show()

def calcExpectedAmax(a, b):
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
    return ((a,b), dAintegration(p)[0], calcExpectedAmax(a,b), numericIntegrationOfModel(solve_dA(a,b), 0.000001, 100), r2)
    
# "patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, numerical integration to 100mm, analytical integral of model"

def area_details(p):
    e = simpledetails(p)
    a1, b1, e1 = fit_best(p)
    a2, b2, e2 = fit_best_uptoAcrit(p)
    #d = data_after_split(p)
    d=p['data']
    dAs = deltaamplitude(d)
    r21 = corr_coef([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a1, b1))
    r22 = corr_coef([a for t,a,da in dAs], [da for t,a,da in dAs], solve_dA(a2, b2))
    area = dAintegration(p)[0]
    n1 = numericIntegrationOfModel(solve_dA(a1, b1), 0.000001, 100)
    n2 = numericIntegrationOfModel(solve_dA(a2, b2), 0.000001, 100)
    c1 = calcExpectedAmax(a1, b1)
    c2 = calcExpectedAmax(a2, b2)
    return (area, (a1, b1, r21, n1, c1), (a2, b2, r22, n2, c2))
    #s = ', '.join([str(n), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']))])

def stringify_area_stuffs(p):
    area, full, partial = area_details(p)
    e = simpledetails(p)
    name = e['fullname']
    if isinstance(name, (int, float, complex)):
        name = int(name)
        
    s = ', '.join([str(name), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']), str(area), str(full[0]), str(full[1]), str(full[2]), str(full[3]), str(full[4]), str(partial[0]), str(partial[1]), str(partial[2]), str(partial[3]), str(partial[4])])
    meta = "patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, r2, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, r2, numerical integration to 100mm, analytical integral of model"
    return meta, s
    
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

foo="""""
def tabulate_data(sheet, outfilename, sampletype):
    g = patientgenerator(sheet, simpledetails, lambda p: p['sampletype'] == sampletype)
    grouped_patients = groupbypatient([p for p in g])
    names = grouped_patients.keys()
    names = sorted(names)
    
    outfile = open(outfilename, 'w')
    outfile.write("patient ID, date-time, sample type, description, t @ Amax (min), Amax (mm), tcrit (min), Acrit (mm), dAmax (mm), tsplit (min), Asplit (mm), empirical area under dA v A curve, a fit to Amax, b fit to Amax, numerical integration to 100mm, analytical integral of model, a fit to Acrit, b fit to Acrit, numerical integration to 100mm, analytical integral of model\n")

    for n in names:
        if isinstance(n, (int, float, complex)):
            n = int(n)
        for e in grouped_patients[n]:
            try:
                grabareastuff
                s = ', '.join([str(n), str(makeexceldatestring(e['date'])), str(e['sampletype']), str(e['description']), str(e['tmax']), str(e['Amax']), str(e['tcrit']), str(e['Acrit']), str(e['dAmax']), str(e['tsplit']), str(e['asplit']))])
                outfile.write(s + "\n")
            except:
                print("Problem with patient: %s, %s, %s, %s" % (str(e['fullname']), makeexceldatestring(e['date']), e['sampletype'], e['description']))
    outfile.close()
"""