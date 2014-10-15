import xlrd
import matplotlib.pyplot as plt

def opensheet(excel_file_name):
    wb = xlrd.open_workbook(excel_file_name)
    return wb.sheets()[0]

def nextpatient(sheet, start_line):
    name = sheet.cell(start_line, 0).value
    i = start_line + 1
    n = sheet.cell(i, 0).value
    d = []
    while n == '':
        d.append((sheet.cell(i,7).value, sheet.cell(i,8).value / 2))
        i = i + 1
        n = sheet.cell(i, 0).value
    return (i, {'fullname': name, 'data': d})

def patientgenerator(sheet, f = lambda x:x):
    last_row = sheet.nrows
    row = 1
    while row < last_row:
        (row, p) = nextpatient(sheet, row)
        yield f(p)

def deltaamplitude(data):
    f_dA = lambda a1,a2: (a2[0], (a2[1] - a1[1]) / (a2[0] - a1[0]))
    return [f_dA(d1, d2)  for (d1, d2) in zip(data[:-1], data[1:])]

def maxdeltaamplitude(p):
    dA = deltaamplitude(p['data'])
    dtot = zip(p['data'][1:], dA)
    m = max(dtot, key=lambda x: x[1][1])
    return (m[0][0], m[0][1], m[1][1])

def plottimeseriesdata(p):
    m = maxdeltaamplitude(p)
    plt.title('Patient ' + str(int(p['fullname'])))
    plt.plot([i[0] for i in p['data']], [i[1] for i in p['data']])
    plt.scatter(m[0], m[1], c='r')
    plt.annotate('Acrit', xy=(m[0],m[1]), xytext=(m[0]+80,m[1]))
    plt.xlabel("Time (sec)")
    plt.ylabel("Blood Fluid Resistance Measurement (units?)")
    plt.show()
