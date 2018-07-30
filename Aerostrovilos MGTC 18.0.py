from tkinter import *
from math import *
win=Tk()
'''entry .s -width 1
entry .t -width 9
grid .s .t -sticky ew
grid columnconf . 0 -weight 10
grid columnconf . 1 -weight 90'''
icon_2 = PhotoImage(file = "logo.png")
win.call('wm','iconphoto',win._w,icon_2)
win.geometry("1200x675")
win.wm_title("Aerostrovilos MGTC 18.0")
rows=0
win.rowconfigure(0,weight=1)
win.columnconfigure(0,weight=1)


canvas=Canvas(win,bg='#e5e7e3')
canvas.grid(sticky='nesw')
mycolor='#%02x%02x%02x' % (64, 204, 208)
#HE
canvas.create_rectangle(150,85,300,150, outline='gray',width=2,fill=mycolor)
canvas.create_rectangle(200,105,250,130, outline='gray',width=2,fill='white')
canvas.create_text(225,120,text="HE")
canvas.create_line(300,85,700,85,width=3,arrow='first')
canvas.create_line(300,150,500,150,width=3,arrow='last')

#COM
canvas.create_rectangle(500,150,575,175,fill='#d60000')
canvas.create_text(535,165,text="COM",font='Ariel')
canvas.create_line(700,85,700,375,width=3,arrow='first')
canvas.create_line(575,175,575,400,width=3,arrow='last')

#T
canvas.create_polygon(700,375,575,400,575,475,700,500,width=2,fill='#ed710b')
canvas.create_rectangle(610,420,660,445,fill='white')
canvas.create_text(635,435,text="T")

#c
canvas.create_line(275,150,275,400,width=3,arrow='first')
canvas.create_polygon(275,400,150,375,150,500,275,475,fill='#0320a0',width=2)
canvas.create_rectangle(240,420,190,445,fill='white')
canvas.create_text(215,435,text="C")
canvas.create_line(600,435,250,435,width=7,arrow='both')

#G
canvas.create_line(700,435,900,435,width=7,arrow='last')
canvas.create_oval(900,485,1000,385,fill='gray')
canvas.create_rectangle(925,420,975,445,fill='white')
canvas.create_text(950,435,text="G")

#extra lines
canvas.create_line(150,120,100,120,width=3,arrow='last')
canvas.create_line(150,435,100,435,width=3,arrow='first')
canvas.create_line(1000,435,1050,435,width=5,arrow='last',fill='red')

#creating legend
canvas.create_text(900,600,text="All temperatures are in kelvin(K) and pressures in kPa",width=250,anchor='w',font='times,12')

#calculations
def calculate():
    deltc.set((pow(pr.get(),(gamma.get()-1)/gamma.get())-1)*to1.get()/nc.get())
    deltci=(deltc.get()*nc.get())
    wci.set(cpa.get()*deltci)
    wc.set(cpa.get()*deltc.get())
    to2.set(deltc.get()+to1.get())
    po2.set(po1.get()*pr.get())
    po3.set(po2.get()*(1-delpxa.get()))
    drop.set(po2.get()-po3.get())
    po4.set(po3.get()*(1-delpcom.get()))
    po6.set(po1.get()/(1-delpke.get()))
    po5.set(po6.get()/(1-delpxg.get()))
    pet.set(po4.get()/po5.get())
    deltt.set((1-pow(1/pet.get(),(gammat.get()-1)/gammat.get()))*nt.get()*to4.get())
    to5.set(to4.get()-deltt.get())
    deltti=deltt.get()/nt.get()

    wti.set(cpg.get()*deltti)
    wt.set(cpg.get()*deltt.get())
    wtc.set(wc.get()/nm.get())
    wtb.set(wt.get()*(1-pb.get()))
    wts.set(wt.get()-wtc.get()-wtb.get())
    wsg.set(wts.get()*nm.get())
    we.set(wsg.get()*ng.get())
    m.set(pe.get()/we.get())
    to3.set(to2.get()+epl.get()*(to5.get()-to2.get()))

    to6.set(to5.get()-cpaex.get()/cpgex.get()*(to3.get()-to2.get()))
    deltxa.set(to3.get()-to2.get())
    deltxg.set(to5.get()-to6.get())
    '''cycle efficiency calculation'''
    deltcom.set(to4.get()-to3.get())
    q.set(m.get()*cpg.get()*deltcom.get())
    mft=(q.get()/lhv.get())
    ft.set(mft/m.get())
    f.set(ft.get()/ncom.get())
    mf.set(f.get()*m.get())
    fphi=0.058
    phi.set(f.get()/fphi)
    ncycle.set(pe.get()*100/mf.get()/lhv.get())
    pf.set(m.get()*we.get())
    pt.set(m.get()*wt.get())
    ptb.set(m.get()*wtb.get())
    pc.set(m.get()*wc.get())
    roundoff()

def roundoff():
    roundlist4=["deltc","wci","wc","to2","po2","po3","drop","po4","po6","po5","pet"
    ,"deltt","to5","wti","wt","wtc","wtb","wts","wsg","we","to3","to6","deltxa",
    "deltxg","deltcom","pf","pt","ptb","pc"]
    roundlist7=["m","q","ft","f","phi","ncycle"]
    for name4 in roundlist4:
        exec("%s.set(round(%s.get(),4))"%(name4,name4))
    for name7 in roundlist7:
        exec("%s.set(round(%s.get(),7))"%(name7,name7))


#enteries and labels
delpke=DoubleVar()
entrydelpke=Entry(canvas,textvariable=delpke)
entrydelpke.config(bd=2,width=10)
canvas.create_window(80,70,window=entrydelpke,)
canvas.create_text(20,70,text=u"\u2206""Pke",anchor='w',font='Times')

po6=DoubleVar()
entrypo6=Entry(canvas,textvariable=po6)
entrypo6.config(bd=2,width=10)
canvas.create_window(80,100,window=entrypo6)
canvas.create_text(20,100,text="Po6",anchor='w',font='Times')

to6=DoubleVar()
entryto6=Entry(canvas,textvariable=to6)
entryto6.config(bd=2,width=10)
canvas.create_window(80,140,window=entryto6)
canvas.create_text(20,140,text="To6",anchor='w',font='Times')

tup1=[("deltxg","entrydeltxg",u"\u2206""Txg"),("delpxg","entrydelpxg",u"\u2206""Pxg"),("epl","entryepl",u"\u03B5")]
h=10
for name,entry,text in tup1:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(250,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(165,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup2=[("delpcom","entrydelpcom",u"\u2206""Pcom"),("ncom","entryncom",u"\u03B7""com")]
h=100
for name,entry,text in tup2:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(590,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(490,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup3=[("delpxa","entrydelpxa",u"\u2206""Pxa"),("drop","entrydrop",""),("cpaex","entrycpaex","Cpaex")]
h=170
for name,entry,text in tup3:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(220,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(125,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup4=[("po3","entrypo3","Po3"),("to3","entryto3","To3"),("deltxa","entrydeltxa",u"\u2206""Txa")]
h=170
for name,entry,text in tup4:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=10)"%(entry))
    exec("canvas.create_window(370,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(300,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup5=[("po5","entrypo5","Po5"),("to5","entryto5","To5"),("cpgex","entrycpgex","Cpgex")]
h=150
for name,entry,text in tup5:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=10)"%(entry))
    exec("canvas.create_window(790,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(710,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup6=[("po2","entrypo2","Po2"),("to2","entryto2","To2")]
h=330
for name,entry,text in tup6:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(220,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(140,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup7=[("po4","entrypo4","Po4"),("to4","entryto4","To4"),("deltcom","entrydeltcom",u"\u2206""Tcom")]
h=320
for name,entry,text in tup7:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(525,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(425,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup8=[("cpg","entrycpg","Cpg"),("gammat","entrygammat",u"\u0263""t")]
h=320
for name,entry,text in tup8:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=10)"%(entry))
    exec("canvas.create_window(650,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(590,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup9=[("po1","entrypo1","Po1"),("to1","entryto1","To1"),("cpa","entrycpa","Cpair"),("gamma","entrygamma",u"\u0263""air")]
h=391
for name,entry,text in tup9:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=10)"%(entry))
    exec("canvas.create_window(90,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(20,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup10=[("pr","entrypr","Pr"),("nc","entrync",u"\u03B7""c"),("deltc","entrydeltc",u"\u2206""Tc"),("wc","entrywc","Wc"),("wci","entrywci","Wci")]
h=520
for name,entry,text in tup10:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(235,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(155,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup11=[("nt","entrynt",u"\u03B7""t"),("pet","entrypet","Pet"),("deltt","entrydeltt",u"\u2206""Tt"),("wt","entrywt","Wt"),("wti","entrywti","Wti")]
h=520
for name,entry,text in tup11:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(645,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(565,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup12=[("pe","entrype","Pe(KW)"),("m","entrym","m(kg/s)"),("ft","entryft","ft"),("f","entryf","f"),
("mf","entrymf","mf(kg/s)"),("q","entryq","Q"),("ncycle","entryncycle",u"\u03B7""cycle")]
h=185
for name,entry,text in tup12:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(975,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(875,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup13=[("nb","entrynb",u"\u03B7""b"),("wtb","entrywtb","Wtb"),("ptb","entryptb","Ptb")]
h=465
for name,entry,text in tup13:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(495,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(410,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

tup14=[("wts","entrywts","Wts"),("wsg","entrywsg","Wsg"),("ng","entryng",u"\u03B7""g"),]
h=465
wid=715
for name,entry,text in tup14:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(wid+80,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(wid,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    wid+=80
    h+=30

tup15=[("we","entrywe","We"),("pf","entrypf","Pf(KW)")]
h=420
for name,entry,text in tup15:
    exec("%s = DoubleVar()" % (name))
    exec("%s=Entry(canvas,textvariable=%s)" % (entry,name))
    exec("%s.config(bd=2,width=15)"%(entry))
    exec("canvas.create_window(1135,%d,window=%s)"%(h,entry))
    exec("canvas.create_text(1035,%d,text='"%h + "%s"%text + "',anchor='w',font='Times')")
    h+=30

lhv=DoubleVar()
entrylhv=Entry(canvas,textvariable=lhv)
entrylhv.config(bd=2,width=15)
canvas.create_window(945,120,window=entrylhv)
canvas.create_text(865,120,text="LHV",anchor='w',font='Times')

wtc=DoubleVar()
entrywtc=Entry(canvas,textvariable=wtc)
entrywtc.config(bd=2,width=10)
canvas.create_window(360,415,window=entrywtc)
canvas.create_text(290,415,text="Wtc",anchor='w',font='Times')

nm=DoubleVar()
entrynm=Entry(canvas,textvariable=nm)
entrynm.config(bd=2,width=10,bg="yellow")
canvas.create_window(515,415,window=entrynm)
canvas.create_text(440,415,text=u"\u03B7""m",anchor='w',font='Times')

pc=DoubleVar()
entrypc=Entry(canvas,textvariable=pc)
entrypc.config(bd=2,width=15)
canvas.create_window(360,610,window=entrypc,anchor='w')
canvas.create_text(290,610,text="Pc(KW)",anchor='w',font='Times')

pt=DoubleVar()
entrypt=Entry(canvas,textvariable=pt)
entrypt.config(bd=2,width=15)
canvas.create_window(760,610,window=entrypt,anchor='w')
canvas.create_text(700,610,text="Pt(KW)",anchor='w',font='Times')

phi=DoubleVar()
entryphi=Entry(canvas,textvariable=phi)
entryphi.config(bd=2,width=15)
canvas.create_window(1140,250,window=entryphi)
canvas.create_text(1080,250,text=u"\u03D5",anchor='w',font='Times')

texttup=[(135,425),(265,265),(400,140),(585,275),(710,245),(135,110)]
j=1
for x,y in texttup:
    canvas.create_text(x,y,text=j,fill="red",font="Ariel")
    j+=1

but=Button(text="Calculate",command=calculate)
canvas.create_window(1100,100,window=but)


pb=DoubleVar()

def usedefault():
    pr.set(4)
    nc.set(0.8)
    to1.set(300)
    gamma.set(1.4)
    cpa.set(1.005)
    po1.set(101.33)
    delpxa.set(0.04)
    delpcom.set(0.05)
    delpke.set(0.05)
    nt.set(0.8)
    to4.set(1200)
    gammat.set(1.323)
    delpxg.set(0.04)
    cpg.set(1.184)
    nm.set(0.99)
    pb.set(0.96)
    nb.set(0.96)
    ng.set(0.92)
    pe.set(100)
    epl.set(0.9)
    ncom.set(0.99)
    lhv.set(47500)
    cpaex.set(1.071)
    cpgex.set(1.1)

bgtup=("pr","nc","to1","gamma","cpa","po1","delpxa","delpcom","delpke","nt","to4","gammat","delpxg",
"cpg","nm","pb","nb","ng","pe","epl","ncom","lhv","cpaex","cpgex")

for item in bgtup:
    if item != "pb":
        exec("entry%s.config(bg='yellow')"%item)
entryncycle.config(bg='#02cae0')

but2=Button(text="Use default values",command=usedefault)
canvas.create_window(1100,50,window=but2)

canvas.addtag_all('all')




win.mainloop()
