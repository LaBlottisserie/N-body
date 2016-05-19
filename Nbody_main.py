# -*- coding: iso-8859-1 -*-

from visual import *
import thread
import numpy as np
from visual.graph import *
import numpy.random as alea
from math import *
import wx
import os
import sys

version="beta 1.1.6"

aide_thread = 0
propos_thread = 0

def galaxie(event):
    global scene, sw, speed
    sw = 1
    scene.delete()
    scene = display(window =w,width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.title = "GALAXIE"
    scene.fov = pi / 2
    scene.scale = (0.0005,0.0005,0.0005)
    G = 1  # constante gravitationnelle
    dt = 1  # pas de temps
    n = 100  # nombre de corps
    m = 2  # masse des corps
    L = 500  # Longueur caractèristique
    speed = 1000
    slid1.SetValue(speed)
    strspeed.SetLabel(str(speed))
    e = 35  # epaisseur du bulbe galactique
    rayon = 10
    mt = 20000

    def color_galaxie(l):
        vitesses = []
        for i in range(len(l)):
            v = sqrt(l[i].v.x ** 2 + l[i].v.y ** 2 + l[i].v.z ** 2)
            vitesses.append(v)
            vmax = max(vitesses)
            vmin = min(vitesses)

        for i in range(len(l)):
            v = sqrt(l[i].v.x ** 2 + l[i].v.y ** 2 + l[i].v.z ** 2)
            l[i].color = ((v - vmin) / (vmax - vmin), 0.3, (vmax-v) / (vmax - vmin))

    def nuageradial(corps, maxi):
        corps.x = alea.normal(0, maxi)
        corps.y = alea.normal(0, maxi)
        corps.z = alea.normal(0, e)

    def vitesseradiale(corps):
        d = sqrt(corps.x * corps.x + corps.y * corps.y + corps.z * corps.z)
        vit = sqrt(100 * n * m * G / d)
        corps.v = vector(-corps.y / d * vit, corps.x / d * vit, 0)
        corps.m = m

    def MAJposition(corps):
        corps.pos = corps.pos + (corps.v) * dt

    def acceleration(C):
        D = np.zeros((n, n))
        A = np.zeros((n, 3))
        for i in range(n):
            for k in range(n):
                corps1 = C[i]
                corps2 = C[k]
                dx = corps2.x - corps1.x
                dy = corps2.y - corps1.y
                dz = corps2.z - corps1.z
                d2 = dx * dx + dy * dy + dz * dz
                if d2 == 0:
                    D[k, i] = 0
                else:
                    D[k, i] = D[i, k] = d2 * sqrt(d2)
        for i in range(n):
            for k in range(n):
                if k != i:
                    corps1 = C[i]
                    corps2 = C[k]
                    A[i, 0] += G * corps2.m * (corps2.x - corps1.x) / D[i, k]
                    A[i, 1] += G * corps2.m * (corps2.y - corps1.y) / D[i, k]
                    A[i, 2] += G * corps2.m * (corps2.z - corps1.z) / D[i, k]
        return A

    def MAJvitesse(corps, mat):
        corps.v.x += mat[0] * dt
        corps.v.y += mat[1] * dt
        corps.v.z += mat[2] * dt

    def creercorps():
        nouveaucorps = sphere(radius=rayon)
        nuageradial(nouveaucorps, L)
        vitesseradiale(nouveaucorps)
        return nouveaucorps

    def nbody():
        ###creation d un astre central massif
        trounoir = sphere(radius=rayon * 4, color=color.yellow)
        trounoir.v = vector(0, 0, 0)
        trounoir.x = 0
        trounoir.y = 0
        trounoir.z = 0
        trounoir.m = mt
        Corps = [trounoir]
        for i in range(n-1):
            Corps.append(creercorps())
        while sw==1:
            M = acceleration(Corps)
            rate(speed)
            for i in range(n):
                corps = Corps[i]
                MAJvitesse(corps, M[i])
                MAJposition(corps)
            color_galaxie(Corps[1:])


    nbody()

def threebody(event):
    global G, scene,sw,speed
    sw = 2
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale=(0.5,0.5,0.5)
    scene.title = "3 corps"
    scene.fov = pi / 2
    # définition des paramètres
    G = 1  # constante gravitationnelle
    dt = 0.001  # pas de temps
    n = 3  # nombre de corps
    m = 1  # masse des corps
    L = 500  # Longueur caract?ristique
    speed = 500  # higher is slower
    slid1.SetValue(speed)
    strspeed.SetLabel(str(speed))
    rayon = 0.05

    def MAJposition(corps):
        corps.pos = corps.pos + (corps.v) * dt

    def acceleration(C):
        D = np.zeros((n, n))
        A = np.zeros((n, 3))
        for i in range(n):
            for k in range(n):
                corps1 = C[i]
                corps2 = C[k]
                dx = corps2.x - corps1.x
                dy = corps2.y - corps1.y
                dz = corps2.z - corps1.z
                d2 = dx * dx + dy * dy + dz * dz
                if i == k:
                    D[k, i] = 0
                elif sqrt(d2) <= 0.01:  # terme correctif pour ?viter l'effet catapulte
                    D[k, i] = 0.01
                else:
                    D[k, i] = D[i, k] = d2 * sqrt(d2)
        for i in range(n):
            for k in range(n):
                if k != i:
                    corps1 = C[i]
                    corps2 = C[k]
                    A[i, 0] += G * corps2.m * (corps2.x - corps1.x) / D[i, k]
                    A[i, 1] += G * corps2.m * (corps2.y - corps1.y) / D[i, k]
                    A[i, 2] += G * corps2.m * (corps2.z - corps1.z) / D[i, k]
        return A

    def MAJvitesse(corps, mat):
        corps.v.x += mat[0] * dt
        corps.v.y += mat[1] * dt
        corps.v.z += mat[2] * dt

    corps0 = sphere(radius=rayon, make_trail=False, color=color.cyan,retain=2100)
    corps1 = sphere(radius=rayon, make_trail=False, color=color.white,retain=2100)
    corps2 = sphere(radius=rayon, make_trail=False, color=color.red,retain=2100)
    corps0.m = 1
    corps1.m = 1
    corps2.m = 1

    ###definition des 3corps
    threebodies = [corps0,corps1,corps2]
    C0 = threebodies[0]
    C1 = threebodies[1]
    C2 = threebodies[2]
    ###vitesses
    C0.v = vector(-0.93240737 / 2, -0.86473146 / 2)
    C1.v = vector(-0.93240737 / 2, -0.86473146 / 2)
    C2.v = vector(0.93240737, 0.86473146, 0)
    ###positions
    C0.x, C0.y, C0.z = -0.97000436, 0.24308753, 0
    C1.x, C1.y, C1.z = 0.97000436, -0.24308753, 0
    C2.x, C2.y, C2.z = 0, 0, 0

    def nbody(Corps):
        ###Vpython settings###
        pid=0
        while sw ==2:
            if pid==2:
                corps0.make_trail = True
                corps0.trail_type = "curve"
                corps1.make_trail = True
                corps1.trail_type = "curve"
                corps2.make_trail = True
                corps2.trail_type = "curve"
            pid+=1
            M = acceleration(Corps)
            rate(speed)
            for i in range(n):
                corps = Corps[i]
                MAJvitesse(corps, M[i])
                MAJposition(corps)

    nbody(threebodies)

def terrelune(event):
    global G, scene, sw,speed
    sw = 3
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale = (0.0000000013, 0.0000000013, 0.0000000013)
    scene.fov = pi / 2
    # définition des paramètres
    G = 6.67384 * 10 ** -11  # constante gravitationnelle
    dt = 600  # pas de temps
    n = 2  # nombre de corps
    m = 1  # masse des corps              #Longueur caractèristique
    speed = 150  # higher is slower
    slid1.SetValue(speed)
    strspeed.SetLabel(str(speed))

    def MAJposition(corps):
        corps.pos = corps.pos + (corps.v) * dt

    def acceleration(C):
        D = np.zeros((n, n))
        A = np.zeros((n, 3))
        for i in range(n):
            for k in range(n):
                corps1 = C[i]
                corps2 = C[k]
                dx = corps2.x - corps1.x
                dy = corps2.y - corps1.y
                dz = corps2.z - corps1.z
                d2 = dx * dx + dy * dy + dz * dz
                if i == k:
                    D[k, i] = 0
                elif sqrt(d2) <=0.01:  # terme correctif pour ?viter l'effet catapulte
                    D[k, i] = 0.01
                else:
                    D[k, i] = D[i, k] = d2 * sqrt(d2)
        for i in range(n):
            for k in range(n):
                if k != i:
                    corps1 = C[i]
                    corps2 = C[k]
                    A[i, 0] += G * corps2.m * (corps2.x - corps1.x) / D[i, k]
                    A[i, 1] += G * corps2.m * (corps2.y - corps1.y) / D[i, k]
                    A[i, 2] += G * corps2.m * (corps2.z - corps1.z) / D[i, k]
        return A

    def MAJvitesse(corps, mat):
        corps.v.x += mat[0] * dt
        corps.v.y += mat[1] * dt
        corps.v.z += mat[2] * dt

    ###definition des 2corps

    Corps0 = sphere(radius=6300000, make_trail=False, color=(0.1, 0.4, 0.8))
    Corps1 = sphere(pos=(0, 384000000, 0), make_trail=False, radius=1700000, color=(0.9, 0.9, 0.9))
    TerreLune = [Corps0, Corps1]
    texte1 = text(text='Terre', align='center', depth=0.1, color=(0.1, 0.4, 0.8), height=20000000, font="Times")
    texte2 = text(text='Lune', align='center', depth=0.1, color=(0.9, 0.9, 0.9), height=20000000, font="Times")
    # masses
    Corps0.m = 5.972 * 10 ** 24
    Corps1.m = 7.35 * 10 ** 22

    ###vitesses
    Corps0.v = vector(0, 0, 0)
    Corps1.v = vector(0, 1052, 0)
    ###positions
    Corps0.x, Corps0.y, Corps0.z = 0, 0, 0
    Corps1.x, Corps1.y, Corps1.z = 384000000, 0, 0

    def nbody(Corps):
        pid = 0
        while sw == 3:
            if pid == 2:
                Corps0.make_trail = True
                Corps0.trail_type = "curve"
                Corps0.retain = 10000
                Corps1.make_trail = True
                Corps1.trail_type = "curve"
                Corps1.retain = 10000
            pid += 1
            texte1.pos = Corps0.pos + (0, 1.6e7, 0)
            texte2.pos = Corps1.pos - (0, 2e7, 0)
            rate(speed)
            M = acceleration(Corps)
            for i in range(n):
                corps = Corps[i]
                MAJvitesse(corps, M[i])
                MAJposition(corps)


    nbody(TerreLune)

def saturne(event):
    global G, scene, sw, speed
    sw = 4
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale = (0.05, 0.05, 0.05)
    scene.fov = pi / 2
    G = 1  # constante gravitationnelle
    dt = 0.00008  # pas de temps
    n = 6000  # nombre de corps
    m = 0.1  # masse des corps
    L = 10  # Longueur caract?ristique
    speed = 1000
    slid1.SetValue(speed)
    strspeed.SetLabel(str(speed))
    rayon = 0.05
    rsaturne = 4
    Msat = 1000000

    def nuageradial(corps, maxi):
        teta = random.uniform(0, 2 * pi)
        p = random.uniform(rsaturne * 1.3, rsaturne * 1.4)  # premier anneau
        q = random.uniform(rsaturne * 1.45, rsaturne * 1.6)  # deuxi?me anneau
        s = random.uniform(rsaturne * 1.8, rsaturne * 2.4)  # troisi?me anneau
        L = [p, q, s]
        d = L[random.randint(0, 3)]
        corps.x = d * cos(teta)
        corps.y = d * sin(teta)
        corps.z = 0

    def vitesseradiale(corps):
        d = sqrt(corps.x * corps.x + corps.y * corps.y + corps.z * corps.z)
        vit = sqrt(Msat * G / d)
        corps.v = vector(-corps.y / d * vit, corps.x / d * vit, 0)
        corps.m = m

    def MAJposition(corps):
        corps.pos = corps.pos + (corps.v) * dt

    def acceleration(C):
        D = []
        A = np.zeros((n, 3))
        for i in range(1, n):
            corps1 = C[i]
            dx = corps1.x
            dy = corps1.y
            dz = corps1.z
            d2 = dx * dx + dy * dy + dz * dz
            D.append(d2 * sqrt(d2))
        for i in range(1, n):
            corps1 = C[i]
            A[i, 0] = G * Msat * (-corps1.x) / D[i - 1]
            A[i, 1] = G * Msat * (-corps1.y) / D[i - 1]
            A[i, 2] = G * Msat * (-corps1.z) / D[i - 1]
        return A

    def MAJvitesse(corps, mat):
        corps.v.x += mat[0] * dt
        corps.v.y += mat[1] * dt
        corps.v.z += mat[2] * dt

    def creercorps():
        nouveaucorps = sphere(radius=rayon, color=color.orange)
        nuageradial(nouveaucorps, L)
        vitesseradiale(nouveaucorps)
        return nouveaucorps

    def nbody():
        ###creation d un astre central massif
        Saturne = sphere(radius=rsaturne, color=color.orange)
        Saturne.v = vector(0, 0, 0)
        Saturne.x = 0
        Saturne.y = 0
        Saturne.z = 0
        Saturne.m = Msat
        Corps = [Saturne]
        for i in range(n):
            Corps.append(creercorps())
        while sw == 4:
            M = acceleration(Corps)
            rate(speed)
            for i in range(n):
                corps = Corps[i]
                MAJvitesse(corps, M[i])
                MAJposition(corps)


    nbody()

def sphere__(event):
    global G, scene, sw, speed
    sw = 5
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale = (0.05, 0.05, 0.05)
    scene.fov = pi / 2
    # définition des paramètres
    G = 1  # constante gravitationnelle
    dt = 0.5  # pas de temps
    n = 120  # nombre de corps
    m = 0.01  # masse des corps
    L = 20  # Longueur caract?ristique
    speed = 200
    slid1.SetValue(speed)
    strspeed.SetLabel(str(speed))
    rayon = 0.5

    def nuagespherique(corps, maxi):
        teta = random.uniform(0, 2 * pi)
        phi = random.uniform(0, 2 * pi)
        corps.x = L * cos(teta) * cos(phi)
        corps.y = L * cos(teta) * sin(phi)
        corps.z = L * sin(teta)

    def MAJposition(corps):
        corps.pos = corps.pos + (corps.v) * dt

    def vitesseinitiale(corps):
        corps.v = vector(0, 0, 0)

    def acceleration(C):
        D = np.zeros((n, n))
        A = np.zeros((n, 3))
        for i in range(n):
            for k in range(n):
                corps1 = C[i]
                corps2 = C[k]
                dx = corps2.x - corps1.x
                dy = corps2.y - corps1.y
                dz = corps2.z - corps1.z
                d2 = dx * dx + dy * dy + dz * dz
                if i == k:
                    D[k, i] = 0
                elif sqrt(d2) <= corps1.radius + corps2.radius:  # terme correctif pour ?viter l'effet catapulte
                    D[k, i] = corps1.radius + corps2.radius
                else:
                    D[k, i] = D[i, k] = d2 * sqrt(d2)
        for i in range(n):
            for k in range(n):
                if k != i:
                    corps1 = C[i]
                    corps2 = C[k]
                    A[i, 0] += G * corps2.m * (corps2.x - corps1.x) / D[i, k]
                    A[i, 1] += G * corps2.m * (corps2.y - corps1.y) / D[i, k]
                    A[i, 2] += G * corps2.m * (corps2.z - corps1.z) / D[i, k]
        return A

    def MAJvitesse(corps, mat):
        corps.v.x += mat[0] * dt
        corps.v.y += mat[1] * dt
        corps.v.z += mat[2] * dt

    def creercorps():
        nouveaucorps = sphere(radius=rayon, color=color.green)
        nuagespherique(nouveaucorps, L)
        vitesseinitiale(nouveaucorps)
        nouveaucorps.m = m
        return nouveaucorps

    def nbody():
        ###creation d un astre central massif
        trounoir = sphere(radius=rayon * 2, color=color.red)
        trounoir.v = vector(0, 0, 0)
        trounoir.x = 0
        trounoir.y = 0
        trounoir.z = 0
        trounoir.m = 1000 * m
        Corps = [trounoir]
        for i in range(n - 1):
            Corps.append(creercorps())
        scene.scale = (0.022,0.022,0.022)
        while sw == 5:
            M = acceleration(Corps)
            rate(speed)
            for i in range(n):
                corps = Corps[i]
                MAJvitesse(corps, M[i])
                MAJposition(corps)


    nbody()

def solaire(event):
    global scene, sw
    sw = 0
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    texte3d1 = text(text='Bientot disponible', align='center', depth=-0.1, color=(0.1, 0.4, 0.8), height=20, font="Times",
                    pos=(0, 15, 0))

def quitter(event):
    os._exit(0)

def reset(event):
    global scene, sw
    sw = 0
    scene.delete()
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale = (0.013, 0.013, 0.013)
    texte3d1 = text(text='N-corps', align='center', depth=-2, color=(0.7, 0.2, 0.2), height=20, font="Sans",pos=(0, -6, 0))
    while sw == 0 : rate(1)

def setspeed(event):
    global speed,slid1
    speed = slid1.GetValue()
    strspeed.SetLabel(str(speed))

def aide(event):
    global aide_thread
    class prop(wx.Frame):
        def __init__(self, titre):
            w2 = wx.Frame.__init__(self, None, 1, title=titre, size=(400, 200), style=(wx.SYSTEM_MENU | wx.CAPTION ) & ~ (wx.RESIZE_BORDER| wx.CLOSE_BOX))
            pnl1 = wx.Panel(self)
            self.Center()
            t4 = wx.StaticText(pnl1, pos=(20, 15), label="Environnement 3D", style=wx.ALIGN_CENTER)
            t4.SetFont(textboldfont)
            t1 = wx.StaticText(pnl1, pos =(20,50),label = "Rotation : Clique droit",style = wx.ALIGN_CENTER)
            t1.SetFont(textfont)
            t2 = wx.StaticText(pnl1, pos =(20,70), label = "Zoom && dezoom : Clique gauche + droit",style =wx.ALIGN_CENTER)
            t2.SetFont(textfont)
            closeBtn = wx.Button(pnl1, label="Fermer",pos=(300, 140))
            closeBtn.Bind(wx.EVT_BUTTON, self.onClose)

        def onClose(self, event):
            global aide_thread
            self.Close()
            aide_thread = 0

    class MonApp(wx.App):
        def OnInit(self):
            fen = prop("Aide")
            fen.Show(True)
            self.SetTopWindow(fen)
            return True
    def afficher():
        app = MonApp()
        app.MainLoop()
    if aide_thread == 0:
        thread.start_new_thread(afficher,())
        aide_thread = 1
        return 0
    else: return 0

def propos(event):
    global propos_thread
    class prop(wx.Frame):
        def __init__(self, titre):
            w2 = wx.Frame.__init__(self, None, 1, title=titre, size=(400, 230), style=(wx.SYSTEM_MENU | wx.CAPTION ) & ~ (wx.RESIZE_BORDER| wx.CLOSE_BOX))
            pnl1 = wx.Panel(self)
            self.Center()
            t4 = wx.StaticText(pnl1, pos=(20, 15), label="N-Corps version " + str(version), style=wx.ALIGN_CENTER)
            t4.SetFont(textboldfont)
            t1 = wx.StaticText(pnl1, pos =(20,50),label = "Developpement : Jean de Peyrecave",style = wx.ALIGN_CENTER)
            t1.SetFont(textfont)
            t2 = wx.StaticText(pnl1, pos =(20,70), label = "Design && interface : Corentin Roche",style =wx.ALIGN_CENTER)
            t2.SetFont(textfont)
            t5 = wx.StaticText(pnl1, pos=(20, 110), label="Powered by VPython and WxPython", style=wx.ALIGN_CENTER)
            t5.SetFont(textfont)
            t6 = wx.StaticText(pnl1, pos=(20, 140), label="Python : "+str(sys.version.split()[0])+"   WxPython : "+str(wx.VERSION_STRING), style=wx.ALIGN_CENTER)
            t6.SetFont(textfont)
            t3 = wx.StaticText(pnl1, pos=(20, 170), label="Copyright © 2016, La Blotisserie Inc.",style =wx.ALIGN_CENTER)
            t3.SetFont(textfont)
            png = wx.Image("logo.png", wx.BITMAP_TYPE_ANY).ConvertToBitmap()
            wx.StaticBitmap(pnl1, -1, png, (260, 15), (png.GetWidth(), png.GetHeight()))
            closeBtn = wx.Button(pnl1, label="Fermer",pos=(290, 160))
            closeBtn.Bind(wx.EVT_BUTTON, self.onClose)

        def onClose(self, event):
            global propos_thread
            self.Close()
            propos_thread = 0


    class MonApp(wx.App):
        def OnInit(self):
            fen = prop("À propos")
            fen.Show(True)
            self.SetTopWindow(fen)
            return True

    def afficher():
        app = MonApp()
        app.MainLoop()
    if propos_thread == 0:
        thread.start_new_thread(afficher,())
        propos_thread = 1
        return 0
    else : return 0

def fenetre():
    global scene,w,slid1,speed,sw,strspeed,titlefont,textboldfont,textfont

    speed = 1

    ws = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_X)
    hs = wx.SystemSettings.GetMetric(wx.SYS_SCREEN_Y)

    titlefont = wx.Font(20, wx.MODERN, wx.NORMAL, wx.BOLD)
    textboldfont = wx.Font(12, wx.SCRIPT, wx.NORMAL, wx.BOLD)
    textfont = wx.Font(10, wx.SCRIPT, wx.NORMAL, wx.NORMAL)

    w = window(title = "N-Corps "+str(version),menus=True,  width=1600, height=870, style=(wx.DEFAULT_FRAME_STYLE) & ~ (wx.RESIZE_BORDER|wx.MAXIMIZE_BOX|wx.CLOSE_BOX),x=int((ws - 1600) / 2), y=int((hs - 850)) / 2)
    pnl = w.panel
    scene = display(window=w, width=1200, height=850, background=(0.1, 0.1, 0.1))
    scene.scale = (0.013, 0.013, 0.013)

    menubar = w.menubar
    menubar.Remove(0)
    menu1 = wx.Menu()
    menu2 = wx.Menu()
    item1 = menu1.Append(-1, 'Quitter', "Quitter l'application")
    w.win.Bind(wx.EVT_MENU, quitter, item1)
    item2 = menu2.Append(-1, 'Ouvrir aide')
    w.win.Bind(wx.EVT_MENU, aide, item2)
    item3 = menu2.Append(-1, 'À propos')
    w.win.Bind(wx.EVT_MENU, propos, item3)
    menubar.Append(menu1, 'Fichier')
    menubar.Append(menu2, 'Aide')


    texte3d1 = text(text='N-corps', align='center', depth=-2, color=(0.7, 0.2, 0.2), height=20, font="Sans", pos=(0, -6, 0))

    titletext = wx.StaticText(pnl, pos=(1320, 30), label="PREREGLAGES")
    titletext.SetFont(titlefont)

    btn1 = wx.Button(pnl,label='Galaxie', pos=(1260, 100))
    btn1.Bind(wx.EVT_BUTTON,galaxie)
    btn2 = wx.Button(pnl,label='3 corps', pos=(1360, 100))
    btn2.Bind(wx.EVT_BUTTON,threebody)
    btn3 = wx.Button(pnl,label='Terre-Lune', pos=(1260, 140))
    btn3.Bind(wx.EVT_BUTTON,terrelune)
    btn4 = wx.Button(pnl, label='Saturne', pos=(1360, 140))
    btn4.Bind(wx.EVT_BUTTON, saturne)
    btn5 = wx.Button(pnl, label='Sphere', pos=(1460, 100))
    btn5.Bind(wx.EVT_BUTTON, sphere__)
    btn5 = wx.Button(pnl, label='Solaire', pos=(1460, 140))
    btn5.Bind(wx.EVT_BUTTON, solaire)
    btnq = wx.Button(pnl,label='Quitter', pos=(1497, 785))
    btnq.Bind(wx.EVT_BUTTON, quitter)
    btnr = wx.Button(pnl,label='Reset', pos=(1400, 785))
    btnr.Bind(wx.EVT_BUTTON, reset)
    btnap = wx.Button(pnl, label='À propos', pos=(1215, 785))
    btnap.Bind(wx.EVT_BUTTON, propos)

    slid1 = wx.Slider(pnl, pos=(1250, 230),size=(280,50),  minValue=1, maxValue=2000)
    slid1.Bind(wx.EVT_SCROLL, setspeed)
    strspeed = wx.StaticText(pnl, pos=(1375, 200), label=str(speed))

    text1 = wx.StaticText(pnl, pos=(1250,200), label="Vitesse de l'animation : ")
    text2 =wx.StaticText(pnl, pos=(1243, 230), label="1")
    text3 = wx.StaticText(pnl, pos=(1530, 230), label="2000")


fenetre()
