import numpy as np
import numpy.linalg as la
import matplotlib.pyplot as plt
import time as tm
import calendar as cld
from pyquaternion import Quaternion
import tlm_transcription2 as tr

SurfaceX = 62.5
SurfaceY = 62.5
SurfaceZ = 62.5

dt = 0.001

kPi = 1
kSi = 1
kBi = 1
kIi = 1

#using namespace std;
#using namespace Eigen;

class Etat:
    def __init__(self):
        self.q = Quaternion(0, 1, 0, 0)
        self.b = np.zeros(3)
        self.omega = np.zeros(3)
        self.qref = Quaternion(0, 0, 1, 0)
        self.lat = 0
        self.longitude = 0
        self.momMag = np.zeros(3)
        self.intensite = np.zeros(3)
        self.controle = True

    def getJulianDate(self):
        t = cld.timegm(tm.gmtime())
        return (t / 86400.0) + 2440587.5

    def  getVecSolRef(self, jd):
        T = (jd - 2451545.)/36525
        lambdaSun = (280.4606184
                     + (36000.77005361 * T))
        MSun = (357.5277233
                + 35999.05034*T)
        lambdaE = (lambdaSun
                   + 1.914666471*sin(MSun*np.pi/180)
                   + 9.918994643*sin(2*MSun*np.pi/180))
        epsilon = (23.439291
                   - 0.0130042*T)
        s = np.array([np.cos(lambdaE*PI/180),
                      np.cos(epsilon*PI/180)*np.sin(lambdaE*PI/180),
                      np.sin(epsilon*PI/180)*np.sin(lambdaE*PI/180)])
        if (s.norm() != 0):
            s = s/la.norm(s)
        return s;

    def getCapteursSolaires(self):
        v = np.array([1, -1, 0, 0, 0, 0])
        return v

    def getVecSol(self):
        vmes = np.zeros(6)
        vecsol = np.zeros(3)
        H = np.eye(3)
        H = H/2
        M = np.array([[1, -1, 0, 0, 0 ,0],
                      [0, 0, 1, -1, 0 ,0],
                      [0, 0, 0, 0, 1,-1]])
        vmes = self.getCapteursSolaires()
        vecsol = np.dot(M, vmes)
        vecsol = np.dot(H, vecsol)
        if(la.norm(vecsol)!= 0):
            vecsol = vecsol/la.norm(vecsol)
        return vecsol

    def getChampRef(self, lat, longitude, om, i, omega, nu):
        fichier = open('tlm.3', 'r')
        x, y, z = 0, 0, 0
        #ifstream fichier("/Users/adrianvalente/Documents/etudes/psc/source1502/ADCS/igrf.txt",ios::in);
        n = round(longitude)*179+90+np.floor(lat)
        for i in range(1, n+1):
            ligne = fichier.readline()
            mag_aux = tr.traiter(tm.decoupage(ligne))
            #ATENCION
        fichier.close()
        theta = omega + nu
        m1 = np.array([[sin(i)*sin(om),cos(i)*cos(om),cos(om)],
		  [-sin(i)*cos(om),-cos(i)*sin(om),sin(om)],
		  [cos(i),sin(i),0]])
        m2 = np.array([[1,0,0],
                       [0,cos(theta),-sin(theta)],
                       [0,sin(theta),cos(theta)]])
        return np.dot(m1,np.dot(m2, mag_aux))

    def getChampM(self):
        return np.array([1, 0, 0])

    def getLat(self):
        return 0

    def getLong(self):
        return 0

    def getOmega(self):
        return np.ones(3)

    def getTLE(self, perih, om, i, nu):
        return

    def getQref(self, perih, om, i, nu):
        m1 = np.array([[np.cos(om), -np.sin(om), 0],
                       [np.sin(om), np.cos(om), 0],
                       [0,       0,      1]])
        theta=perih+nu
        m2 = np.array([[1, 0, 0],
                       [0, np.cos(i), -np.sin(i)],
                       [0, np.sin(i), np.cos(i)]])
                      
        m3 = np.array([[cos(theta), -sin(theta), 0],
                       [sin(theta), cos(theta), 0],
                       [0, 0, 1]])

        m4 = np.array([[0,-1,0],
                       [1,0,0],
                       [0,0,1]])
                      
        mtot = np.dot(m1,np.dot(m2,np.dot(m3, m4)))
    
        q = Quaternion(0)
        q[0] = sqrt(1.+mtot[0,0]+mtot[1,1]+mtot[2,2])/2.
        if (q[0]!=0):
            q[1] = (mtot[2,1]-mtot[1,2])/(4.*q[3])
            q[2] = (mtot[0,2]-mtot[2,0])/(4*q[3])
            q[3] = (mtot[1,0]-mtot[0,1])/(4*q[3])
        else: #Si la trace est nulle, il faut trouver l'axe de rotation...
            if (mtot[0,0]>0):
                q[1] = 1
                q[2] = 0
                q[3] = 0
            elif (mtot[1,1]>0):
                q[2] = 1
                q[1] = 0
                q[0] = 0
            else:
                q[3] = 1
                q[2] = 0
                q[1] = 0
        return q;

    def actualiser(self):
           #Initialisation : on recupere les donnees des capteurs et des tables
        julianDate = self.getJulianDate()
        lat = self.getLat()   #A FAIRE : recuperer les donnees du GPS
        longitude = self.getLong()
        omega = self.getOmega() #A FAIRE : recuperer vitesse de rotation du gyro
                                            #Vecteur rotation instantanée du satellite par rapport au ref geocentrique, exprimé dans le référentiel du satellite
            #TLE
            #A faire : fonction pour recuperer les donnees des TLE
        perih=0#argument du perihelie
        om=0#Longitude du point ascendant croisant le plan equatorial
        i=0#Inclinaison de l'orbite par rapport au plan equatorial
        nu=90#Anomalie (pos du satellite)
            #getTLE(&perih, &om, &i, &nu);
        #Passage en radians
        om = om*np.pi/180
        perih = perih*np.pi/180
        nu = nu*np.pi/180
        i = i*np.pi/180
        qref = self.getQref(perih, om, i, nu)
            #Champs Mags
        champRef = self.getChampRef(lat, longitude,om,i,perih,nu) #dans le ref GEOCENTRIQUE
        champRef = q.conjug(champRef)  #On le passe immediatement dans le ref SATELLITE
        champRefNorme = champRef
        if (la.norm(champRef) != 0 ):
            champRefNorme =  champRef/la.norm(champRef)
        champM = self.getChampM()#A FAIRE : Champ mesuré dans le ref SATELLITE
        champMNorme = champM
        if (la.norm(champM) != 0):
            champMNorme = champM/la.norm(champM)
            #Capteurs solaires. Les vecteurs sont directement normalises
        vecSolRef = self.getVecSolRef(julianDate)#Vecteur solaire dans le ref GEOCENTRIQUE
        vecSolRef = q.conjug(vecSolRef)
        if (la.norm(vecSolRef) != 0):
            vecSolRef = vecSolRef/la.norm(vecSolRef)
        vecSol = self.getVecSol() #Vecteur solaire mesuré dans le ref SATELLITE
        if (la.norm(vecSol) != 0):
            vecSol /= la.norm(vecSol)
            #Calcul des quaternions (dans physic.m)
        GrandOmega = kBi*np.cross(champMNorme,champRefNorme) + kSi*np.cross(vecSol, vecSolRef)  #Intermediaire de calcul (omega dans le rapport de Valentin)
        vtmp = omega - b + kPi*GrandOmega
        qtmp = Quaternion(vtmp(0),vtmp(1),vtmp(2),0.)
        q += (0.5*q.prod(qtmp) + (1-q.norm()*q.norm())*q)*dt
        b -= kIi*GrandOmega*dt   
            #Controle (l. 540)	calcul du moment genere
        dq = Quaternion()
        m = np.zeros(3)
        if (controle==true and champM.norm()!=0):
            qtmp = q.inv()
            dq = qref.prod(qtmp)
            dq13 = np.array([dq(0),dq(1),dq(2)])
            m = (-0.000048*champM.cross(omega-b-omegaref) - 0.0000003*champM.cross(dq13))/(champM.norm()*champM.norm())
    
            #Calcul intensite l.325
        intensite = np.array([m(0)/SurfaceX, m(1)/SurfaceY, m(2)/SurfaceZ])
        return intensite;

#-----------------------------------------------------------------------------

e = Etat()
