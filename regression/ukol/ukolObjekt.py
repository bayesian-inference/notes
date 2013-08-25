import numpy as np
#import matplotlib.pyplot as plt 


from scipy import random, linalg

MATRIXSIZE = 50
MAXONESIZE = 17
RANDOMSPREAD=30
RANDOMSIZE=30
REALCOUNT = 1000
EXAMPLESIZE = 5
VELIKOST_KOLEM_JEDNOTKOVE = 1
PRIOR_MEAN_NONZERO=1000
RANDOMDATA_DEV=8
RANDOMDATA_KONST=30
PRIOR_B=10

assert(MATRIXSIZE>MAXONESIZE*2)

class TypMnozinaFci:
    SINY_POSUN = 2
    SINY_FREKVPLUSPOSUN = 3
    POLYNOMY=4

def fce(typMnozinaFci, n, x):
    assert(len(n.shape)==2)
    assert(typMnozinaFci<=4)
    assert(typMnozinaFci>=2)
    if (typMnozinaFci== TypMnozinaFci.SINY_POSUN):
        return fce_siny_posun(n,x)
    if (typMnozinaFci == TypMnozinaFci.SINY_FREKVPLUSPOSUN):
        return fce_siny_frekvposun(n,x)
    if (typMnozinaFci==TypMnozinaFci.POLYNOMY):
        return fce_polynom(n,x)
    assert(False)

def fce_siny_posun(n,x):
    half=MATRIXSIZE/2
    #assert(half==25)
    
    assert(len(n.shape) == 2)

    jenula = (n==0).squeeze()

    res = np.sin(x+(np.pi/half)*n)
    if jenula.shape==(): jenula = np.array([jenula])
    res[jenula]=1
    assert(res.shape[0]==n.shape[0])
    assert(res.shape[1]==x.shape[0])
    return res


def fce_siny_frekvposun(n,x):
    #cast = MATRIXSIZE*2/5
    #cast=20
    cast=MATRIXSIZE/3

    posuny = n%cast
    verzeFrekv = n/cast
    verzeFrekv = verzeFrekv+1
    #posuny=float(posuny)
    #print div
    #print posuny
    #assert(false)
    jenula = (n==0).squeeze()
    if jenula.shape==(): jenula = np.array([jenula])
    res = np.sin((x+(np.pi/10.0)*posuny)*verzeFrekv)    
    res[jenula]=1
    #print res
    #assert(false)
    assert(res.shape[0]==n.shape[0])
    assert(res.shape[1]==x.shape[0])
    return res

def fce_polynom(n,x):
    res=x**(n/10.0)
    if (len(res.shape)==1):
        res = np.array([res])
    #print "res - " 
    #print res.shape
    #print "n - "
    #print n.shape
    #print "x - " 
    #print x.shape
    assert(res.shape[0]==n.shape[0])
    assert(res.shape[1]==x.shape[0])
    return res
    

class TypPriorMean:
     VSE_NULA = 1
     VSE_TIS = 2
     ZACATEK = 3
     KONEC=4

def getPriorMean(typPriorMean, matrixSize, maxOnSize):
    assert(typPriorMean<=4)
    #assert(matrixSize>=5)
    #assert(maxOnSize>=5)
    if (typPriorMean == TypPriorMean.VSE_NULA or typPriorMean == TypPriorMean.KONEC):
        prvniPul = [0]*(matrixSize/2)
    else :
        prvniPul = [PRIOR_MEAN_NONZERO]*maxOnSize+[0]*(matrixSize/2-maxOnSize)
    
    if (typPriorMean == TypPriorMean.VSE_NULA or typPriorMean == TypPriorMean.ZACATEK):
        druhaPul = [0]*(matrixSize/2)
    else :
        druhaPul = [0]*(matrixSize/2-maxOnSize)+[PRIOR_MEAN_NONZERO]*maxOnSize
    celek = prvniPul+druhaPul
    assert(len(prvniPul)==matrixSize/2)
    assert(len(druhaPul)==matrixSize/2)
    assert(len(celek)==matrixSize)
    return celek

class TypDat:
    KOLEM_NULY=0
    KOLEM_KONSTANTY=1   
    KONST_KRAT_SIN_1=2
    KONST_KRAT_SIN_2=3

def trenovaciData(typDat, kolik):
    assert(typDat >= 0 and typDat <= 3)
    if typDat==TypDat.KOLEM_NULY:
        res = nahodnaData_kolemNuly(kolik)
    elif typDat == TypDat.KOLEM_KONSTANTY:
        res = nahodnaData_konst(kolik)
    elif typDat == TypDat.KONST_KRAT_SIN_1:
        res = nahodnaData_sin(kolik,1.0)
    elif typDat == TypDat.KONST_KRAT_SIN_2:
        res = nahodnaData_sin(kolik,2.0)
    assert(res.shape==(kolik,))
    return res

def nahodnaData_kolemNuly(kolik):
    return np.random.normal(0,RANDOMDATA_DEV,kolik)

def nahodnaData_konst(kolik):
    return np.random.normal(RANDOMDATA_KONST, RANDOMDATA_DEV, kolik)

def nahodnaData_sin(kolik,i):
    randomdata = np.random.normal(1,RANDOMDATA_DEV,kolik)
    randX = np.linspace(0, RANDOMSPREAD, RANDOMSIZE)
   # randomdata = np.sin(i*1.0*randX)
    randomdata = abs(randomdata) *np.sin(i*1.0*randX)
    return randomdata

class TypPriorCov:
   JEDNOTKOVA = 1
   #JEDNICKY = 2
   KOLEM_JEDNOTKOVE = 3
   JEDNICKY_VLEVO_NAHORE = 4
   JEDNICKY_VPRAVO_DOLE = 5
   RANDOM= 6

def getPriorCov(typPriorCov, matrixSize, maxOnSize, koef):
        assert(typPriorCov<=6)
        #assert(matrixSize>=5)
        #assert(maxOnSize>=5)
        half=matrixSize/2
        if (typPriorCov == TypPriorCov.JEDNOTKOVA):
            priorCov = np.identity(matrixSize)
        #elif (typPriorCov == TypPriorCov.JEDNICKY):
        #    priorCov = np.ones((matrixSize, matrixSize))
        elif (typPriorCov == TypPriorCov.KOLEM_JEDNOTKOVE):
            priorCov = almost_diag(matrixSize,1)+almost_diag(matrixSize,2)+almost_diag(matrixSize,3)
            priorCov = priorCov.dot(priorCov.transpose())
        elif (typPriorCov == TypPriorCov.JEDNICKY_VLEVO_NAHORE
             or typPriorCov == TypPriorCov.JEDNICKY_VPRAVO_DOLE):
            jednicky = np.ones((maxOnSize,maxOnSize))
            nulyPravo = np.zeros((maxOnSize, matrixSize-maxOnSize))
            nulyDole = np.zeros((matrixSize-maxOnSize,matrixSize))
            if (typPriorCov==TypPriorCov.JEDNICKY_VLEVO_NAHORE):
                dorovnano = np.hstack((jednicky, nulyPravo))
                priorCov = np.vstack((dorovnano, nulyDole)) 
                for i in range(0,maxOnSize): priorCov[i,i]+=1
                for i in range(maxOnSize,matrixSize): priorCov[i,i]=10e-30
            else:
                dorovnano = np.hstack((nulyPravo,jednicky))
                priorCov = np.vstack((nulyDole,dorovnano)) 
                for i in range(0,maxOnSize): 
                    j=matrixSize-i-1
                    priorCov[j,j]+=1
                for i in range(maxOnSize,matrixSize): 
                    j=matrixSize-i-1
                    priorCov[j,j]=10e-30
                
                #print priorCov.shape

        elif (typPriorCov == TypPriorCov.JEDNICKY_VPRAVO_DOLE):
            priorCov = np.vstack((np.zeros((maxOnSize,matrixSize)),np.hstack((np.zeros((maxOnSize,maxOnSize-maxOnSize)),np.ones((maxOnSize,maxOnSize)))))) 
            for i in range(0,matrixSize): priorCov=1
        elif (typPriorCov == TypPriorCov.RANDOM):
            priorCov = random.rand(matrixSize,matrixSize)
            priorCov = priorCov.dot(priorCov.transpose())

        assert (np.linalg.inv(priorCov).shape == (matrixSize, matrixSize))
        return priorCov*koef

def almost_diag(size,far):
    res = np.zeros((size,size))
    for i in range(0,size):
        for j in range(0,size):
            if (abs(i-j)<=far):
                res[i,j]=VELIKOST_KOLEM_JEDNOTKOVE
    return res

class UkolInference:
    def __init__(self, typPriorCov, covKoef, typPriorMean, typMnozinaFci, typDat, sigma=0, gammaA=0):
        assert(sigma != 0 or gammaA != 0)

        self.priorCov = getPriorCov(typPriorCov, MATRIXSIZE, MAXONESIZE, covKoef)
        self.priorMean = getPriorMean(typPriorMean, MATRIXSIZE, MAXONESIZE)
        self.typMnozinaFci = typMnozinaFci

        self.linX = np.linspace(0, RANDOMSPREAD, REALCOUNT)    
        self.randX = np.linspace(0, RANDOMSPREAD, RANDOMSIZE)
        self.typDat=typDat
        self.randY=self.trenovaciData(RANDOMSIZE)
        Phi=self.phi(self.randX).transpose()
        if (sigma != 0):
            self.jeNormalGamma=0    
            nasobek= (sigma**(-1))*Phi.transpose().dot(Phi);
            self.posteriorCov = np.linalg.inv(np.linalg.inv(self.priorCov)+nasobek) 
            nasobekA=np.linalg.inv(self.priorCov).dot(self.priorMean)
            nasobekB=(sigma**(-1))*Phi.transpose().dot(self.randY)
            self.posteriorMean = self.posteriorCov.dot(nasobekA+nasobekB)   
        else:
            self.jeNormalGamma=1
            nasobekA= np.linalg.inv(self.priorCov)+Phi.transpose().dot(Phi);
            nasobekB= np.linalg.inv(self.priorCov).dot(self.priorMean)+Phi.transpose().dot(self.randY);
            self.posteriorCov =np.linalg.inv(nasobekA) 
            self.posteriorMean = self.posteriorCov.dot(nasobekB)
            prv = np.array(self.priorMean).transpose().dot(np.linalg.inv(self.priorCov)).dot(self.priorMean)
            dru = self.randY.transpose().dot(self.randY)
            tri = self.posteriorMean.transpose().dot(np.linalg.inv(self.posteriorCov)).dot(self.posteriorMean)
            self.posteriorB=PRIOR_B+0.5*( prv + dru - tri)
            self.posteriorA = gammaA + 0.5*(self.randY.shape[0])
        assert(self.posteriorMean.shape==np.array(self.priorMean).shape)
        assert(self.posteriorCov.shape==self.priorCov.shape)
        

    def fce(self, n,x):
        assert(len(n.shape) == 2)
        res= fce(self.typMnozinaFci,n,x)
        assert(res.shape==(n.shape[0],x.shape[0]))
        return res

    def phi(self,x):
        phi = self.fce(np.array([range(0,MATRIXSIZE)]).transpose(),x) 
        return phi

    def trenovaciData(self,kolik):
        res= trenovaciData(self.typDat, kolik)
        assert(res.shape==(kolik,))
        return res
        
    def fce_s_vahami(self,vahy,x):
        phi = self.phi(x)
        res = vahy.dot(phi)
        return res.squeeze()

    def maluj_base_fce(self, kam):
        #for i in range(1,20):
        for i in range(0,MATRIXSIZE):
            #print self.linX.shape
            #print self.fce(np.array([[i]]),self.linX).shape
            
            kam.plot(self.linX, self.fce(np.array([[i]]),self.linX).transpose())
        
    def maluj_nahodne_priory(self,kam):
        for i in range(0,EXAMPLESIZE):
            randomly=np.random.multivariate_normal(self.priorMean, self.priorCov,1)
            #print randomly
            kam.plot(self.linX, self.fce_s_vahami(randomly,self.linX))
        #assert(False)
    
    def maluj_nahodna_data(self,kam):
        kam.plot(self.randX, self.randY, 'bo')
    
    def maluj_nahodne_posteriory(self, kam):
        for i in range(0,EXAMPLESIZE):
            if (self.jeNormalGamma==0) :
                randomly=np.random.multivariate_normal(self.posteriorMean, self.posteriorCov,1)
            
            else:
                print 1.0/(self.posteriorA)
                randomly_lambda = np.random.gamma(self.posteriorB, 1.0/(self.posteriorA),1)        
                randomly_lambda = randomly_lambda[0]
                randomly=np.random.multivariate_normal(self.posteriorMean, (randomly_lambda**(-1) )* self.posteriorCov,1)
            kam.plot(self.linX, self.fce_s_vahami(randomly,self.linX))
        
            
#ukol = UkolInference(TypPriorCov.KOLEM_JEDNOTKOVE, 1000, TypPriorMean.KONEC, TypMnozinaFci.POLYNOMY, TypDat.KOLEM_NULY, sigma=10000)
##ukol = UkolInference(TypPriorCov.RANDOM, TypPriorMean.VSE_NULA, TypMnozinaFci.POLYNOMY, TypDat.KONST_KRAT_SIN_1, sigma=10)
##ukol = UkolInference(TypPriorCov.JEDNOTKOVA, TypPriorMean.VSE_NULA, TypMnozinaFci.POLYNOMY, TypDat.KOLEM_NULY, sigma=10)
#
##ukol.maluj_nahodne_priory(plt)
##ukol.maluj_base_fce(plt)
#ukol.maluj_nahodne_posteriory(plt)
##ukol.maluj_nahodne_priory(plt)
#ukol.maluj_nahodna_data(plt)

#plt.show()
