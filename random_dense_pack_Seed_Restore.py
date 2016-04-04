#Imports
"Math and other non-woo imports"
import os.path, csv, math, os
"Woo imports"
import woo, woo.dem, woo.core, woo.pack, woo.log, woo.gl, math
from woo.core import *
from woo.dem import *
from minieigen import *
from woo import qt
from woo import gl

#Is simulation being reloaded from elsewhere?
global restore
restore = False
if os.path.isfile("PackSeed100_1.bin"):
   restore = True

#Parameter Values ( Simulation )
"Frequency of both condition checking and data recording."
Frequency            = 1
"Value to determine limit of movement"
global limit
limit = 0.3
"Number of steps to record"
global steps
steps = 10000


#Parameters for Simulation Tree
"Boxsize"
global Boxsize
Boxsize  = 10	
"Poisson"
global poisson
poisson = 0.3
"Distfactor"
distfact = 1.2
"Porosity"
global poros
poros = 0.01
"Assumed weight of a sphere with r = 1"
global weight
weight=41092.0319091


#Parameter Values ( Physics & Material )
global youngs_mod
Damping = 0.3
youngs_mod = 1.0e8
density = 1000
shear_frictangle = 36
tan_Phi = math.tan( math.radians( shear_frictangle) )
KtdivKn = 1/(2*(1+poisson)) 
breakKn = 0.04 
alpha0 = 1/(2*math.pi)
alpha1 = 1/(4*math.pi)
Alpha = [ alpha0, alpha1] 
beta0 = 0.4
beta1 = 1/math.sqrt(2*math.pi)
beta2 = 1/math.sqrt(4*math.pi)
Beta                      = [ beta0, beta1, beta2]
sliding_angle = 36
Mu                        = math.tan( math.radians( sliding_angle) )


#Functions
"Used to create Scene with specified parameters"
def sceneset( damp , yng , dns , tphi , ktdkn , brkn , alp , bet , u , distfact ):
     
      "Creating the Scene, to use throughout the simulation"

      Scene = woo.master.scene = woo.core.Scene(fields=[DemField()], engines=DemField.minimalEngines(model=woo.models.ContactModelSelector(name='ice', damping=damp , distFactor = distfact )))  

      "Ice material from recieved parameters"
      global Ice
      Ice = IceMat(young=yng, density=dns, tanPhi=tphi, ktDivKn=ktdkn, breakN=brkn, alpha=alp, beta=bet, mu=u)
      "Bond properties"
      Scene.lab.cp2.bonds0 = 15
      Scene.lab.cp2.bonds1 = 0     

      return Scene
      
      
"Sets up simulation with the specidied parameters"
def simmake( boxsize , poros ):
      "Generating the packing"
      genr = MinMaxSphereGenerator( dRange = ( 2 , 2 ) )		
      pred = woo.pack.inAlignedBox( ( -3 , -boxsize/2 , -boxsize/2 ) , ( boxsize + 3 , boxsize/2 , boxsize/2 ) )
      pack = woo.pack.randomDensePack2( pred , genr , porosity = poros )
      pack.cellSize = ( 0 , 0 , 0 )
      pack.toDem( Scene , Scene.dem , mat = Ice )

      "Blocking Rotation"
      for p in Scene.dem.par:
            p.blocked = 'XYZ'
            
      "Adding the walls"
#      Scene.dem.par.add([ woo.utils.wall( position = ( 0 , 0 , 0 ) , axis = 0 , mat = Ice ) , woo.utils.wall( position = ( boxsize , 0 , 0 ) , axis = 0 , mat = Ice ) ], nodes = True|False )
#      Scene.dem.par.add([ woo.utils.wall( position = ( 0 , -boxsize/2 , 0 ) , axis = 1 , mat = Ice ) , woo.utils.wall( position = ( 0 , boxsize/2 , 0 ) , axis = 1 , mat = Ice ) ], nodes = True|False )
#      Scene.dem.par.add([ woo.utils.wall( position = ( 0 , 0 , -boxsize/2 ) , axis = 2 , mat = Ice ) , woo.utils.wall( position = ( 0 , 0 , boxsize/2 ) , axis = 2 , mat = Ice ) ], nodes = True|False )
      Scene.dem.par.add(woo.triangulated.box( dim = ((boxsize+6),boxsize,boxsize), center = ((boxsize/2),0,0), which=(0,1,1,0,1,1) ) )
      Scene.dem.par.add(woo.triangulated.box( dim = ((boxsize+6),boxsize,boxsize), center = ((boxsize/2),0,0), which=(1,0,0,1,0,0) ) )
      "Collecting nodes and returning wall"
#      Scene.dem.collectNodes()		#Added to Scene.dem.par.add above (nodes = True|False)
      return Scene.dem.par[-1]      
     
"Imposes force on specified wall/body"
def impose( par ):

      
      "This imposes a lateral movement which stretches the mass of particles"
      par.impose = woo.dem.Local6Dofs( values = ( 0.1 , 0 , 0 , 0 , 0 , 0 ) , whats = ( 1 , 0 , 0 , 0 , 0 , 0 ) )

      return par
      

"Checks conditions of the simulation for recording data"
def condition1():

      "Used for finding particles of interest"
      temp1 = 5.2
      temp2 = 4.8
      i = 1
      while i<136:			#136, 3020
            x, y, z = Scene.dem.par[i].shape.nodes[0].pos 
            if x < temp1:
                  if x > temp2:
                        if y < 1:
                              if y > -1:
                                    miny = i
                                    print (x,',', y ,',', z)
                                    print (miny)
            i = i + 1
      print (miny)
 
      if abs( particle.pos[0] - Boxsize ) > 0:
            customrunner.command = 'condition2()'
      
def condition2():

                  global dataLog
                  
                  "Initializing data log"
                  "Generating name"
                  fileName= 'settle=' + str(poros) + '_expansion_'
                  fileNum=0
                  while os.path.exists(fileName+str(fileNum)+'.csv'):fileNum+=1
                  "Creating and opening file"
                  csvfile = open(fileName+str(fileNum)+'.csv', mode='w')
                  dataLog=csv.writer(csvfile, delimiter=',',quotechar='"')
                  
                  "Writing initial row of information"
                  dataLog.writerow([ 'Force' , 'dL' , 'NPar' , 'F/A' , 'dL/L' , 'E_eff' , 'E_param' , 'dW' , 'dW/W' , 'dH' , 'dH/H' , 'nu_param' , 'nu_z' , 'nu_y' , 'nu_calc' ])

                  x, delH1i, z = Scene.dem.par[100].shape.nodes[0].pos
                  x, delH2i, z = Scene.dem.par[101].shape.nodes[0].pos
                  global Ho
                  Ho = (abs(delH1i) + abs(delH2i))

                  x, y, delW1i = Scene.dem.par[100].shape.nodes[0].pos
                  x, y, delW2i = Scene.dem.par[101].shape.nodes[0].pos
                  print (delW1i)
                  print (delW2i)
                  global Wo
                  Wo = (abs(delW1i) + abs(delW2i))
		  #print (Wo)

                  customrunner.command = 'collect(' + str(Scene.step) + ')'
                  
  
"Collects data"
def collect( stepstart ):
      
      Force = particle.f[0]
      delL = particle.pos[0] - Boxsize
      x, delH1, z = Scene.dem.par[100].shape.nodes[0].pos
      x, delH2, z = Scene.dem.par[101].shape.nodes[0].pos
      delH = (abs(delH2) + abs(delH1)) - Ho
      x, y, delW1 = Scene.dem.par[100].shape.nodes[0].pos
      x, y, delW2 = Scene.dem.par[101].shape.nodes[0].pos
      delW = (abs(delW2) + abs(delW1)) - Wo
      NPar = len(Scene.dem.par)-2
      FdivA = Force/(Boxsize*Boxsize)
      E_eff = FdivA/(delL/Boxsize)	
      EpsZ = delW/Wo
      EpsY = delH/Ho
      EpsX = delL/(Boxsize+6)	
      nu_z = EpsZ/EpsX
      nu_y = EpsY/EpsX
      
      dataLog.writerow([Force , delL , NPar , FdivA , EpsX , E_eff , youngs_mod , delW , EpsZ , delH , EpsY , poisson , nu_z , nu_y , (0.5*(nu_z+nu_y)) ])
      "Stopping Scene after recording x steps"      
      if (Scene.step - stepstart)  > steps:
            Scene.stop()
      
      
#Main Code
"Creating the Scene"
if restore:
   global Scene, particle
   Scene = woo.master.scene = Scene.load('PackSeed100_1.bin')
   idRead = open('ParticleId.txt',mode='r')
   particle = Scene.dem.par[eval(idRead.read())]
   idRead.close()
   particle = impose( Scene.dem.par[-1] )			#Imposes motion
   #fileName= 'settle=' + str(poros) + '_expansion_'
   #fileNum=0
   #while os.path.exists(fileName+str(fileNum)+'.csv'):
   #   fileNum+=1
   #csvfile = open(fileName+str(fileNum)+'.csv', mode='a')
   #dataLog=csv.writer(csvfile, delimiter=',',quotechar='"')
   customrunner = PyRunner( Frequency , 'condition1()' )
   Scene.engines += [customrunner]
else:
   Scene = sceneset( Damping , youngs_mod , density , tan_Phi , KtdivKn , breakKn , Alpha , Beta , Mu , distfact )
   "Creating the simulation and imposing forces"
   global particle
   particle = impose( simmake( Boxsize , poros ) )		#Creates pack and imposes motion
   Scene.save("PackSeed100_1.bin")
   idWrite = open('ParticleId.txt',mode='w')
   particleIdWrite = idWrite.write(str(particle.id))
   idWrite.close()
   "Adding custom runner"
   customrunner = PyRunner( Frequency , 'condition1()' )
   Scene.engines += [customrunner]
