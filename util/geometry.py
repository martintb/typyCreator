from ..cy.lib import maxDist,minDist
import numpy as np
def rotateTo(positions,rot_from,rot_to,transPos=0):
  rot_from = np.array(rot_from)
  rot_to = np.array(rot_to)
  if all(np.array(rot_from) == np.array(rot_to)):
    R = np.eye(3)
  else:
    #rot_from: 3 element vector defining vector in pos used for the rotation
    #rot_to: 3 element vector defining the final location of rot_from after rotation
    rot_to_norm=(rot_to[0]**2.0+rot_to[1]**2.0+rot_to[2]**2.0)**(0.5)
    rot_to = [rot_to[0]/rot_to_norm,rot_to[1]/rot_to_norm,rot_to[2]/rot_to_norm]
    rot_from_norm=(rot_from[0]**2.0+rot_from[1]**2.0+rot_from[2]**2.0)**(0.5)
    rot_from = [rot_from[0]/rot_from_norm,rot_from[1]/rot_from_norm,rot_from[2]/rot_from_norm]
    
    #calculate rotation axis and rotation angle
    rot_axis = np.cross(rot_from,rot_to)/np.linalg.norm(np.cross(rot_from,rot_to))
    rot_angle = np.arccos(np.dot(rot_to,rot_from))

    #calculate quanterion elements
    a = np.cos(rot_angle/2.)
    b = rot_axis[0]*np.sin(rot_angle/2.)
    c = rot_axis[1]*np.sin(rot_angle/2.)
    d = rot_axis[2]*np.sin(rot_angle/2.)

    #construct rotation matrix
    R1 = np.array([a**2+b**2-c**2-d**2,2*(b*c-a*d),2*(b*d+a*c)])
    R2 = np.array([2*(b*c+a*d),a**2+c**2-b**2-d**2,2*(c*d-a*b)])
    R3 = np.array([2*(b*d-a*c),2*(c*d+a*b),a**2+d**2-b**2-c**2])
    R = np.vstack((R1,R2,R3))
  
  #calculate rotated positions
  rot_pos = np.dot(positions, R.T)

  transVec=transPos-rot_pos[0]
  return (rot_pos+transVec)

def spiralPos(natoms, bondLength=1.4):
  d=bondLength
  a=d; b=0.5*d; c=b*(3.**(0.5))
  baseHex = [[a,0.,0.],[b,c,0.],[-b,c,0.],[-a,0.,0.],[-b,-c,0.],[b,-c,0.]]

  stepDir=0
  pos=np.array([[0.,0.,0.]])
  for i in range(natoms-1):
    r1=np.array([pos[i]])
    hexDir=(stepDir+2)%6
    r2=np.array([pos[i]+baseHex[hexDir]])
    hexDir=(stepDir+1)%6
    r3=np.array([pos[i]+baseHex[hexDir]])
    hexDir=(stepDir)%6
    r4=np.array([pos[i]+baseHex[hexDir]])
    if minDist(pos,r2,0,False)[2]>1e-9:
      pos = np.append(pos,r2,axis=0)
      stepDir = (stepDir+2)%6
    elif minDist(pos,r3,0,False)[2]>1e-9:
      pos = np.append(pos,r3,axis=0)
      stepDir=(stepDir+1)%6
    else:
      pos = np.append(pos,r4,axis=0)
  return pos

def linearPos(natoms, bondLength=1.4):
  return np.array([[0.,0.,float(i)*bondLength] for i in range(natoms)])

def compactPos(natoms, bondLength=1.0, lead=5):
  d=bondLength
  a=d; b=0.5*d; c=b*(3.**(0.5))
  baseHex = np.array([[a,0.,0.],[b,c,0.],[-b,c,0.],[-a,0.,0.],[-b,-c,0.],[b,-c,0.]])
  HexNum=0
  
  pos=[]
  pos = [[a,0.,float(i)*bondLength] for i in range(lead)]
  newPos=np.array([pos[-1]+baseHex[5]])
  pos = np.append(pos,newPos,axis=0)
  pos=np.array(pos)
  stepDir=0
  while len(pos)<=natoms:
    newPos=np.array([pos[-1]+baseHex[(stepDir+1)%6]])
    if minDist(pos,newPos,0,False)[2]>1e-9:
      pos=np.append(pos,newPos,axis=0)
      stepDir = (stepDir+1)%6
    else:
      newPos=np.array([pos[-1]+[0,0,d]])
      pos=np.append(pos,newPos,axis=0)

      newPos=np.array([[a,0,pos[-1,2]]])
      pos=np.append(pos,newPos,axis=0)

      newPos=np.array([pos[-1]+baseHex[stepDir]])
      pos = np.append(pos,newPos,axis=0)
      stepDir=(stepDir-5)%6

  return pos[:natoms]

def compactPos2(natoms, bondLength=1.4, lead=6, layerSize=12):
  pos = [[0.,0.,float(i)*bondLength] for i in range(lead)]
  nLayers = int(np.ceil((natoms-lead)/float(layerSize)))
  flip=True
  z=pos[-1][2]+bondLength
  for i in range(nLayers):
    newPos = getSpiralPos(natoms=layerSize,bondLength=bondLength)
    newPos = [[a[0],a[1],a[2]+z] for a in newPos]
    if flip:
      pos.extend(newPos[::1])
    else:
      pos.extend(newPos[::-1])
    flip = not flip
    z+=bondLength
  return np.array(pos[:natoms])

def compactPos3(natoms, bondLength=1.4, growth=4):
  pos = [[0.,0.,0.]]
  flip=True
  z=pos[-1][2]+bondLength
  a=range(1,natoms*growth,growth) 
  b = [i for j in zip(a,a,a) for i in j]
  for i in range(natoms):
    newPos = getSpiralPos(natoms=b[i],bondLength=bondLength)
    newPos = [[a[0],a[1],a[2]+z] for a in newPos]
    if flip:
      pos.extend(newPos[::1])
    else:
      pos.extend(newPos[::-1])
    flip = not flip
    z+=bondLength
  return np.array(pos[:natoms])
