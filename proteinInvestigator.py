'''
Created on Jan 21, 2014

@author: ola
'''


from Bio.PDB import PDBIO, PDBParser
from Bio.PDB.PDBExceptions import PDBException
from Bio.PDB.Polypeptide import PPBuilder, Polypeptide
from Bio.PDB.Vector import calc_dihedral
from argparse import ArgumentParser
from math import pi
import argparse
import sys



parser=PDBParser()

def argp():
    ''' parse arguments given by user in command line'''
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-pdb", dest="pdbFile", required=True, help="pdbFile", type=str)
    args = parser.parse_args()
    return args

def unifyAngle(angle): #angle in degrees
    i=1
    while angle <-180 or angle>180:
        if angle >180:
            angle=angle-180*i
        if angle <-180:
            angle=angle+180*i
    return angle

class _PPBuilder2(PPBuilder):
    def __init__(self):
        PPBuilder.__init__(self)
   
    def build_peptides(self, entity, aa_only=1):
        """Build and return a list of Polypeptide objects.
       
         @param entity: polypeptides are searched for in this object
         @type entity: L{Structure}, L{Model} or L{Chain}
       
         @param aa_only: if 1, the residue needs to be a standard AA
         @type aa_only: int
         """
        is_connected=self._is_connected
        accept=self._accept
        level=entity.get_level()
        # Decide which entity we are dealing with
        if level=="S":
            model=entity[0]
            chain_list=model.get_list()
        elif level=="M":
            chain_list=entity.get_list()
        elif level=="C":
            chain_list=[entity]
        else:
            raise PDBException("Entity should be Structure, Model or Chain.")
        pp_list=[]
        for chain in chain_list:
            chain_it=iter(chain)
            try:
                prev_res = next(chain_it)
                while not accept(prev_res, aa_only):
                    prev_res = next(chain_it)
            except StopIteration:
                #No interesting residues at all in this chain
                continue
            pp=None
            for next_res in chain_it:
                if accept(prev_res, aa_only) and accept(next_res, aa_only) and is_connected(prev_res, next_res):
                    if pp is None:
                        pp=Polypeptide2()
                        pp.append(prev_res)
                        pp_list.append(pp)
                    pp.append(next_res)
                else:
                    #Either too far apart, or one of the residues is unwanted.
                    #End the current peptide
                    pp=None
                prev_res=next_res
        return pp_list
       
       
class Polypeptide2(Polypeptide):
    def __init__(self):
        Polypeptide.__init__(self)
       
    def get_phi_psi_list(self):
        """Return the list of phi/psi dihedral angles."""
        ppl=[]
        pplDict={}
        lng=len(self)
        for i in range(0, lng):
            res=self[i]
            try:
                n=res['N'].get_vector()
                ca=res['CA'].get_vector()
                c=res['C'].get_vector()
            except:
                # Some atoms are missing
                # Phi/Psi cannot be calculated for this residue
                ppl.append((None, None))
                pplDict.update({res:(None, None)})
                res.xtra["PHI"]=None
                res.xtra["PSI"]=None
                continue
            # Phi
            if i>0:

                rp=self[i-1]
                try:
                    cp=rp['C'].get_vector()
                    phi=calc_dihedral(cp, n, ca, c)
                except:
                    phi=None
            else:
                # No phi for residue 0!
                phi=None
            # Psi
            if i<(lng-1):
                rn=self[i+1]
                try:
                    nn=rn['N'].get_vector()
                    psi=calc_dihedral(n, ca, c, nn)
                except:
                    psi=None
            else:
                # No psi for last residue!
                psi=None
            ppl.append((phi, psi))
            pplDict.update({res:(phi, psi)})
            # Add Phi/Psi to xtra dict of residue
            res.xtra["PHI"]=phi
            res.xtra["PSI"]=psi
        return pplDict 

class protein_investigation():
    def __init__(self,resNr,chainName,pdbFile):
        self.resNr=resNr
        self.chainName=chainName
        self.pdbFile=pdbFile
        self.result =False
        #structure=self.structure
        #self.structure_file = structure_file
        self.structure=parser.get_structure(pdbFile,pdbFile)
   
    def ifChainInStructure(self):
        result=False
        self.chain=[chain for chain in self.structure[0].child_list if chain.id==self.chainName]
        if self.chain:
            result= True
        return result
   
       
    def if_structure_in_file(self):
        #try:
        structure=parser.get_structure(self.pdbFile,self.pdbFile)
        if len([a for a in structure.get_residues()])!=0:
            self.result= True
           
        #print self.result
        return self.result   
   
    def findResInStructure(self):
        self.residue=False
     
        b=False
        if (self.if_structure_in_file() and self.ifChainInStructure()) ==True:
       
            structure=parser.get_structure(self.pdbFile,self.pdbFile)
       
            try:
                self.residue =[a for a in structure.get_residues() if a.id[1]==int(self.resNr)][0]
                b=True
            except:
                pass
      
        if b==False:
            print 'can not find residue in structure'
            sys.exit()
        else: 
            #print b 
            #print self.residue
            return self.residue
       
    def CalcDihedral(self,atom1,atom2,atom3,atom4):
        if self.findResInStructure():
            vector1=atom1.get_vector()
            vector2=atom2.get_vector()
            vector3=atom3.get_vector()
            vector4=atom4.get_vector()
            #print vector1, vector2,vector3,vector4
            angle=calc_dihedral(vector1, vector2, vector3, vector4)*180/pi
            #print angle
            return angle
        else:
            return False
           
 
    def defineSs(self,polipeptide): #define secondary structure
        phiPsidict=polipeptide.get_phi_psi_list()
       
       
        dictSs={}
        for (residue,phiPsi) in [(k,v) for (k,v) in phiPsidict.iteritems()]:
            (phi,psi)=phiPsi
            ss='U'
            if phi !=None:
                #phi=phi*180
                phi=unifyAngle(phi*180)
           
            if psi !=None:
                #psi=psi*180
                psi=unifyAngle(psi*180)
               
            if -145<phi<-135 and 130<psi<140:
                ss='b_anti'
                #print ss
            if  -125<phi<-115 and 110<psi<120:
                ss='b_par'
                #print residue.resname, residue.id[1], ss
            if -180<phi<-100 and 100<psi<180:
                ss='betha'
                #print residue.resname,residue.id[1], ss
            if -135<phi<-45 and-50<psi<40:
                ss='alpha'
                #print residue.resname, residue.id[1],ss
            #else:
                #ss='U'
                #print residue.resname,residue.id[1], ss
           
            print residue.resname,residue.id[1],ss,phi,psi
           
        pass
           
ppb=_PPBuilder2()
def main():
   
   
    #print unifyAngle(1200)
    #sys.exit()
    pdbFile=argp().pdbFile
    structure=parser.get_structure(pdbFile,pdbFile)
    chainA=[a for a in structure.get_chains() if a.id=='A'][0]
    polipeptide=ppb.build_peptides(chainA,1)[0]
    #print polipeptide
    #sys.exit()
    resNr=238
    chainName='A'
    #peptide = [res for res in structure.get_residues()][0:4]
   
    #print peptide
    #[atom1,atom2,atom3,atom4]=[a for a in structure.get_atoms()][0:4]
    #print [atom1,atom2,atom3,atom4]
    #print  [atom1.get_coord(),atom2.get_coord(),atom3.get_coord(),atom4.get_coord()] 
    #residue=parser.get_structure(pdbFile,pdbFile)
    PI= protein_investigation(resNr,chainName,pdbFile)
    #PI.CalcDihedral(atom1, atom2, atom3, atom4)
    PI.defineSs(polipeptide)
   
   

if __name__ == '__main__':
    main()


if __name__ == '__main__':
    pass