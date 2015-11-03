import ase.units as unit
import numpy as np

class fCMtoAll:
    # could me made into a utiliy class with a more general name
    # if (/when) more things are added
    def __init__(self,f_iCM=None,tau_cm=None,atoms=None,mp=3):
        self.f_iCM = f_iCM
        self.tau_cm = tau_cm
        self.atoms = atoms
        self.mp = mp
        self.nm = len(atoms)/mp
        self.f_all = None
        sq2=1/np.sqrt(2)
        oh=0.5
        mF=np.array([-sq2, oh, -oh, 0.00, sq2, sq2, sq2, oh, -oh,  \
                     -oh, -sq2, -oh, oh, -sq2, oh, -sq2, 0.0, sq2, \
                     sq2, sq2, 0.0, oh, -oh, -sq2, -oh, oh, -sq2])
        self.mF=mF.reshape((3,3,3),order='F') 
        mB=np.array([-sq2, 0.0, sq2, 0.5, sq2, 0.5, -0.5,    \
                      sq2, -0.5,-0.5, 0.5, -sq2, -sq2, -sq2, \
                      0.0, -0.5, 0.5, sq2,sq2, 0.5, -0.5,    \
                      sq2, -0.5, 0.5, 0., -sq2, -sq2])
        self.mB=mB.reshape((3,3,3),order='F') 

    def cm_coords(self):
        mp = self.mp
        mm = self.atoms
        rcm = np.zeros((self.nm,3))
        rO  = np.zeros((self.nm,3))
        rH1 = np.zeros((self.nm,3))
        rH2 = np.zeros((self.nm,3))
        for i in range(0,len(mm),mp):
            pos = mm[i:i+3].get_positions()
            cm = mm[i:i+3].get_center_of_mass()
            rcm[i/mp,:] = cm
            rO[i/mp,:]  = pos[0] - rcm[i/mp] #r3
            rH1[i/mp,:] = pos[1] - rcm[i/mp] #r1
            rH2[i/mp,:] = pos[2] - rcm[i/mp] #r2
            # or why not just keep it molwise: ThisMolInCM = pos - cm ...
        return rcm, rO, rH1, rH2 
    

    def inv6(self,r1,r2,r3): 
        # Dunno 100% what this does so cant exchange with numpy oneliner ...
        # Has been checked against original in mdUtil2.f, giving same results! 
        r1a = np.zeros((3,3))
        r2a = np.zeros((3,3))
        r3a = np.zeros((3,3))
        y = np.zeros((6,6))
        mF = self.mF
        a = r1[0]
        b = r1[1]
        c = r1[2]

        d = r2[0]
        e = r2[1]
        f = r2[2]

        g = r3[0]
        h = r3[1]
        i = r3[2]

        det = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*f*h + \
               d*f*h - a*b*i + 2*b*d*i - d*e*i - b*g*i + e*g*i + a*h*i - d*h*i

        flag = 0
        if det < 0.1 and det > -0.1:
            for iRot in range(3):
                for ii in range(3):
                    r1a[ii,iRot] = 0.0
                    r2a[ii,iRot] = 0.0
                    r3a[ii,iRot] = 0.0
                    for jj in range(3):
                        r1a[ii,iRot] = r1a[ii,iRot] + mF[ii,jj,iRot] * r1[jj]
                        r2a[ii,iRot] = r2a[ii,iRot] + mF[ii,jj,iRot] * r2[jj]
                        r3a[ii,iRot] = r3a[ii,iRot] + mF[ii,jj,iRot] * r3[jj]

            detOld = 0.0
            for iRot in range(3):
                a = r1a[0,iRot]
                b = r1a[1,iRot]
                c = r1a[2,iRot]
                
                d = r2a[0,iRot]
                e = r2a[1,iRot]
                f = r2a[2,iRot]
                
                g = r3a[0,iRot]
                h = r3a[1,iRot]
                i = r3a[2,iRot]

                det = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*f*h \
                      +d*f*h - a*b*i + 2*b*d*i - d*e*i - b*g*i + e*g*i + a*h \
                      *i -d*h*i

                if abs(det) > abs(detOld):
                   detOld = det
                   flag = iRot

            a = r1a[0,flag]
            b = r1a[1,flag]
            c = r1a[2,flag]
               
            d = r2a[0,flag]
            e = r2a[1,flag]
            f = r2a[2,flag]
               
            g = r3a[0,flag]
            h = r3a[1,flag]
            i = r3a[2,flag]
            det = detOld

      #end of if dets 0.1

        if det < 1.e-10 and det > -1.e-10:
            print "Det:"
            print det
            raise ValueError('Error in inv6(). Det = 0.')

        y[0,0] = -a*f*h + d*f*h + b*d*i - d*e*i - b*g*i + e*g*i + a*h*i  - d*h*i
        y[0,1] =  a*f*g - d*f*g - a*d*i + d*d*i
        y[0,2] = -b*d*d + a*d*e + b*d*g - a*e*g
        y[0,3] = -a*d + d*d + a*g - d*g
        y[0,4] = -b*d + d*e + b*g - e*g
        y[0,5] = -a*f + d*f + a*i - d*i

        y[1,0] = -b*c*h + c*e*h + b*b*i - b*e*i
        y[1,1] =  b*c*g - c*e*g - a*b*i + b*d*i - b*g*i + e*g*i + a*h*i - d*h*i 
        y[1,2] = -b*b*d + a*b*e + b*d*h - a*e*h
        y[1,3] = -a*b + b*d + a*h - d*h
        y[1,4] = -b*b + b*e + b*h - e*h
        y[1,5] = -b*c + c*e + b*i - e*i

        y[2,0] = -c*f*h + b*f*i + c*h*i - b*i*i
        y[2,1] = c*f*g - c*d*i - f*g*i + d*i*i
        y[2,2] = c*d*e - b*d*f - c*e*g + d*f*h + b*d*i - d*e*i + e*g*i - d*h*i 
        y[2,3] = -c*d + c*g + d*i - g*i
        y[2,4] = -b*f + f*h + b*i - h*i
        y[2,5] = -c*f + c*i + f*i - i*i

        y[3,0] = -b*c*d + c*d*e + a*b*f - b*d*f + b*c*g - c*e*g - a*b*i + b*d*i
        y[3,1] = -a*f*g + d*f*g + a*d*i - d*d*i
        y[3,2] = b*d*d - a*d*e - b*d*g + a*e*g
        y[3,3] = a*d - d*d - a*g + d*g
        y[3,4] = b*d - d*e - b*g + e*g
        y[3,5] = a*f - d*f - a*i + d*i

        y[4,0] = b*c*h - c*e*h - b*b*i + b*e*i
        y[4,1] = -b*c*d + c*d*e + a*b*f - b*d*f - a*f*h + d*f*h + b*d*i - d*e*i
        y[4,2] = b*b*d - a*b*e - b*d*h + a*e*h
        y[4,3] = a*b - b*d - a*h + d*h
        y[4,4] = b*b - b*e - b*h + e*h
        y[4,5] = b*c - c*e - b*i + e*i

        y[5,0] = c*f*h - b*f*i - c*h*i + b*i*i
        y[5,1] = -c*f*g + c*d*i + f*g*i - d*i*i
        y[5,2] = -b*c*d + a*b*f + b*c*g - a*f*h - a*b*i + b*d*i - b*g*i + a*h*i
        y[5,3] = c*d - c*g - d*i + g*i
        y[5,4] = b*f - f*h - b*i + h*i
        y[5,5] = c*f - c*i - f*i + i*i

        for jj in range(6):
            for ii in range(6):
                y[ii,jj] = y[ii,jj] / det

        return y,flag 
  
    def mol_force3(self, r1, r2, r3, ftot, tt):
        f1 = np.zeros((3))
        f2 = np.zeros((3))
        f3 = np.zeros((3))
        b = np.zeros((6))
        c = np.zeros((6))
        f1a = np.zeros((3))
        f2a = np.zeros((3))
        f3a = np.zeros((3))

        mF = self.mF
        mB = self.mB

        y,flag = self.inv6(r1,r2,r3)

        if flag > 0:
            for ii in range(3):
                f1a[ii] = 0.0
                f2a[ii] = 0.0
                for jj in range(3):
                    f1a[ii] = f1a[ii] + mF[ii,jj,flag] * ftot[jj]
                    f2a[ii] = f2a[ii] + mF[ii,jj,flag] * tt[jj]
            for ii in range(3):
                ftot[ii] = f1a[ii]
                tt[ii] = f2a[ii]
        # end of if

        for i in range(3):
            b[i] = ftot[i]
            b[i+3] = tt[i]
            
        for i in range(6):
            c[5] = 0.0
            for j in range(6):
                c[i] = c[i] + y[i,j] * b[j]

        # New forces...

        f1[0] = c[0]
        f1[1] = 0.0
        f1[2] = c[2]

        f2[0] = 0.0
        f2[1] = c[1]
        f2[2] = c[5]

        f3[0] = c[3]
        f3[1] = c[4]
        f3[2] = 0.0

        # Here we bring the forces to the original orientation of the space. 
        if flag > 0:
            for ii in range(3):
                f1a[ii] = 0.0
                f2a[ii] = 0.0
                f3a[ii] = 0.0
                for jj in range(3):
                    f1a[ii] = f1a[ii] + mB[ii,jj,flag] * f1[jj]
                    f2a[ii] = f2a[ii] + mB[ii,jj,flag] * f2[jj]
                    f3a[ii] = f3a[ii] + mB[ii,jj,flag] * f3[jj]
            for ii in range(3):
                f1[ii] = f1a[ii]
                f2[ii] = f2a[ii]
                f3[ii] = f3a[ii]
        return f1,f2,f3

    def distribute(self): 
        """ loop over each mol and get out the forces """
        mp = self.mp
        f_iCM = self.f_iCM
        tau_cm = self.tau_cm
        atoms = self.atoms

        rcm, rO, rH1, rH2 = self.cm_coords()
        fa = np.zeros((len(atoms),3))
        for mol, ftot in enumerate(f_iCM):
            tau_cm = self.tau_cm
            r1 = rH1[mol]
            r2 = rH2[mol]
            r3 = rO[mol]
            tt = tau_cm[mol]
            f1,f2,f3 = self.mol_force3(r1,r2,r3,ftot,tt)
            fa[mol*mp,:]   = r3 #OHH ... 
            fa[mol*mp+1,:] = r1
            fa[mol*mp+2,:] = r2
        return fa     
    

