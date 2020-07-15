from __future__ import print_function, absolute_import, division
import _rsumf
import f90wrap.runtime
import logging

class Real_Space_Electrostatic_Sum(f90wrap.runtime.FortranModule):
    """
    Module real_space_electrostatic_sum
    
    
    Defined at ../source/real_space_electrostatic_sum.f90 lines 22-322
    
    """
    @staticmethod
    def energy(a1, a2, a3, n, rx, ry, rz, z, rc, rd):
        """
        e = energy(a1, a2, a3, n, rx, ry, rz, z, rc, rd)
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 lines 29-123
        
        Parameters
        ----------
        a1 : float array
        a2 : float array
        a3 : float array
        n : int
        rx : float array
        ry : float array
        rz : float array
        z : float array
        rc : float
        rd : float
        
        Returns
        -------
        e : float
        
        """
        e = _rsumf.f90wrap_energy(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, \
            rc=rc, rd=rd)
        return e
    
    @staticmethod
    def add_energy_corr_term(zi, qi, rho, rd):
        """
        corr_term = add_energy_corr_term(zi, qi, rho, rd)
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 lines 125-139
        
        Parameters
        ----------
        zi : float
        qi : float
        rho : float
        rd : float
        
        Returns
        -------
        corr_term : float
        
        """
        corr_term = _rsumf.f90wrap_add_energy_corr_term(zi=zi, qi=qi, rho=rho, rd=rd)
        return corr_term
    
    @staticmethod
    def force(a1, a2, a3, n, rx, ry, rz, z, rc, rd, fx, fy, fz):
        """
        force(a1, a2, a3, n, rx, ry, rz, z, rc, rd, fx, fy, fz)
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 lines 141-208
        
        Parameters
        ----------
        a1 : float array
        a2 : float array
        a3 : float array
        n : int
        rx : float array
        ry : float array
        rz : float array
        z : float array
        rc : float
        rd : float
        fx : float array
        fy : float array
        fz : float array
        
        """
        _rsumf.f90wrap_force(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, rc=rc, \
            rd=rd, fx=fx, fy=fy, fz=fz)
    
    @staticmethod
    def stress(a1, a2, a3, n, rx, ry, rz, z, rc, rd, s):
        """
        stress(a1, a2, a3, n, rx, ry, rz, z, rc, rd, s)
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 lines 210-294
        
        Parameters
        ----------
        a1 : float array
        a2 : float array
        a3 : float array
        n : int
        rx : float array
        ry : float array
        rz : float array
        z : float array
        rc : float
        rd : float
        s : float array
        
        """
        _rsumf.f90wrap_stress(a1=a1, a2=a2, a3=a3, n=n, rx=rx, ry=ry, rz=rz, z=z, rc=rc, \
            rd=rd, s=s)
    
    @staticmethod
    def invert_3x3(a, b):
        """
        invert_3x3(a, b)
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 lines 296-322
        
        Parameters
        ----------
        a : float array
        b : float array
        
        """
        _rsumf.f90wrap_invert_3x3(a=a, b=b)
    
    @property
    def dp(self):
        """
        Element dp ftype=integer pytype=int
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 line 24
        
        """
        return _rsumf.f90wrap_real_space_electrostatic_sum__get__dp()
    
    @property
    def pi(self):
        """
        Element pi ftype=real(dp) pytype=float
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 line 25
        
        """
        return _rsumf.f90wrap_real_space_electrostatic_sum__get__pi()
    
    @property
    def sqrt_pi(self):
        """
        Element sqrt_pi ftype=real(dp) pytype=float
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 line 26
        
        """
        return _rsumf.f90wrap_real_space_electrostatic_sum__get__sqrt_pi()
    
    @property
    def one_third(self):
        """
        Element one_third ftype=real(dp) pytype=float
        
        
        Defined at ../source/real_space_electrostatic_sum.f90 line 27
        
        """
        return _rsumf.f90wrap_real_space_electrostatic_sum__get__one_third()
    
    def __str__(self):
        ret = ['<real_space_electrostatic_sum>{\n']
        ret.append('    dp : ')
        ret.append(repr(self.dp))
        ret.append(',\n    pi : ')
        ret.append(repr(self.pi))
        ret.append(',\n    sqrt_pi : ')
        ret.append(repr(self.sqrt_pi))
        ret.append(',\n    one_third : ')
        ret.append(repr(self.one_third))
        ret.append('}')
        return ''.join(ret)
    
    _dt_array_initialisers = []
    

real_space_electrostatic_sum = Real_Space_Electrostatic_Sum()

