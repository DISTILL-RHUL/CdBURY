import numpy
from enum import Enum,auto

from typing import Self

class CadmiumFluxes:
    def __init__(self,input,oxyhydroxide,carbonate,clay,organic,sulphide):
        self.input = input
        self.oxyhydroxide = oxyhydroxide
        self.carbonate = carbonate
        self.clay = clay
        self.organic = organic
        self.sulphide = sulphide
    def __add__(self, other):
        if isinstance(other, CadmiumFluxes):
            return CadmiumFluxes(self.input+other.input,
                                self.oxyhydroxide+other.oxyhydroxide,
                                self.carbonate+other.carbonate,
                                self.clay+other.clay,
                                self.organic+other.organic,
                                self.sulphide+other.sulphide)
        raise TypeError("Unsupported operand type(s) for +: 'CadmiumFluxes' and '{}'".format(type(other).__name__))
    def __mul__(self,other):
        if isinstance(other, (int, float)):
            return CadmiumFluxes(self.input*other,
                                self.oxyhydroxide*other,
                                self.carbonate*other,
                                self.clay*other,
                                self.organic*other,
                                self.sulphide*other)
        raise TypeError("Unsupported operand type(s) for *: 'CadmiumFluxes' and '{}'".format(type(other).__name__))

class CadmiumMassIsotope:
    def __init__(self,mass,isotope_delta):
        self.standard = 2.30416
        self.mass = mass
        self.isotope_delta = isotope_delta
    @staticmethod
    def from_masses(mass_114,mass_110):
        standard = 2.30416

        mass = mass_114 + mass_110
        isotope_delta = (((mass_114/mass_110)/(standard))-1.0)*1000.0
        return CadmiumMassIsotope(mass,isotope_delta)
    @property
    def isotope_mass(self):
        return self.mass * self.isotope_delta
    @property
    def isotope_ratio(self):
        return ((self.isotope_delta/1000)+1)*self.standard
    @property
    def isotope_114_fraction(self):
        return self.isotope_ratio/(1+self.isotope_ratio)
    @property
    def isotope_110_fraction(self):
        return 1-self.isotope_114_fraction
    @property
    def mass_114(self):
        return self.isotope_114_fraction*self.mass
    @property
    def mass_110(self):
        return self.isotope_110_fraction*self.mass
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            total_mass = self.mass
            if total_mass>0:
                return CadmiumMassIsotope.from_masses(self.mass_114*other,self.mass_110*other)
            else:
                return CadmiumMassIsotope(0.0,self.isotope_delta)
        raise TypeError("Unsupported operand type(s) for *: 'CadmiumMassIsotope' and '{}'".format(type(other).__name__))
    def __add__(self, other):
        if isinstance(other, CadmiumMassIsotope):
            total_mass = self.mass + other.mass
            if total_mass>0:
                return CadmiumMassIsotope.from_masses(self.mass_114+other.mass_114,self.mass_110+other.mass_110)
            else:
                return CadmiumMassIsotope(0.0,self.isotope_delta)
        raise TypeError("Unsupported operand type(s) for +: 'CadmiumMassIsotope' and '{}'".format(type(other).__name__))

class Time:
    def __init__(self,start,step,stop,units="yr"):
        self.start = start
        self.step = step
        self.stop = stop
        self.reversed = start>stop
        self.units = self.translate_units(units)
    @property
    def length(self):
        return abs(int((self.stop - self.start)/self.step))
    @property
    def padded_length(self):
        return len(self.as_padded_array())
    @property
    def step_in_yr(self):
        return abs(self.step*self.units)
    @property
    def beginning(self):
        match self.reversed:
            case True:
                return self.stop
            case False:
                return self.start
    @property
    def end(self):
        match self.reversed:
            case True:
                return self.start
            case False:
                return self.stop
    def as_array(self:Self) -> numpy.ndarray:
        match self.reversed:
            case True:
                return numpy.arange(self.start,self.stop,-self.step)
            case False:
                return numpy.arange(self.start,self.stop,self.step)
            case _:
                raise ValueError("Invalid time configuration")
    def as_padded_array(self):
        return numpy.append(self.as_array(),[self.stop,self.stop])
    @staticmethod
    def translate_units(units):
        if units=="yr":
            return 1.0
        elif units=="Myr":
            return 1e6
    
class DynamicTime:
    def __init__(self,start,initial_step,stop,units="yr",max_step=None,min_step=None):
        self.start = start
        self.initial_step = initial_step
        self.stop = stop
        self.units = Time.translate_units(units)
        self.reversed = start>stop
        self.steps = []
        self.max_step = max_step if max_step is not None else initial_step*10
        self.min_step = min_step if min_step is not None else initial_step/10
    def as_array_using_initial(self):
        match self.reversed:
            case True:
                return numpy.arange(self.start,self.stop,-self.initial_step)
            case False:
                return numpy.arange(self.start,self.stop,self.initial_step)
    def as_array_after(self):
        match self.reversed:
            case True:
                return self.start - numpy.cumsum(self.steps)
            case False:
                return self.start + numpy.cumsum(self.steps)
    def inside(self,value):
        match self.reversed:
            case True:
                return value<=self.start and value>self.stop
            case False:
                return value<self.stop and value>=self.start
    def increase_step(self,step):
        return step if step>=self.max_step else step*2
    def decrease_step(self,step):
        return step if step<=self.min_step else step/2
    @property
    def initial_step_in_yr(self):
        return self.initial_step*self.units

class MultiTime:
    def __init__(self,steps,units="yr"):
        self.steps = steps
        self.reversed = bool(steps[0]>steps[-1])
        self.units = Time.translate_units(units)
    @property
    def length(self):
        return len(self.steps)
    def as_array(self: Self) -> numpy.ndarray:
        return numpy.array(self.steps)
    def as_padded_array(self):
        return numpy.append(self.as_array(),[self.steps[-1],self.steps[-1]])
    @property
    def beginning(self):
        match self.reversed:
            case True:
                return self.steps[-1]
            case False:
                return self.steps[0]
    @property
    def end(self):
        match self.reversed:
            case True:
                return self.steps[0]
            case False:
                return self.steps[-1]
    @property
    def initial_step_in_yr(self):
        return numpy.abs(self.steps[0]-self.steps[1])*self.units

class FluxType(Enum):
    input = 1
    oxyhydroxide = 2
    carbonate = 3
    clay = 4
    organic = 5
    sulphide = 6
    
class Perturbation:
    def __init__(self,time,type:FluxType,isotope_flux):
        self.time = time
        self.type = type
        self.isotope_flux = isotope_flux

