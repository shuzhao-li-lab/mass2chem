from __future__ import print_function

from .inchi_converter import convert_inchi_to_formula
from .adducts import positive_mode_adducts, negative_mode_adducts
from . import molmass


class FormulaBase(object):
  H_ADDUCT = 1.0073

  def __init__(self, formula=None):
    self._formula = formula

  @property
  def formula(self):
    return molmass.Formula(self._formula)

  def __str__(self):
    return self.formula

  @property
  def monoisotopic_mass(self):
    """
    >>> a = FormulaBase('C4H4Na2O4')
    >>> print('%.9f' % a.monoisotopic_mass)
    161.990497957
    """
    return self.formula.isotope.mass

  def positive_adduct(self, what):
    """
    >>> a = FormulaBase('C4H4Na2O4')
    >>> print('%.9f' % a.positive_adduct('M+H'))
    162.997797957
    """

    class StupidPythonClosure(object):
      mz = None

    obj = StupidPythonClosure()

    def func(name, mul, off):
      if name == what:
        obj.mz = self.monoisotopic_mass * mul + off

    positive_mode_adducts(func)
    return obj.mz

  def negative_adduct(self, what):
    """
    >>> a = FormulaBase('C4H4Na2O4')
    >>> print('%.9f' % a.negative_adduct('M-H'))
    160.983197957
    """

    class StupidPythonClosure(object):
      mz = None
    obj = StupidPythonClosure()

    def func(name, mul, off):
      if name == what:
        obj.mz = self.monoisotopic_mass * mul + off
    negative_mode_adducts(func)
    return obj.mz


# For backwards compatability
Formula_Base = FormulaBase


class InChI(FormulaBase):
  """
  >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
  >>> print('%.9f' % a.monoisotopic_mass)
  161.990497957
  """

  def __init__(self, inchi, formula=None):
    if not inchi.startswith("InChI="):
      raise ValueError("Invalid InChI %s, must start with InChI=" % (inchi,))

    self.__inchi = inchi
    inchi = inchi[6:]
    self.__layers = inchi.split('/')[1:]

    super(InChI, self).__init__(formula)

  @property
  def inchi(self):
    """
    >>> a = InChI('InChI=1S/H2O/h1H2')
    >>> print(a.inchi)
    InChI=1S/H2O/h1H2
    """
    return self.__inchi

  @property
  def formula(self):
    """
    >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
    >>> print(a.formula)
    C4H4Na2O4
    >>> a = InChI(
    ...   'InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2',
    ...   'H2O',
    ... )
    >>> print(a.formula)
    H2O
    """
    if self._formula is None:
      self._formula = convert_inchi_to_formula(self.__inchi)
    return molmass.Formula(self._formula)

  @property
  def spectrum(self):
    """
    >>> a = InChI('InChI=1S/C4H6O4.2Na/c5-3(6)1-2-4(7)8;;/h1-2H2,(H,5,6)(H,7,8);;/q;2*+1/p-2')
    >>> spc = a.spectrum
    >>> print(spc)
    Relative mass    Fraction %      Intensity
    161.99050         94.816904     100.000000
    162.99385          4.102055       4.326291
    162.99471          0.144473       0.152370
    162.99677          0.043621       0.046005
    163.99474          0.779393       0.821997
    163.99721          0.066550       0.070188
    163.99807          0.006250       0.006592
    163.99893          0.000041       0.000044
    164.00013          0.001887       0.001990
    164.00099          0.000017       0.000018
    164.99810          0.033719       0.035562
    164.99896          0.000891       0.000939
    165.00063          0.000569       0.000601
    165.00117          0.000140       0.000148
    165.00229          0.000002       0.000002
    165.00435          0.000001       0.000001
    165.99899          0.002402       0.002534
    166.00145          0.000547       0.000577
    166.00232          0.000019       0.000020
    166.00439          0.000004       0.000004
    167.00234          0.000052       0.000055
    167.00321          0.000001       0.000001
    167.00481          0.000001       0.000001
    168.00324          0.000002       0.000003
    168.00570          0.000001       0.000001
    >>> for v in sorted(spc.values()): print('%.4f %.6f' % (v[0], v[1]))
    161.9905 0.948169
    162.9939 0.041021
    162.9947 0.001445
    162.9968 0.000436
    163.9947 0.007794
    163.9972 0.000666
    163.9981 0.000063
    163.9989 0.000000
    164.0001 0.000019
    164.0010 0.000000
    164.9981 0.000337
    164.9990 0.000009
    165.0006 0.000006
    165.0012 0.000001
    165.0023 0.000000
    165.0043 0.000000
    165.9990 0.000024
    166.0015 0.000005
    166.0023 0.000000
    166.0032 0.000000
    166.0044 0.000000
    166.0053 0.000000
    167.0023 0.000001
    167.0032 0.000000
    167.0048 0.000000
    167.0055 0.000000
    168.0032 0.000000
    168.0057 0.000000
    """
    return self.formula.high_resolution_spectrum(3, minfract=1e-6)


if __name__ == '__main__':
  import doctest
  doctest.testmod()
