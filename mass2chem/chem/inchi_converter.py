import pybel


def convert_inchi_to_formula(inchi_string):
  """
  Converts InChI to formula. Depends on/uses openbabel.
  """
  # We always cast strings because some applications return unicode (such as Django model fields)
  return pybel.readstring('inchi', str(inchi_string)).formula
