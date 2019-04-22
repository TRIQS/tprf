################################################################################
#
# TRIQS: a Toolbox for Research in Interacting Quantum Systems
#
# Copyright (C) 2017 by Hugo U.R. Strand
#
# TRIQS is free software: you can redistribute it and/or modify it under the
# terms of the GNU General Public License as published by the Free Software
# Foundation, either version 3 of the License, or (at your option) any later
# version.
#
# TRIQS is distributed in the hope that it will be useful, but WITHOUT ANY
# WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
# FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
# details.
#
# You should have received a copy of the GNU General Public License along with
# TRIQS. If not, see <http://www.gnu.org/licenses/>.
#
################################################################################

import inspect
import numpy as np

# ----------------------------------------------------------------------
class ParameterCollection(object):

    """ Helper class for handing collections of parameters.

    Parameters
    ----------

    kwargs : dict
        Key-word argument list of parameters.

    Examples
    --------

    A ``ParameterCollection`` has any number of attributes, accessible with the
    dot operator.

    >>> p = ParameterCollection(beta=10., U=1.0, t=1.0)
    >>> print p
    U = 1.0
    beta = 10.0
    t = 1.0
    >>> print p.beta
    10.0
    >>> p.W = 1.2
    >>> print p
    U = 1.0
    W = 1.2
    beta = 10.0
    t = 1.0
    
    and can be stored and loaded to/from TRIQS hdf archives.

    >>> from pytriqs.archive import HDFArchive
    >>> with HDFArchive('data.h5', 'w') as arch: arch['p'] = p
    >>> with HDFArchive('data.h5', 'r') as arch: p_ref = arch['p']
    >>> print p_ref
    U = 1.0
    beta = 10.0
    t = 1.0
    
    """
    
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)

    def items(self):
        return self.__dict__.items()

    def keys(self):
   	return self.__dict__.keys()

    def dict(self):
        return self.__dict__

    def __getitem__(self, key):
        return self.__dict__[key]

    def __reduce_to_dict__(self):
        return self.__dict__

    def _clean_bools(self):
        """ Fix for bug in Triqs that cast bool to numpy.bool_ 
        here we cast all numpy.bools_ to plain python bools """
        
        for key, value in self.items():
            if type(value) == np.bool_:
                self.dict()[key] = bool(value)

    def convert_keys_from_string_to_python(self, dict_key):
        """ pytriqs.archive.HDFArchive incorrectly mangles tuple keys to string
        running this on the affected dict tries to revert this by running eval
        on the string representation. UGLY FIX... """

        d = self.dict()[dict_key]
        d_fix = {}
        for key, value in d.items():
            d_fix[eval(key)] = value            
        self.dict()[dict_key] = d_fix
    
    def grab_attribs(self, obj, keys):
        for key in keys:
            val = getattr(obj, key)
            self.dict()[key] = val

    @classmethod
    def __factory_from_dict__(cls, name, d):
        ret = cls()
        ret.__dict__.update(d)
        ret._clean_bools()
        return ret

    def __str__(self):
        out = ''
        keys = np.sort(self.__dict__.keys()) # sort keys
        for key in keys:
            value = self.__dict__[key]
            if type(value) is ParameterCollection:
                pc_list = str(value).splitlines()
                pc_txt = ''.join([ key + '.' + row + '\n' for row in pc_list ])
                out += pc_txt
            else:
                str_value = str(value)

                # Cut things that take more than five rows
                str_value_lines = str_value.splitlines()
                max_lines = 10
                if len(str_value_lines) > max_lines:
                    str_value = '\n'.join(str_value_lines[:max_lines] + ['...'])
                
                out += ''.join([key, ' = ', str_value]) + '\n'
        return out

    def get_my_name(self):
        ans = []
        frame = inspect.currentframe().f_back
        tmp = dict(frame.f_globals.items() + frame.f_locals.items())
        for k, var in tmp.items():
            if isinstance(var, self.__class__):
                if hash(self) == hash(var):
                    ans.append(k)
        return ans


# -- Register ParameterCollection in Triqs hdf_archive_schemes

from pytriqs.archive.hdf_archive_schemes import register_class 
register_class(ParameterCollection)

# ----------------------------------------------------------------------
class ParameterCollections(object):

    r""" Helper class for handing a series of collections of parameters.

    Parameters
    ----------

    objects : list
        List of ``ParameterCollection`` instances.

    Examples
    --------

    The ``ParameterCollections`` class makes it easy to get vectors of 
    parameters from a list of ``ParameterCollection`` objects.

    >>> p1 = ParameterCollection(beta=10., U=1.0, t=1.0)
    >>> p2 = ParameterCollection(beta=5., U=2.0, t=1.337)
    >>> ps = ParameterCollections(objects=[p1, p2])
    >>> print ps.beta
    [10.  5.]
    >>> print ps.U
    [1. 2.]

    """
    
    def __init__(self, objects=[]):
        self.objects = objects

    def append(self, obj):
        self.objects.append(obj)
        
    def sort_on(self, attr):
        val = self.getattr_from_objects(attr)
        sidx = np.argsort(val)
        self.set_sorted_order(sidx)

    def set_sorted_order(self, sorted_idx):
        sidx = np.array(sorted_idx)
        self.objects = list(np.array(self.objects)[sidx])

    def getattr_from_objects(self, attr):
        return np.array([getattr(o, attr) for o in self.objects ])
    
    def __getattr__(self, attr):
        return self.getattr_from_objects(attr)

    def __reduce_to_dict__(self):
        return self.__dict__

    @classmethod
    def __factory_from_dict__(cls, name, d):
        ret = cls(d['objects'])
        return ret    

# ----------------------------------------------------------------------
# -- Register ParameterCollection in Triqs hdf_archive_schemes

from pytriqs.archive.hdf_archive_schemes import register_class 
register_class(ParameterCollections)
