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

class ParameterCollection(object):

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
