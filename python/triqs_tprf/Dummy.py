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

class Dummy(object):
    """ Dummy class that is writable to a HDFArchive. """
    
    def __init__(self):
        pass

    def __reduce_to_dict__(self):
        return self.__dict__

    @classmethod
    def __factory_from_dict__(cls, name, d):
        ret = cls()
        ret.__dict__.update(d)
        return ret
        
from h5.formats import register_class 
register_class(Dummy)
