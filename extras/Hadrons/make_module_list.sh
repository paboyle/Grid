#!/usr/bin/env bash

echo 'modules_cc =\' > modules.inc
find Modules -name '*.cc' -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
echo 'modules_hpp =\' >> modules.inc
find Modules -name '*.hpp' -type f -print | sed 's/^/  /;$q;s/$/ \\/' >> modules.inc
echo '' >> modules.inc
rm -f Modules.hpp
echo "/*************************************************************************************

Grid physics library, www.github.com/paboyle/Grid 

Source file: extras/Hadrons/Modules.hpp

Copyright (C) 2015
Copyright (C) 2016
Copyright (C) 2017

Author: Antonin Portelli <antonin.portelli@me.com>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation; either version 2 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License along
with this program; if not, write to the Free Software Foundation, Inc.,
51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

See the full license in the file \"LICENSE\" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
" > Modules.hpp
for f in `find Modules -name '*.hpp'`; do
	echo "#include <Grid/Hadrons/${f}>" >> Modules.hpp
done
