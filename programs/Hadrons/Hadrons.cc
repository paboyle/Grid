/*
 * Hadrons.cc, part of Grid
 *
 * Copyright (C) 2015 Antonin Portelli
 *
 * LatCore is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * LatCore is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with LatCore.  If not, see <http://www.gnu.org/licenses/>.
 */

#include <Hadrons/Application.hpp>

using namespace std;
using namespace Hadrons;

int main(int argc, char *argv[])
{
    Application application(argc, argv);
    
    application.run();
    
    return EXIT_SUCCESS;
}
