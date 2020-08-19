#!/usr/bin/python3

import re
import argparse
import sys

# Grid for A64FX
#
# * should align std::vector to (multiples of) cache block size = 256 bytes

# place benchmark runtime in cycles here !
measured_cycles = 690 #1500 #775 #1500


# command line parser
parser = argparse.ArgumentParser(description="Dslash generator.")
parser.add_argument("--single", action="store_true", default="False")
parser.add_argument("--double", action="store_true", default="True")
parser.add_argument("--debug", action="store_true", default="False")
parser.add_argument("--gridbench", action="store_true", default="False")
args = parser.parse_args()

print(args)

ASM_LOAD_CHIMU = True       # load chimu
ASM_LOAD_GAUGE = True       # load gauge
ASM_LOAD_TABLE = True       # load table
ASM_STORE = True            # store result

# Disable all loads and stores in asm for benchmarking purposes
#DISABLE_ASM_LOAD_STORE = True
DISABLE_ASM_LOAD_STORE = False

if DISABLE_ASM_LOAD_STORE:
    ASM_LOAD_CHIMU = True       # load chimu
    ASM_LOAD_GAUGE = True       # load gauge
    ASM_LOAD_TABLE = True       # load table
    ASM_STORE = False            # store result

# Alternative implementation using PROJ specific loads works,
# but be careful with predication

ALTERNATIVE_LOADS = False
#ALTERNATIVE_LOADS = not ALTERNATIVE_LOADS      # True

# Alternative register mapping,
# must use with my_wilson4.h and my_wilson4pf.h

ALTERNATIVE_REGISTER_MAPPING = False
#ALTERNATIVE_REGISTER_MAPPING = not ALTERNATIVE_REGISTER_MAPPING

if ALTERNATIVE_REGISTER_MAPPING == True:
    ALTERNATIVE_LOADS = False

# use movprfx
MOVPRFX = False
MOVPRFX = not MOVPRFX


PREFETCH = False
PREFETCH = not PREFETCH # True

PRECISION = 'double'   # DP by default
PRECSUFFIX = 'A64FXd'
if args.single == True:
    PRECISION = 'single'
    PRECSUFFIX = 'A64FXf'

_DEBUG = False #True       # insert debugging output
if args.debug == True:
    _DEBUG = True

GRIDBENCH = False
if args.gridbench == True:
    GRIDBENCH = True

print("PRECISION                    = ", PRECISION)
print("DEBUG                        = ", _DEBUG)
print("ALTERNATIVE_LOADS            = ", ALTERNATIVE_LOADS)
print("ALTERNATIVE_REGISTER_MAPPING = ", ALTERNATIVE_REGISTER_MAPPING)
print("MOVPRFX                      = ", MOVPRFX)
print("DISABLE_ASM_LOAD_STORE       = ", DISABLE_ASM_LOAD_STORE)
print("GRIDBENCH                    = ", GRIDBENCH)

print("")

#sys.exit(0)


#_DEBUG = True       # insert debugging output

FETCH_BASE_PTR_COLOR_OFFSET = 2  # offset for scalar plus signed immediate addressing
STORE_BASE_PTR_COLOR_OFFSET = 2

# 64-bit gp register usage !!! armclang 20.0 complains about the register choice !!!
# table address: x30
# data address:  x29
# store address: x28
# debug address: r8

# Max performance of complex FMA using FCMLA instruction
# is 25% peak.
#
# Issue latency of FCMLA is 2 cycles.
# Need 2 FCMLA instructions for complex FMA.
# Complete complex FMA takes 4 cycles.
# Peak throughput is 4 * 8 Flops DP = 32 Flops DP in 4 cycles.
# A64FX FMA throughput is 4 * 8 * 2 * 2 = 132 Flops DP in 4 cycles.
# -> 25% peak FMA
#
# In:  3x 512 bits = 192 bytes
# Out: 1x 512 bits = 64 bytes
# Tot: 4x 512 bits = 256 bytes
#
# 256 bytes * 2.2 GHz = 563.2 GB/s (base 10), 524 GB/s (base 2)

OPT = """
* interleave prefetching and compute in MULT_2SPIN
* could test storing U's in MULT_2SPIN to L1d for cache line update
* structure reordering: MAYBEPERM after MULT_2SPIN ?
"""

filename = 'XXX'
LEGAL = """/*************************************************************************************

    Grid physics library, www.github.com/paboyle/Grid

    Source file: {}

    Copyright (C) 2020

Author: Nils Meyer <nils.meyer@ur.de>

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

    See the full license in the file "LICENSE" in the top level distribution directory
*************************************************************************************/
/*  END LEGAL */
"""

class Register:

    def __init__(self, variable, asmreg='X', predication=False):
        global d
        x = 'Y'
        if predication == False:
            x = asmreg # + d['asmsuffix']
        else:
            x = asmreg
        self.asmreg = x
        self.asmregwithsuffix = asmreg + d['asmsuffix']
        self.asmregbyte = asmreg + '.b'
        self.name = variable
        self.asmname = variable
        self.asmnamebyte = variable + '.b'
        self.predication = predication

        d['registers'] += 1

    def define(self, statement):
        global d
        d['C'] += F'#define {self.name} {statement}'
        #d['A'] += F'#define {self.name} {statement}'

    def declare(self, predication=False):
        global d

        if self.predication == False:
            d['C'] += F'    Simd {self.name};               \\\n'

            predtype = 'svfloat64_t'
            if PRECISION == 'single':
                predtype = 'svfloat32_t'

            d['I'] += F'    {predtype} {self.name};        \\\n'
        else:
            d['I'] += F'    svbool_t {self.name};        \\\n'
        #d['A'] += F'#define {self.name} {self.asmreg} \n'

    def loadpredication(self, target='A'):
        global d
        if (target == 'A'):
            d['A'] += F'    "ptrue {self.asmregwithsuffix} \\n\\t" \\\n'
            d['asmclobber'].append(F'"{self.asmreg}"')

    def loadtable(self, t):
        global d
        d['load'] += d['factor']
        gpr = d['asmtableptr']

        cast = 'uint64_t'
        #asm_opcode = 'ld1d'
        #if PRECISION == 'single':
        #   asm_opcode = 'ld1w'
        #    cast = 'uint32_t'
        asm_opcode = 'ldr'
        if PRECISION == 'single':
            asm_opcode = 'ldr'
            cast = 'uint32_t'

        d['I'] += F'    {self.name} = svld1(pg1, ({cast}*)&lut[{t}]);  \\\n'

        # using immediate index break-out works
        if asm_opcode == 'ldr':
            # ldr version
            d['A'] += F'    "{asm_opcode} {self.asmreg}, [%[tableptr], %[index], mul vl] \\n\\t" \\\n'
        else:
            # ld1 version
            d['A'] += F'    "{asm_opcode} {{ {self.asmregwithsuffix} }}, {pg1.asmreg}/z, [%[tableptr], %[index], mul vl] \\n\\t" \\\n'

        d['asminput'].append(F'[tableptr] "r" (&lut[0])')
        d['asminput'].append(F'[index] "i" ({t})')
        d['asmclobber'].append(F'"memory"')
        d['asmclobber'].append(F'"cc"')

    def load(self, address, target='ALL', cast='float64_t', colors=3, offset=FETCH_BASE_PTR_COLOR_OFFSET):
        global d
        d['load'] += d['factor']
        indices = re.findall(r'\d+', address)
        index = (int(indices[0]) - offset) * colors + int(indices[1])

        #asm_opcode = 'ld1d'
        #if PRECISION == 'single':
        #asm_opcode = 'ld1w'
        #    cast = 'float32_t'

        asm_opcode = 'ldr'
        if PRECISION == 'single':
            asm_opcode = 'ldr'
            cast = 'float32_t'

        gpr = d['asmfetchbaseptr']
        intrinfetchbase = d['intrinfetchbase']
        if (target in ['ALL', 'C']):
            d['C'] += F'    {self.name} = {address};        \\\n'
        if (target in ['ALL', 'I']):
#            d['I'] += F'    {self.name} = svldnt1(pg1, ({cast}*)({intrinfetchbase} + {index} * 64));  \\\n'
            d['I'] += F'    {self.name} = svld1(pg1, ({cast}*)({intrinfetchbase} + {index} * 64));  \\\n'
        if (target in ['ALL', 'A']):
            if asm_opcode == 'ldr':
                d['A'] += F'    "{asm_opcode} {self.asmreg}, [%[fetchptr], {index}, mul vl] \\n\\t" \\\n'
            else:
                d['A'] += F'    "{asm_opcode} {{ {self.asmregwithsuffix} }}, {pg1.asmreg}/z, [%[fetchptr], {index}, mul vl] \\n\\t" \\\n'

    def store(self, address, cast='float64_t', colors=3, offset=STORE_BASE_PTR_COLOR_OFFSET):
        global d
        d['store'] += d['factor']
        indices = re.findall(r'\d+', address)
        index = (int(indices[0]) - offset) * colors + int(indices[1])

        #asm_opcode = 'stnt1d'
        #if PRECISION == 'single':
        #    asm_opcode = 'stnt1w'
        #    cast = 'float32_t'
        asm_opcode = 'str'
        if PRECISION == 'single':
            asm_opcode = 'str'
            cast = 'float32_t'

        intrinstorebase = d['intrinstorebase']

        d['C'] += F'    {address} = {self.name};        \\\n'
        #d['I'] += F'    svstnt1(pg1, ({cast}*)({intrinstorebase} + {index} * 64), {self.name});  \\\n'
        d['I'] += F'    svst1(pg1, ({cast}*)({intrinstorebase} + {index} * 64), {self.name});  \\\n'
        if asm_opcode == 'str':
            d['A'] += F'    "{asm_opcode} {self.asmreg}, [%[storeptr], {index}, mul vl] \\n\\t" \\\n'
        else:
            d['A'] += F'    "{asm_opcode} {{ {self.asmregwithsuffix} }}, {pg1.asmreg}, [%[storeptr], {index}, mul vl] \\n\\t" \\\n'

    def movestr(self, str):
        global d
        #d['move'] += d['factor']
        d['I'] += F'    {self.name} = {str};        \\\n'

    def move(self, op1):
        global d
        d['move'] += d['factor']
        d['C'] += F'    {self.name} = {op1.name};   \\\n'
        d['I'] += F'    {self.name} = {op1.name};        \\\n'
        d['A'] += F'    "mov {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix} \\n\\t" \\\n'

    # a = a + b , a = b + c
    def add(self, op1, op2=None):
        global d
        d['add'] += d['factor']
        if op2 is None:
            d['C'] += F'    {self.name} = {self.name} + {op1.name};   \\\n'
            d['I'] += F'    {self.name} = svadd_x(pg1, {self.name}, {op1.name}); \\\n'
            d['A'] += F'    "fadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op1.asmregwithsuffix} \\n\\t"  \\\n'
        else:
            d['C'] += F'    {self.name} = {op1.name} + {op2.name};    \\\n'
            d['I'] += F'    {self.name} = svadd_x(pg1, {op1.name}, {op2.name});  \\\n'
            d['A'] += F'    "fadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix} \\n\\t"  \\\n'

    # a = a -b , a = b - c
    def sub(self, op1, op2=None):
        global d
        d['sub'] += d['factor']
        if op2 is None:
            d['C'] += F'    {self.name} = {self.name} - {op1.name};    \\\n'
            d['I'] += F'    {self.name} = svsub_x(pg1, {self.name}, {op1.name}); \\\n'
            d['A'] += F'    "fsub {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op1.asmregwithsuffix} \\n\\t" \\\n'
        else:
            d['C'] += F'    {self.name} = {op1.name} - {op2.name};     \\\n'
            d['I'] += F'    {self.name} = svsub_x(pg1, {op1.name}, {op2.name});  \\\n'
            d['A'] += F'    "fsub {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix} \\n\\t" \\\n'

    # a = a * b , a = b * c
    def mul(self, op1, op2):
        global d
        d['mul'] += 2 * d['factor']
        d['C'] += F'    {self.name} = {op1.name} * {op2.name};                \\\n'
        d['I'] += F'    {self.name} = __svzero({self.name}); \\\n'
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 0); \\\n'
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 90); \\\n'
        d['A'] += F'    "mov {self.asmregwithsuffix} , 0 \\n\\t" \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 0 \\n\\t" \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'

    def mul0(self, op1, op2, op3=None, constructive=False):
        global d
        d['mul'] += d['factor']

        # no movprfx intrinsics support
        if constructive == True:
            d['movprfx'] += d['factor']
            d['I'] += F'    {self.name} = svcmla_x(pg1, {op1.name}, {op2.name}, {op3.name}, 0); \\\n'
            d['A'] += F'    "movprfx {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix} \\n\\t" \\\n'
            d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op2.asmregwithsuffix}, {op3.asmregwithsuffix}, 0 \\n\\t" \\\n'
        else:
            d['C'] += F'    {self.name} = {op1.name} * {op2.name};                \\\n'
            d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 0); \\\n'
            d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 0 \\n\\t" \\\n'

    def mul1(self, op1, op2):
        global d
        d['mul'] += d['factor']
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 90); \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'

    def mac(self, op1, op2):
        global d
        d['mac'] += 2 * d['factor']
        d['C'] += F'    {self.name} = {self.name} + {op1.name} * {op2.name};    \\\n'
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 0); \\\n'
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 90); \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 0 \\n\\t" \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'

    def mac0(self, op1, op2):
        global d
        d['mac'] += d['factor']
        d['C'] += F'    {self.name} = {self.name} + {op1.name} * {op2.name};    \\\n'
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 0); \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 0 \\n\\t" \\\n'

    def mac1(self, op1, op2):
        global d
        d['mac'] += d['factor']
        d['I'] += F'    {self.name} = svcmla_x(pg1, {self.name}, {op1.name}, {op2.name}, 90); \\\n'
        d['A'] += F'    "fcmla {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'

    def zero(self, zeroreg=False):
        d['zero'] += d['factor']
        d['C'] += F'    {self.name} = 0; \\\n'
        #d['I'] += F'    {self.name} = __svzero({self.name}); \\\n'   only armclang

        if PRECISION == 'double':
            d['I'] += F'    {self.name} = svdup_f64(0.); \\\n'
        else:
            d['I'] += F'    {self.name} = svdup_f32(0.); \\\n'

        if zeroreg == True:
            d['A'] += F'    "fmov {self.asmregwithsuffix} , 0 \\n\\t" \\\n'
        else:
            #using mov z, zero0      issue 1c, FLA, latency 6c
            #d['A'] += F'    "mov {self.asmregwithsuffix} , {zero0.asmregwithsuffix} \\n\\t" \\\n'

            #using mov z, 0      issue 1c, FLA, latency 6c
            d['A'] += F'    "fmov {self.asmregwithsuffix} , 0 \\n\\t" \\\n'

            #using xor z, z, z   issue 0.5c, FL*, latency 4c
            #d['A'] += F'    "eor  {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {self.asmregwithsuffix} \\n\\t" \\\n'

            #using and z, z, zero0  issue 0.5c, FL*, latency 4c
            #d['A'] += F'    "and {self.asmregwithsuffix}, {self.asmregwithsuffix} , {zero0.asmregwithsuffix} \\n\\t" \\\n'

            #using sub z, z, z    issue 0.5c, FL*, latency 9c
            #d['A'] += F'    "sub  {self.asmregwithsuffix}, {self.asmregwithsuffix}, {self.asmregwithsuffix} \\n\\t" \\\n'

    # without table
    def timesI(self, op1, tempreg=None, tablereg=None):
        global d
        d['timesI'] += d['factor']
        d['C'] += F'    {self.name} = timesI({op1.name});    \\\n'
        # correct if DEBUG enabled, wrong if DEBUG disabled; no idea what's causing this
        #table.load('table2', target='I', cast='uint64_t')
        #d['I'] += F'    {self.name} = svtbl({op1.name}, {tablereg.name});  \\\n'
        #d['I'] += F'    {self.name} = svneg_x(pg2, {self.name});  \\\n'
        # timesI using trn tested, works but tbl should be faster
        d['I'] += F'    {tempreg.name} = svtrn2({op1.name}, {op1.name});   \\\n'
        d['I'] += F'    {tempreg.name} = svneg_x(pg1, {tempreg.name});   \\\n'
        d['I'] += F'    {self.name} = svtrn1({tempreg.name}, {op1.name});   \\\n'
        d['A'] += F'    "trn2 {tempreg.asmregwithsuffix}, {op1.asmregwithsuffix}, {op1.asmregwithsuffix} \\n\\t"   \\\n'
        d['A'] += F'    "fneg {tempreg.asmregwithsuffix}, {pg1.asmreg}/m, {tempreg.asmregwithsuffix} \\n\\t"   \\\n'
        d['A'] += F'    "trn1 {self.asmregwithsuffix}, {tempreg.asmregwithsuffix}, {op1.asmregwithsuffix} \\n\\t"   \\\n'

    def addTimesI(self, op1, op2=None, constructive=False):
        global d
        d['addTimesI'] += d['factor']

        if op2 is None:
            d['C'] += F'    {self.name} = {self.name} + timesI({op1.name});    \\\n'
        else:
            d['C'] += F'    {self.name} = {op1.name} + timesI({op2.name});    \\\n'

        # no movprfx intrinsics support
        if constructive == True:
            d['movprfx'] += d['factor']
            d['I'] += F'    {self.name} = svcadd_x(pg1, {op1.name}, {op2.name}, 90);   \\\n'
            d['A'] += F'    "movprfx {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix} \\n\\t" \\\n'
            d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'
        else:
            if op2 is None:
                d['C'] += F'    {self.name} = {self.name} + timesI({op1.name});    \\\n'
                d['I'] += F'    {self.name} = svcadd_x(pg1, {self.name}, {op1.name}, 90);   \\\n'
                d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op1.asmregwithsuffix}, 90 \\n\\t" \\\n'
            else:
                d['C'] += F'    {self.name} = {op1.name} + timesI({op2.name});    \\\n'
                d['I'] += F'    {self.name} = svcadd_x(pg1, {op1.name}, {op2.name}, 90);   \\\n'
                d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 90 \\n\\t" \\\n'

    def subTimesI(self, op1, op2=None, constructive=False):
        global d
        d['subTimesI'] += d['factor']

        # no movprfx intrinsics support
        if constructive == True:
            d['movprfx'] += d['factor']
            d['I'] += F'    {self.name} = svcadd_x(pg1, {op1.name}, {op2.name}, 270);   \\\n'
            d['A'] += F'    "movprfx {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix} \\n\\t" \\\n'
            d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op2.asmregwithsuffix}, 270 \\n\\t" \\\n'
        else:
            if op2 is None:
                d['C'] += F'    {self.name} = {self.name} - timesI({op1.name});    \\\n'
                d['I'] += F'    {self.name} = svcadd_x(pg1, {self.name}, {op1.name}, 270);   \\\n'
                d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {self.asmregwithsuffix}, {op1.asmregwithsuffix}, 270 \\n\\t" \\\n'
            else:
                d['C'] += F'    {self.name} = {op1.name} - timesI({op2.name});    \\\n'
                d['I'] += F'    {self.name} = svcadd_x(pg1, {op1.name}, {op2.name}, 270);   \\\n'
                d['A'] += F'    "fcadd {self.asmregwithsuffix}, {pg1.asmreg}/m, {op1.asmregwithsuffix}, {op2.asmregwithsuffix}, 270 \\n\\t" \\\n'

    # timesMinusI is not used, def is probably wrong !!!!  OPTIMIZATION with table
    def timesMinusI(self, op1):
        global d
        d['timesMinusI'] += d['factor']
        d['C'] += F'    {self.name} = timesMinusI({self.name});    \\\n'
        d['I'] += F'    {self.name} = svtrn1({op1.name}, {op1.name});   \\\n'
        d['I'] += F'    {self.name} = svneg_x(pg1, {self.name});   \\\n'
        d['I'] += F'    {self.name} = svtrn1({op1.name}, {self.name});   \\\n'

    def permute(self, dir, tablereg=None):
        global d
        d['permutes'] += d['factor']

        d['C'] += F'    permute{dir}({self.name}, {self.name});    \\\n'

        d['I'] += F'    {self.name} = svtbl({self.name}, {tablereg.name});    \\\n'
        d['A'] += F'    "tbl {self.asmregwithsuffix}, {{ {self.asmregwithsuffix} }}, {tablereg.asmregwithsuffix} \\n\\t"  \\\n'

        # if dir == 0:
        #     d['I'] += F'    {self.name} = svext({self.name}, {self.name}, 4); \\\n'
        #     # this might not work, see intrinsics assembly
        #     # d['A'] += F'    ext {self.name}, {self.name}, {self.name}, #4 \\\n'
        #     # use registers directly
        #     d['A'] += F'    "ext {self.asmregbyte}, {self.asmregbyte}, {self.asmregbyte}, 32 \\n\\t" \\\n'
        #
        # elif dir in [1, 2]:
        #     d['I'] += F'    {self.name} = svtbl({self.name}, {tablereg.name});    \\\n'
        #     d['A'] += F'    "tbl {self.asmregwithsuffix}, {{ {self.asmregwithsuffix} }}, {tablereg.asmregwithsuffix} \\n\\t"  \\\n'

    def debug(self):
        global d
        typecast = d['cfloat']
        gpr = d['asmdebugptr']
        vregs = d['asmclobberlist']
        if (d['debug'] == True):
            d['C'] += F'std::cout << "{self.name} -- " <<  {self.name} << std::endl; \\\n'

            d['I'] += F'svst1(pg1, ({typecast}*)&debugreg.v, {self.name}); \\\n'
            d['I'] += F'std::cout << "{self.name} -- " <<  debugreg << std::endl; \\\n'
            #d['I'] += F'std::cout << "{self.name} -- " <<  {self.name} << std::endl; \\\n'

            d['A'] += F'asm ( \\\n'
            d['A'] += F'    " DMB SY \\n\\t " " DSB SY \\n\\t " " ISB SY \\n\\t " \\\n'   # memory barrier
            d['A'] += F'    "str {self.asmreg}, [%[ptr]] \\n\\t" \\\n'
            d['A'] += F'    " DMB SY \\n\\t " " DSB SY \\n\\t " " ISB SY \\n\\t " \\\n'   # memory barrier
            d['A'] += F'    : "=m" (debugreg.v)  \\\n'
            d['A'] += F'    : [ptr] "r" (&debugreg.v)  \\\n'
            d['A'] += F'    : "p5", "cc", "memory" \\\n'
            d['A'] += F'); \\\n'
            d['A'] += F'std::cout << "{self.name} -- " <<  debugreg << std::endl; \\\n'
            # this form of addressing is not valid!
            #d['A'] += F'    "str {self.asmreg}, %[ptr] \\n\\t" \\\n'
# end Register

def define(s, target='ALL'):
    x = F'#define {s}  \n'
    global d
    if (target in ['ALL', 'C']):
        d['C'] += x
    if (target in ['ALL', 'I']):
        d['I'] += x
    if (target in ['ALL', 'A']):
        d['A'] += x

def definemultiline(s):
    x = F'#define {s}  \\\n'
    global d
    d['C'] += x
    d['I'] += x
    d['A'] += x

def write(s, target='ALL'):
    x = F'{s}\n'
    global d
    if (target in ['ALL', 'C']):
        d['C'] += x
    if (target in ['ALL', 'I']):
        d['I'] += x
    if (target in ['ALL', 'A']):
        d['A'] += x

def curlyopen():
    write(F'{{ \\')

def curlyclose():
    write(F'}}')

def newline(target='ALL'):
    global d

    if target == 'A':
        if d['A'][-2:] == '\\\n':
            d['A'] = d['A'][:-2] + '\n\n'
    else:
        if d['C'][-2:] == '\\\n':
            d['C'] = d['C'][:-2] + '\n\n'
        if d['I'][-2:] == '\\\n':
            d['I'] = d['I'][:-2] + '\n\n'
        if d['A'][-2:] == '\\\n':
            d['A'] = d['A'][:-2] + '\n\n'

# load the base pointer for fetches
def fetch_base_ptr(address, target='A'):
    global d
    #d['load'] += d['factor']

    # DEBUG
    #colors=3
    #indices = re.findall(r'\d+', address)
    #index = (int(indices[0]) - FETCH_BASE_PTR_COLOR_OFFSET) * colors + int(indices[1])
    #print(F'{address}                      (base)')

    vregs = d['asmclobberlist']
    if target == 'A':
        d['asminput'].append(F'[fetchptr] "r" ({address})')
        d['asmclobber'].extend(vregs)
        d['asmclobber'].append(F'"memory"')
        d['asmclobber'].append(F'"cc"')
    if target == 'I':
        #print("intrinfetchbase = ", address)
        d['intrinfetchbase'] = address

# load the base pointer for stores
def store_base_ptr(address, target='A'):
    global d
    #d['load'] += d['factor']
    gpr = d['asmstorebaseptr']
    vregs = d['asmclobberlist']
    if target == 'A':
        d['asminput'].append(F'[storeptr] "r" ({address})')
        d['asmclobber'].extend(vregs)
        d['asmclobber'].append(F'"memory"')
        d['asmclobber'].append(F'"cc"')
    if target == 'I':
        d['intrinstorebase'] = address

def prefetch_L1(address, offset):
    global d
    multiplier = 4  # offset in CL, have to multiply by 4
    policy = "PLDL1STRM"     # weak
    #policy = "PLDL1KEEP"     # strong

    d['I'] += F'    svprfd(pg1, (int64_t*)({address} + {offset * multiplier * 64}), SV_{policy}); \\\n'
    d['A'] += F'    "prfd {policy}, {pg1.asmreg}, [%[fetchptr], {offset * multiplier}, mul vl] \\n\\t" \\\n'

def prefetch_L2(address, offset):
    global d
    multiplier = 4  # offset in CL, have to multiply by 4
    policy = "PLDL2STRM"     # weak
    #policy = "PLDL2KEEP"     # strong

    d['I'] += F'    svprfd(pg1, (int64_t*)({address} + {offset * multiplier * 64}), SV_{policy}); \\\n'
    d['A'] += F'    "prfd {policy}, {pg1.asmreg}, [%[fetchptr], {offset * multiplier}, mul vl] \\n\\t" \\\n'
    #d['A'] +=

def prefetch_L2_store(address, offset):
    global d
    multiplier = 4  # offset in CL, have to multiply by 4
    policy = "PSTL2STRM"     # weak
    #policy = "PSTL2KEEP"     # strong

    d['I'] += F'    svprfd(pg1, (int64_t*)({address} + {offset * multiplier * 64}), SV_{policy}); \\\n'
    d['A'] += F'    "prfd {policy}, {pg1.asmreg}, [%[fetchptr], {offset * multiplier}, mul vl] \\n\\t" \\\n'

def prefetch_L1_store(address, offset):
    global d
    multiplier = 4  # offset in CL, have to multiply by 4
    policy = "PSTL1STRM"     # weak
    #policy = "PSTL2KEEP"     # strong

    d['I'] += F'    svprfd(pg1, (int64_t*)({address} + {offset * multiplier * 64}), SV_{policy}); \\\n'
    d['A'] += F'    "prfd {policy}, {pg1.asmreg}, [%[fetchptr], {offset * multiplier}, mul vl] \\n\\t" \\\n'


def asmopen():
    #write('asm volatile ( \\', target='A')
    write('asm ( \\', target='A')

    # DEBUG
    #write(F'    " DMB SY \\n\\t " " DSB SY \\n\\t " " ISB SY \\n\\t " \\', target='A')   # memory barrier
    #write('asm volatile ( \\', target='A')

def asmclose():
    global d

    #print(d['asminput'])

    asmin  = d['asminput']
    asmin_s = ''
    if len(asmin) > 0:
        asmin = list(dict.fromkeys(asmin))   # remove duplicates
        #print(asmin)
        for el in asmin:
            asmin_s += el + ','
        asmin_s = asmin_s[:-1]
        #print("-> ", asmin_s)

    d['asminput'] = []

    asmout = d['asmoutput']
    asmout_s = ''
    if len(asmout) > 0:
        asmout = list(dict.fromkeys(asmout))   # remove duplicates
        for el in asmout:
            asmout_s += el + ','
        asmout_s = asmout_s[:-1]

    d['asmoutput'] = []

    # DEBUG  put all regs into clobber by default
    d['asmclobber'].extend(d['asmclobberlist'])

    asmclobber = d['asmclobber']
    asmclobber_s = ''
    #print(asmclobber)
    if len(asmclobber) > 0:
        asmclobber = list(dict.fromkeys(asmclobber))   # remove duplicates
        for el in asmclobber:
            asmclobber_s += el + ','
        asmclobber_s = asmclobber_s[:-1]

    d['asmclobber'] = []

    # DEBUG
    #write(F'    " DMB SY \\n\\t " " DSB SY \\n\\t " " ISB SY \\n\\t " \\', target='A')   # memory barrier


    write(F'    : {asmout_s} \\', target='A')
    write(F'    : {asmin_s} \\', target='A')
    write(F'    : {asmclobber_s} \\', target='A')
    write('); \\', target='A')

# --------------------------------------------------------------------------------

# string of vector registers to be used in clobber list
#clobberlist = ['"p0"']
clobberlist = ['"p5"']
clobberlist.append('"cc"')
for i in range(0, 32):
    clobberlist.append(F'"z{i}"')

d = {
'debug': _DEBUG,
'C': '',
'I': '',
'A': '',
'asmsuffix': '.d',      # double precision by default
'cfloat': 'float64_t',
'registers': 0,
'load': 0,
'store': 0,
'move': 0,
'movprfx': 0,
'zero': 0,
'add': 0,
'sub': 0,
'mul': 0,
'mac': 0,
'permutes': 0,
'neg': 0,
'addTimesI': 0,
'subTimesI': 0,
'timesI': 0,
'timesMinusI': 0,
'flops': 0,
'factor': 1,                    # multiplicity
'asmtableptr': 'x30',
'asmfetchbaseptr': 'x29',
'asmstorebaseptr': 'x28',
'asmdebugptr': 'r12',
'asminput': [],
'asmoutput': [],
'asmclobber': [],
'asmclobberlist': clobberlist,
'intrinfetchbase': '',
'intrinstorebase': '',
'cycles_LOAD_CHIMU': 0,
'cycles_PROJ': 0,
'cycles_PERM': 0,
'cycles_MULT_2SPIN': 0,
'cycles_RECON': 0,
'cycles_RESULT': 0,
'cycles_ZERO_PSI': 0,
'cycles_PREFETCH_L1': 0,
'cycles_PREFETCH_L2': 0
}

if PRECISION == 'single':
    d['asmsuffix'] = '.s'
    d['cfloat'] = 'float32_t'

# --------------------------------------------------------------------------------
# Grid
# --------------------------------------------------------------------------------

# Variables / Registers
result_00 = Register('result_00', asmreg='z0')
result_01 = Register('result_01', asmreg='z1')
result_02 = Register('result_02', asmreg='z2')
result_10 = Register('result_10', asmreg='z3')
result_11 = Register('result_11', asmreg='z4')
result_12 = Register('result_12', asmreg='z5')
result_20 = Register('result_20', asmreg='z6')
result_21 = Register('result_21', asmreg='z7')
result_22 = Register('result_22', asmreg='z8')
result_30 = Register('result_30', asmreg='z9')
result_31 = Register('result_31', asmreg='z10')
result_32 = Register('result_32', asmreg='z11')       # 12 Regs
Chi_00	  = Register('Chi_00', asmreg='z12')
Chi_01	  = Register('Chi_01', asmreg='z13')
Chi_02	  = Register('Chi_02', asmreg='z14')
Chi_10	  = Register('Chi_10', asmreg='z15')
Chi_11	  = Register('Chi_11', asmreg='z16')
Chi_12	  = Register('Chi_12', asmreg='z17')          # 6
UChi_00   = Register('UChi_00', asmreg='z18')
UChi_01   = Register('UChi_01', asmreg='z19')
UChi_02   = Register('UChi_02', asmreg='z20')
UChi_10   = Register('UChi_10', asmreg='z21')
UChi_11   = Register('UChi_11', asmreg='z22')
UChi_12   = Register('UChi_12', asmreg='z23')         # 6
U_00	  = Register('U_00', asmreg='z24')
U_10	  = Register('U_10', asmreg='z25')
U_20	  = Register('U_20', asmreg='z26')
U_01	  = Register('U_01', asmreg='z27')
U_11	  = Register('U_11', asmreg='z28')
U_21      = Register('U_21', asmreg='z29')            # 6    -> 30 Registers

table0    = Register('table0', asmreg='z30')
zero0     = Register('zero0', asmreg='z31')           # 2    -> 32 Registers
# can't overload temp1 / table due to type mismatch using intrinsics :(
# typecasting SVE intrinsics variables is not allowed

pg1       = Register('pg1', predication=True, asmreg='p5')
#pg2       = Register('pg2', predication=True, asmreg='p1')

# Overloaded with Chi_* and UChi_*
Chimu_00  = Register('Chimu_00', asmreg=Chi_00.asmreg)
Chimu_01  = Register('Chimu_01', asmreg=Chi_01.asmreg)
Chimu_02  = Register('Chimu_02', asmreg=Chi_02.asmreg)
Chimu_10  = Register('Chimu_10', asmreg=Chi_10.asmreg)
Chimu_11  = Register('Chimu_11', asmreg=Chi_11.asmreg)
Chimu_12  = Register('Chimu_12', asmreg=Chi_12.asmreg)
if ALTERNATIVE_REGISTER_MAPPING == False:
    Chimu_20  = Register('Chimu_20', asmreg=UChi_00.asmreg)
    Chimu_21  = Register('Chimu_21', asmreg=UChi_01.asmreg)
    Chimu_22  = Register('Chimu_22', asmreg=UChi_02.asmreg)
    Chimu_30  = Register('Chimu_30', asmreg=UChi_10.asmreg)
    Chimu_31  = Register('Chimu_31', asmreg=UChi_11.asmreg)
    Chimu_32  = Register('Chimu_32', asmreg=UChi_12.asmreg)        # 12 Registers
else: # wilson4.h
    Chimu_20  = Register('Chimu_20', asmreg=U_00.asmreg)
    Chimu_21  = Register('Chimu_21', asmreg=U_10.asmreg)
    Chimu_22  = Register('Chimu_22', asmreg=U_20.asmreg)
    Chimu_30  = Register('Chimu_30', asmreg=U_01.asmreg)
    Chimu_31  = Register('Chimu_31', asmreg=U_11.asmreg)
    Chimu_32  = Register('Chimu_32', asmreg=U_21.asmreg)

# debugging output
def debugall(msg=None, group='ALL'):
    global d
    if (d['debug'] == False):
        return
    write(F'std::cout << std::endl << "DEBUG -- {msg}" << std::endl; \\')
    if (group in ['ALL', 'result']):
        result_00.debug()
        result_01.debug()
        result_02.debug()
        result_10.debug()
        result_11.debug()
        result_12.debug()
        result_20.debug()
        result_21.debug()
        result_22.debug()
        result_30.debug()
        result_31.debug()
        result_32.debug()
    if (group in ['ALL', 'Chi']):
        Chi_00.debug()
        Chi_01.debug()
        Chi_02.debug()
        Chi_10.debug()
        Chi_11.debug()
        Chi_12.debug()
    if (group in ['ALL', 'UChi']):
        UChi_00.debug()
        UChi_01.debug()
        UChi_02.debug()
        UChi_10.debug()
        UChi_11.debug()
        UChi_12.debug()
    if (group in ['ALL', 'U']):
        U_00.debug()
        U_10.debug()
        U_20.debug()
        U_01.debug()
        U_11.debug()
        U_21.debug()
    if (group in ['ALL', 'Chimu']):
        Chimu_00.debug()
        Chimu_01.debug()
        Chimu_02.debug()
        Chimu_10.debug()
        Chimu_11.debug()
        Chimu_12.debug()
        Chimu_20.debug()
        Chimu_21.debug()
        Chimu_22.debug()
        Chimu_30.debug()
        Chimu_31.debug()
        Chimu_32.debug()

# --------------------------------------------------------------------------------
# Output
# --------------------------------------------------------------------------------

if ALTERNATIVE_LOADS == True:
    define(F'LOAD_CHIMU_0213_PLUG    LOAD_CHIMU_0213_{PRECSUFFIX}')
    define(F'LOAD_CHIMU_0312_PLUG    LOAD_CHIMU_0312_{PRECSUFFIX}')
    define(F'LOAD_CHIMU(x)')
else:
    #define(F'LOAD_CHIMU_{PRECSUFFIX}(x)           LOAD_CHIMU_INTERLEAVED_{PRECSUFFIX}(x)')
    define(F'LOAD_CHIMU(base)               LOAD_CHIMU_INTERLEAVED_{PRECSUFFIX}(base)')

if PREFETCH:
    define(F'PREFETCH_CHIMU_L1(A)           PREFETCH_CHIMU_L1_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PREFETCH_GAUGE_L1(A)           PREFETCH_GAUGE_L1_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PREFETCH_CHIMU_L2(A)           PREFETCH_CHIMU_L2_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PREFETCH_GAUGE_L2(A)           PREFETCH_GAUGE_L2_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PF_GAUGE(A)')
    define(F'PREFETCH_RESULT_L2_STORE(A)    PREFETCH_RESULT_L2_STORE_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PREFETCH_RESULT_L1_STORE(A)    PREFETCH_RESULT_L1_STORE_INTERNAL_{PRECSUFFIX}(A)')
    define(F'PREFETCH1_CHIMU(A)             PREFETCH_CHIMU_L1(A)')
#    define(F'PREFETCH1_CHIMU(A)')
    define(F'PREFETCH_CHIMU(A)              PREFETCH_CHIMU_L1(A)')
#    define(F'PREFETCH_CHIMU(A)')
else:
    define(F'PREFETCH_CHIMU_L1(A)')
    define(F'PREFETCH_GAUGE_L1(A)')
    define(F'PREFETCH_CHIMU_L2(A)')
    define(F'PREFETCH_GAUGE_L2(A)')
    define(F'PF_GAUGE(A)')
    define(F'PREFETCH1_CHIMU(A)')
    define(F'PREFETCH_CHIMU(A)')
    define(F'PREFETCH_RESULT_L2_STORE(A)')

# standard defines
define(F'LOCK_GAUGE(A)')
define(F'UNLOCK_GAUGE(A)')
define(F'MASK_REGS                      DECLARATIONS_{PRECSUFFIX}')
define(F'SAVE_RESULT(A,B)               RESULT_{PRECSUFFIX}(A); PREFETCH_RESULT_L2_STORE(B)')
define(F'MULT_2SPIN_1(Dir)              MULT_2SPIN_1_{PRECSUFFIX}(Dir)')
define(F'MULT_2SPIN_2                   MULT_2SPIN_2_{PRECSUFFIX}')
define(F'LOAD_CHI(base)                 LOAD_CHI_{PRECSUFFIX}(base)')
# don't need zero psi, everything is done in recons
#define(F'ZERO_PSI                       ZERO_PSI_{PRECSUFFIX}')
define(F'ADD_RESULT(base,basep)         LOAD_CHIMU(base); ADD_RESULT_INTERNAL_{PRECSUFFIX}; RESULT_{PRECSUFFIX}(base)')
# loads projections
define(F'XP_PROJ                        XP_PROJ_{PRECSUFFIX}')
define(F'YP_PROJ                        YP_PROJ_{PRECSUFFIX}')
define(F'ZP_PROJ                        ZP_PROJ_{PRECSUFFIX}')
define(F'TP_PROJ                        TP_PROJ_{PRECSUFFIX}')
define(F'XM_PROJ                        XM_PROJ_{PRECSUFFIX}')
define(F'YM_PROJ                        YM_PROJ_{PRECSUFFIX}')
define(F'ZM_PROJ                        ZM_PROJ_{PRECSUFFIX}')
define(F'TM_PROJ                        TM_PROJ_{PRECSUFFIX}')
# recons
define(F'XP_RECON                       XP_RECON_{PRECSUFFIX}')
define(F'XM_RECON                       XM_RECON_{PRECSUFFIX}')
define(F'XM_RECON_ACCUM                 XM_RECON_ACCUM_{PRECSUFFIX}')
define(F'YM_RECON_ACCUM                 YM_RECON_ACCUM_{PRECSUFFIX}')
define(F'ZM_RECON_ACCUM                 ZM_RECON_ACCUM_{PRECSUFFIX}')
define(F'TM_RECON_ACCUM                 TM_RECON_ACCUM_{PRECSUFFIX}')
define(F'XP_RECON_ACCUM                 XP_RECON_ACCUM_{PRECSUFFIX}')
define(F'YP_RECON_ACCUM                 YP_RECON_ACCUM_{PRECSUFFIX}')
define(F'ZP_RECON_ACCUM                 ZP_RECON_ACCUM_{PRECSUFFIX}')
define(F'TP_RECON_ACCUM                 TP_RECON_ACCUM_{PRECSUFFIX}')
# new permutes
define(F'PERMUTE_DIR0                   0')
define(F'PERMUTE_DIR1                   1')
define(F'PERMUTE_DIR2                   2')
define(F'PERMUTE_DIR3                   3')
define(F'PERMUTE                        PERMUTE_{PRECSUFFIX};')
# load table
#define(F'MAYBEPERM(A,perm)              if (perm) {{ A ; }}')
if PRECISION == 'double':
    define(F'LOAD_TABLE(Dir)                if (Dir == 0) {{ LOAD_TABLE0; }} else if (Dir == 1) {{ LOAD_TABLE1; }} else if (Dir == 2) {{ LOAD_TABLE2; }}')
    define(F'MAYBEPERM(Dir,perm)            if (Dir != 3) {{ if (perm) {{ PERMUTE; }} }}')
else:
    define(F'LOAD_TABLE(Dir)                if (Dir == 0) {{ LOAD_TABLE0; }} else if (Dir == 1) {{ LOAD_TABLE1 }} else if (Dir == 2) {{ LOAD_TABLE2; }} else if (Dir == 3) {{ LOAD_TABLE3; }}')
    define(F'MAYBEPERM(A,perm)              if (perm) {{ PERMUTE; }}')



write('// DECLARATIONS')
definemultiline(F'DECLARATIONS_{PRECSUFFIX}')
# debugging register
if d['debug'] == True:
    write('    Simd debugreg; \\')
# perm tables
if PRECISION == 'double':
    write('    const uint64_t lut[4][8] = { \\')
    write('        {4, 5, 6, 7, 0, 1, 2, 3}, \\')  #0 = swap register halves
    write('        {2, 3, 0, 1, 6, 7, 4, 5}, \\')  #1 = swap halves of halves
    write('        {1, 0, 3, 2, 5, 4, 7, 6}, \\')  #2 = swap re/im
    write('        {0, 1, 2, 4, 5, 6, 7, 8} };\\')  #3 = identity
else:
    write('    const uint32_t lut[4][16] = { \\')
    write('        {8, 9, 10, 11, 12, 13, 14, 15, 0, 1, 2, 3, 4, 5, 6, 7}, \\')  #0 = swap register halves
    write('        {4, 5, 6, 7, 0, 1, 2, 3, 12, 13, 14, 15, 8, 9, 10, 11}, \\')  #1 = swap halves of halves
    write('        {2, 3, 0, 1, 6, 7, 4, 5, 10, 11, 8, 9, 14, 15, 12, 13}, \\')  #2 = swap halves of halves of halves
    write('        {1, 0, 3, 2, 5, 4, 7, 6, 9, 8, 11, 10, 13, 12, 15, 14} }; \\') #3 = swap re/im

#newline(target='A')
result_00.declare()
result_01.declare()
result_02.declare()
result_10.declare()
result_11.declare()
result_12.declare()
result_20.declare()
result_21.declare()
result_22.declare()
result_30.declare()
result_31.declare()
result_32.declare()     # 12
Chi_00.declare()
Chi_01.declare()
Chi_02.declare()
Chi_10.declare()
Chi_11.declare()
Chi_12.declare()        # 6
UChi_00.declare()
UChi_01.declare()
UChi_02.declare()
UChi_10.declare()
UChi_11.declare()
UChi_12.declare()       # 6
U_00.declare()
U_10.declare()
U_20.declare()
U_01.declare()
U_11.declare()
U_21.declare()          # 6   -> 30 regs

# all predications true
pg1.declare()
if PRECISION == 'double':
    pg1.movestr('svptrue_b64()')
else:
    pg1.movestr('svptrue_b32()')

# tables
if PRECISION == 'double':
    write('    svuint64_t table0; \\', target='I')   #      -> 31 regs
else:
    write('    svuint32_t table0; \\', target='I')   #      -> 31 regs

zero0.declare()

# zero register
asmopen()
zero0.zero(zeroreg=True)
asmclose()
newline()

define('Chimu_00 Chi_00', target='I')
define('Chimu_01 Chi_01', target='I')
define('Chimu_02 Chi_02', target='I')
define('Chimu_10 Chi_10', target='I')
define('Chimu_11 Chi_11', target='I')
define('Chimu_12 Chi_12', target='I')
if ALTERNATIVE_REGISTER_MAPPING == False:
    define('Chimu_20 UChi_00', target='I')
    define('Chimu_21 UChi_01', target='I')
    define('Chimu_22 UChi_02', target='I')
    define('Chimu_30 UChi_10', target='I')
    define('Chimu_31 UChi_11', target='I')
    define('Chimu_32 UChi_12', target='I')
else: # wilson4.h
    define('Chimu_20 U_00', target='I')
    define('Chimu_21 U_10', target='I')
    define('Chimu_22 U_20', target='I')
    define('Chimu_30 U_01', target='I')
    define('Chimu_31 U_11', target='I')
    define('Chimu_32 U_21', target='I')
newline()


d['cycles_RESULT'] += 12
write('// RESULT')
definemultiline(F'RESULT_{PRECSUFFIX}(base)')
if ASM_STORE:
    curlyopen()
    #write('    SiteSpinor & ref(out[ss]); \\')
    asmopen()
    #pg1.loadpredication()
    #store_base_ptr("&ref[0][0]")
    #store_base_ptr(F"&ref[{STORE_BASE_PTR_COLOR_OFFSET}][0]")
    store_base_ptr(F"base + {STORE_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='I')
    store_base_ptr(F"base + {STORE_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='A')
    result_00.store("ref[0][0]")
    result_01.store("ref[0][1]")
    result_02.store("ref[0][2]")
    result_10.store("ref[1][0]")
    result_11.store("ref[1][1]")
    result_12.store("ref[1][2]")
    result_20.store("ref[2][0]")
    result_21.store("ref[2][1]")
    result_22.store("ref[2][2]")
    result_30.store("ref[3][0]")
    result_31.store("ref[3][1]")
    result_32.store("ref[3][2]")
    asmclose()
    debugall('RESULT', group='result')
    curlyclose()
newline()

# prefetch spinors from memory into L2 cache
d['factor'] = 0
d['cycles_PREFETCH_L2'] += 0 * d['factor']
write('// PREFETCH_CHIMU_L2 (prefetch to L2)')
definemultiline(F'PREFETCH_CHIMU_L2_INTERNAL_{PRECSUFFIX}(base)')
curlyopen()
fetch_base_ptr(F"base")
asmopen()
#pg1.loadpredication()
#fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
fetch_base_ptr(F"base", target='A')
prefetch_L2(F"base", 0)
prefetch_L2(F"base", 1)
prefetch_L2(F"base", 2)
asmclose()
curlyclose()
newline()

# prefetch spinors from memory into L1 cache
d['factor'] = 0
d['cycles_PREFETCH_L1'] += 0 * d['factor']
write('// PREFETCH_CHIMU_L1 (prefetch to L1)')
definemultiline(F'PREFETCH_CHIMU_L1_INTERNAL_{PRECSUFFIX}(base)')
curlyopen()
fetch_base_ptr(F"base")
asmopen()
#pg1.loadpredication()
fetch_base_ptr(F"base", target='A')
prefetch_L1(F"base", 0)
prefetch_L1(F"base", 1)
prefetch_L1(F"base", 2)
asmclose()
curlyclose()
newline()

# prefetch gauge from memory into L2 cache
d['factor'] = 0
d['cycles_PREFETCH_L2'] += 0 * d['factor']
write('// PREFETCH_GAUGE_L2 (prefetch to L2)')
definemultiline(F'PREFETCH_GAUGE_L2_INTERNAL_{PRECSUFFIX}(A)')
curlyopen()
if GRIDBENCH:   # referencing differs in Grid and GridBench
    write('    const auto & ref(U[sUn][A]); uint64_t baseU = (uint64_t)&ref + 3 * 3 * 64; \\')
else:
    write('    const auto & ref(U[sUn](A)); uint64_t baseU = (uint64_t)&ref + 3 * 3 * 64; \\')
asmopen()
#pg1.loadpredication()
#fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
fetch_base_ptr(F"baseU", target='A')
prefetch_L2(F"baseU", -1)
prefetch_L2(F"baseU", 0)
prefetch_L2(F"baseU", 1)
prefetch_L2(F"baseU", 2)
prefetch_L2(F"baseU", 3)
prefetch_L2(F"baseU", 4)
prefetch_L2(F"baseU", 5)
prefetch_L2(F"baseU", 6)
prefetch_L2(F"baseU", 7)
#prefetch_L2(F"baseU", 8)
asmclose()
curlyclose()
newline()

# prefetch gauge from memory into L1 cache
d['factor'] = 0
d['cycles_PREFETCH_L1'] += 0 * d['factor']
write('// PREFETCH_GAUGE_L1 (prefetch to L1)')
definemultiline(F'PREFETCH_GAUGE_L1_INTERNAL_{PRECSUFFIX}(A)')
curlyopen()
if GRIDBENCH:   # referencing differs in Grid and GridBench
    write('    const auto & ref(U[sU][A]); uint64_t baseU = (uint64_t)&ref; \\')
else:
    write('    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \\')
asmopen()
#pg1.loadpredication()
#fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
fetch_base_ptr(F"baseU", target='A')
prefetch_L1(F"baseU", 0)
prefetch_L1(F"baseU", 1)
prefetch_L1(F"baseU", 2)
asmclose()
curlyclose()
newline()

d['factor'] = 0
write('// LOAD_CHI')
definemultiline(F'LOAD_CHI_{PRECSUFFIX}(base)')
if ASM_LOAD_CHIMU:
    curlyopen()
    #write('    const SiteSpinor & ref(in[offset]); \\')
    asmopen()
    #fetch_base_ptr(F"base + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='I')
    #fetch_base_ptr(F"base + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='A')
    fetch_base_ptr(F"base", target='I')
    fetch_base_ptr(F"base", target='A')

    Chi_00.load("ref[0][0]", offset=0)
    Chi_01.load("ref[0][1]", offset=0)
    Chi_02.load("ref[0][2]", offset=0)
    Chi_10.load("ref[1][0]", offset=0)
    Chi_11.load("ref[1][1]", offset=0)
    Chi_12.load("ref[1][2]", offset=0)
    asmclose()
    debugall('LOAD_CHI', group='Chi')
    curlyclose()
newline()



d['factor'] = 8
# 12 loads = 12 issues, load latency = 8+1 cycles
# (not perfectly clear to me from docs)
d['cycles_LOAD_CHIMU'] += 11 * d['factor']
write('// LOAD_CHIMU')
definemultiline(F'LOAD_CHIMU_INTERLEAVED_{PRECSUFFIX}(base)')
if ASM_LOAD_CHIMU:
    curlyopen()
    #write('    const SiteSpinor & ref(in[offset]); \\')
    asmopen()
    pg1.loadpredication()
    #fetch_base_ptr("&ref[0][0]")
    #fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
    fetch_base_ptr(F"base + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='I')
    fetch_base_ptr(F"base + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='A')
    # Chimu_00.load("ref[0][0]")
    # Chimu_01.load("ref[0][1]")
    # Chimu_02.load("ref[0][2]")
    # Chimu_10.load("ref[1][0]")
    # Chimu_11.load("ref[1][1]")
    # Chimu_12.load("ref[1][2]")
    # Chimu_20.load("ref[2][0]")
    # Chimu_21.load("ref[2][1]")
    # Chimu_22.load("ref[2][2]")
    # Chimu_30.load("ref[3][0]")
    # Chimu_31.load("ref[3][1]")
    # Chimu_32.load("ref[3][2]")

    Chimu_00.load("ref[0][0]")  # minimum penalty for all directions
    Chimu_30.load("ref[3][0]")
    Chimu_10.load("ref[1][0]")
    Chimu_20.load("ref[2][0]")

    Chimu_01.load("ref[0][1]")
    Chimu_31.load("ref[3][1]")
    Chimu_11.load("ref[1][1]")
    Chimu_21.load("ref[2][1]")

    Chimu_02.load("ref[0][2]")
    Chimu_32.load("ref[3][2]")
    Chimu_12.load("ref[1][2]")
    Chimu_22.load("ref[2][2]")
    asmclose()
    debugall('LOAD_CHIMU', group='Chimu')
    curlyclose()
newline()

# alternative load chimu: dirac order 0213
# placed into asm (...)
d['factor'] = 0
d['cycles_LOAD_CHIMU'] += 11 * d['factor']
write('// LOAD_CHIMU_0213')
definemultiline(F'LOAD_CHIMU_0213_{PRECSUFFIX}')
if ASM_LOAD_CHIMU:
    curlyopen()
    write('    const SiteSpinor & ref(in[offset]); \\')
    asmopen()
    pg1.loadpredication()
    fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
    Chimu_00.load("ref[0][0]")  # reordered
    Chimu_20.load("ref[2][0]")

    Chimu_01.load("ref[0][1]")
    Chimu_21.load("ref[2][1]")

    Chimu_02.load("ref[0][2]")
    Chimu_22.load("ref[2][2]")

    Chimu_10.load("ref[1][0]")
    Chimu_30.load("ref[3][0]")

    Chimu_11.load("ref[1][1]")
    Chimu_31.load("ref[3][1]")

    Chimu_12.load("ref[1][2]")
    Chimu_32.load("ref[3][2]")
    asmclose()
    debugall('LOAD_CHIMU_0213', group='Chimu')
    curlyclose()
newline()

# alternative load chimu: dirac order 0312
# placed into asm (...)
d['factor'] = 0
d['cycles_LOAD_CHIMU'] += 11 * d['factor']
write('// LOAD_CHIMU_0312')
definemultiline(F'LOAD_CHIMU_0312_{PRECSUFFIX}')
if ASM_LOAD_CHIMU:
    curlyopen()
    write('    const SiteSpinor & ref(in[offset]); \\')
    asmopen()
    pg1.loadpredication()
    fetch_base_ptr(F"&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]")
    Chimu_00.load("ref[0][0]")  # reordered
    Chimu_30.load("ref[3][0]")

    Chimu_01.load("ref[0][1]")
    Chimu_31.load("ref[3][1]")

    Chimu_02.load("ref[0][2]")
    Chimu_32.load("ref[3][2]")

    Chimu_10.load("ref[1][0]")
    Chimu_20.load("ref[2][0]")

    Chimu_11.load("ref[1][1]")
    Chimu_21.load("ref[2][1]")

    Chimu_12.load("ref[1][2]")
    Chimu_22.load("ref[2][2]")
    asmclose()
    debugall('LOAD_CHIMU_0312', group='Chimu')
    curlyclose()
newline()

d['factor'] = 2
d['cycles_PERM'] += 1 * d['factor']
write('// LOAD_TABLE0')
definemultiline(F'LOAD_TABLE0')
asmopen()
table0.loadtable(0)
asmclose()
newline()

d['factor'] = 2
d['cycles_PERM'] += 1 * d['factor']
write('// LOAD_TABLE1')
definemultiline(F'LOAD_TABLE1')
asmopen()
table0.loadtable(1)
asmclose()
newline()

d['factor'] = 2
d['cycles_PERM'] += 1 * d['factor']
write('// LOAD_TABLE2')
definemultiline(F'LOAD_TABLE2')
asmopen()
table0.loadtable(2)
asmclose()
newline()

d['factor'] = 0
d['cycles_PERM'] += 1 * d['factor']
write('// LOAD_TABLE3')
definemultiline(F'LOAD_TABLE3')
asmopen()
table0.loadtable(3)
asmclose()
newline()

d['factor'] = 2     # factor is 2
d['cycles_PERM'] += 6 * d['factor']
write('// PERMUTE')
definemultiline(F'PERMUTE_{PRECSUFFIX}')
debugall('PERM PRE', group='Chi')
asmopen()
#table0.loadtable(2)
Chi_00.permute(2, table0)
Chi_01.permute(2, table0)
Chi_02.permute(2, table0)
Chi_10.permute(2, table0)
Chi_11.permute(2, table0)
Chi_12.permute(2, table0)
asmclose()
debugall('PERM POST', group='Chi')
newline()

write('// LOAD_GAUGE')
definemultiline(F'LOAD_GAUGE')
if GRIDBENCH:   # referencing differs in Grid and GridBench
    write('    const auto & ref(U[sU][A]); uint64_t baseU = (uint64_t)&ref; \\')
else:
    write('    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \\')
curlyopen()
asmopen()
pg1.loadpredication()
fetch_base_ptr(F"baseU + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='I')
if ASM_LOAD_GAUGE:
    fetch_base_ptr(F"baseU + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='A')
    U_00.load("ref[0][0]")
    U_10.load("ref[1][0]")
    U_20.load("ref[2][0]")
    U_01.load("ref[0][1]")
    U_11.load("ref[1][1]")
    U_21.load("ref[2][1]")
asmclose()
curlyclose()
newline()

d['factor'] = 8    # MULT_2SPIN executes 1 time per direction = 8 times total
# assume all U loads are hidden
# FCMLA issue latency = 2 cycles
# measurement: latency = 16 cycles if FULLY pipelined !?
# spec says 6+6+9 cycles
# 6 rounds of FCMLA, each with 6 FCMLA -> 21 - 6*2 = 9
d['cycles_MULT_2SPIN'] += 6 * 21 * d['factor']
write('// MULT_2SPIN')
definemultiline(F'MULT_2SPIN_1_{PRECSUFFIX}(A)')
curlyopen()
#write('    const auto & ref(U[sU][A]); \\')
if GRIDBENCH:   # referencing differs in Grid and GridBench
    write('    const auto & ref(U[sU][A]); uint64_t baseU = (uint64_t)&ref; \\')
else:
    write('    const auto & ref(U[sU](A)); uint64_t baseU = (uint64_t)&ref; \\')
asmopen()
#pg1.loadpredication()
#fetch_base_ptr("&ref[0][0]")
fetch_base_ptr(F"baseU + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='I')
fetch_base_ptr(F"baseU + {FETCH_BASE_PTR_COLOR_OFFSET} * 3 * 64", target='A')
#fetch_base_ptr(F"(uint64_t)&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]", target='I')
#fetch_base_ptr(F"(uint64_t)&ref[{FETCH_BASE_PTR_COLOR_OFFSET}][0]", target='A')
#fetch_base_ptr(F"&ref[0][{FETCH_BASE_PTR_COLOR_OFFSET}]")
if ASM_LOAD_GAUGE:
    U_00.load("ref[0][0]")
    U_10.load("ref[1][0]")
    U_20.load("ref[2][0]")
    U_01.load("ref[0][1]")
    U_11.load("ref[1][1]")
    U_21.load("ref[2][1]")

if MOVPRFX == False:
    UChi_00.zero()  # implementation specific
    UChi_10.zero()
    UChi_01.zero()
    UChi_11.zero()
    UChi_02.zero()
    UChi_12.zero()

    # round 1
    UChi_00.mul0(U_00, Chi_00) # FCMLA latency is 6+6+9 cycles
    UChi_10.mul0(U_00, Chi_10)
    UChi_01.mul0(U_10, Chi_00)
    UChi_11.mul0(U_10, Chi_10)
    UChi_02.mul0(U_20, Chi_00)
    UChi_12.mul0(U_20, Chi_10)
else:
    # round 1
    UChi_00.mul0(zero0, U_00, Chi_00, constructive=True) # FCMLA latency is 6+6+9 cycles
    UChi_10.mul0(zero0, U_00, Chi_10, constructive=True)
    UChi_01.mul0(zero0, U_10, Chi_00, constructive=True)
    UChi_11.mul0(zero0, U_10, Chi_10, constructive=True)
    UChi_02.mul0(zero0, U_20, Chi_00, constructive=True)
    UChi_12.mul0(zero0, U_20, Chi_10, constructive=True)

# round 2
UChi_00.mul1(U_00, Chi_00)
UChi_10.mul1(U_00, Chi_10)
UChi_01.mul1(U_10, Chi_00)
UChi_11.mul1(U_10, Chi_10)
UChi_02.mul1(U_20, Chi_00)
UChi_12.mul1(U_20, Chi_10)  # Chi_00 and Chi_10 available from here

if ASM_LOAD_GAUGE:
    U_00.load("ref[0][2]")      # U_00, U_10, U_20 overloaded
    U_10.load("ref[1][2]")      # early load
    U_20.load("ref[2][2]")      # A -->
asmclose()
debugall('MULT_2SPIN_1', group='UChi')
curlyclose()
newline()

write('// MULT_2SPIN_BACKEND')
definemultiline(F'MULT_2SPIN_2_{PRECSUFFIX}')
curlyopen()
asmopen()
# round 3
UChi_00.mac0(U_01, Chi_01)  # armclang separates fcmla(..., 0) and
UChi_10.mac0(U_01, Chi_11)  #                    fcmla(..., 90)
UChi_01.mac0(U_11, Chi_01)  # autonomously using intrinsics
UChi_11.mac0(U_11, Chi_11)
UChi_02.mac0(U_21, Chi_01)
UChi_12.mac0(U_21, Chi_11)
# round 4
UChi_00.mac1(U_01, Chi_01)
UChi_10.mac1(U_01, Chi_11)
UChi_01.mac1(U_11, Chi_01)
UChi_11.mac1(U_11, Chi_11)
UChi_02.mac1(U_21, Chi_01)
UChi_12.mac1(U_21, Chi_11)
# round 5
UChi_00.mac0(U_00, Chi_02)  # <-- A
UChi_10.mac0(U_00, Chi_12)
UChi_01.mac0(U_10, Chi_02)
UChi_11.mac0(U_10, Chi_12)
UChi_02.mac0(U_20, Chi_02)
UChi_12.mac0(U_20, Chi_12)
# round 6
UChi_00.mac1(U_00, Chi_02)
UChi_10.mac1(U_00, Chi_12)
UChi_01.mac1(U_10, Chi_02)
UChi_11.mac1(U_10, Chi_12)
UChi_02.mac1(U_20, Chi_02)
UChi_12.mac1(U_20, Chi_12)
asmclose()
debugall('MULT_2SPIN_2', group='UChi')
curlyclose()
newline()


#//      hspin(0)=fspin(0)+timesI(fspin(3));
#//      hspin(1)=fspin(1)+timesI(fspin(2));
d['factor'] = 1
# FCADD issue latency = 1, latency is 6+9
d['cycles_PROJ'] += 15 * d['factor']
write('// XP_PROJ')
definemultiline(F'XP_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0312_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.addTimesI(Chimu_00, Chimu_30)
Chi_01.addTimesI(Chimu_01, Chimu_31)
Chi_02.addTimesI(Chimu_02, Chimu_32)
Chi_10.addTimesI(Chimu_10, Chimu_20)
Chi_11.addTimesI(Chimu_11, Chimu_21)
Chi_12.addTimesI(Chimu_12, Chimu_22)
asmclose()
debugall('XP_PROJ', group='Chi')
curlyclose()
newline()

#//      fspin(0)=hspin(0);
#//      fspin(1)=hspin(1);
#//      fspin(2)=timesMinusI(hspin(1));
#//      fspin(3)=timesMinusI(hspin(0));
# does not occur in GridBench
d['factor'] = 0
d['cycles_RECON'] += 15 * d['factor']
write('// XP_RECON')
definemultiline(F'XP_RECON_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
if MOVPRFX == False:
    result_20.zero()
    result_21.zero()
    result_22.zero()
    result_30.zero()
    result_31.zero()
    result_32.zero()

    result_20.subTimesI(UChi_10)
    result_21.subTimesI(UChi_11)
    result_22.subTimesI(UChi_12)
    result_30.subTimesI(UChi_00)
    result_31.subTimesI(UChi_01)
    result_32.subTimesI(UChi_02)
else:
    result_20.subTimesI(zero0, UChi_10, constructive=True)
    result_21.subTimesI(zero0, UChi_11, constructive=True)
    result_22.subTimesI(zero0, UChi_12, constructive=True)
    result_30.subTimesI(zero0, UChi_00, constructive=True)
    result_31.subTimesI(zero0, UChi_01, constructive=True)
    result_32.subTimesI(zero0, UChi_02, constructive=True)

result_00.move(UChi_00) # don't reorder !
result_01.move(UChi_01)
result_02.move(UChi_02)
result_10.move(UChi_10)
result_11.move(UChi_11)
result_12.move(UChi_12)

# result_00.add(UChi_00)   # faster than move?
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
asmclose()
debugall('XP_RECON', group='result')
newline()


d['factor'] = 1
# FCADD issue latency = 1, latency is 6+9
d['cycles_RECON'] += 15 * d['factor']
write('// XP_RECON_ACCUM')
definemultiline(F'XP_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_20.subTimesI(UChi_10)
# result_21.subTimesI(UChi_11)
# result_22.subTimesI(UChi_12)
# result_30.subTimesI(UChi_00)
# result_31.subTimesI(UChi_01)
# result_32.subTimesI(UChi_02)
#
# result_00.add(UChi_00) # reordered
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)

result_30.subTimesI(UChi_00)    # reordered
result_00.add(UChi_00)

result_31.subTimesI(UChi_01)
result_01.add(UChi_01)

result_32.subTimesI(UChi_02)
result_02.add(UChi_02)

result_20.subTimesI(UChi_10)
result_10.add(UChi_10)

result_21.subTimesI(UChi_11)
result_11.add(UChi_11)

result_22.subTimesI(UChi_12)
result_12.add(UChi_12)
asmclose()
debugall('XP_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
# add/sub issue latency = 1, latency is 9
d['cycles_PROJ'] += 9 * d['factor']
write('// YP_PROJ')
definemultiline(F'YP_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0312_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.sub(Chimu_00, Chimu_30)
Chi_01.sub(Chimu_01, Chimu_31)
Chi_02.sub(Chimu_02, Chimu_32)
Chi_10.add(Chimu_10, Chimu_20)
Chi_11.add(Chimu_11, Chimu_21)
Chi_12.add(Chimu_12, Chimu_22)
asmclose()
debugall('YP_PROJ', group='Chi')
curlyclose()
newline()

d['factor'] = 1
# FCADD issue latency = 1, latency is 6+9
d['cycles_PROJ'] += 15 * d['factor']
write('// ZP_PROJ')
definemultiline(F'ZP_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0213_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.addTimesI(Chimu_00, Chimu_20)
Chi_01.addTimesI(Chimu_01, Chimu_21)
Chi_02.addTimesI(Chimu_02, Chimu_22)
Chi_10.subTimesI(Chimu_10, Chimu_30)
Chi_11.subTimesI(Chimu_11, Chimu_31)
Chi_12.subTimesI(Chimu_12, Chimu_32)
asmclose()
debugall('ZP_PROJ', group='Chi')
curlyclose()
newline()

d['factor'] = 1
# add/sub issue latency = 1, latency is 9
d['cycles_PROJ'] += 9 * d['factor']
write('// TP_PROJ')
definemultiline(F'TP_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0213_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.add(Chimu_00, Chimu_20)
Chi_01.add(Chimu_01, Chimu_21)
Chi_02.add(Chimu_02, Chimu_22)
Chi_10.add(Chimu_10, Chimu_30)
Chi_11.add(Chimu_11, Chimu_31)
Chi_12.add(Chimu_12, Chimu_32)
asmclose()
debugall('TP_PROJ', group='Chi')
curlyclose()
newline()

#//      hspin(0)=fspin(0)-timesI(fspin(3));
#//      hspin(1)=fspin(1)-timesI(fspin(2));

d['factor'] = 1
# FCADD issue latency = 1, latency is 6+9
d['cycles_PROJ'] += 15 * d['factor']
write('// XM_PROJ')
definemultiline(F'XM_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0312_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.subTimesI(Chimu_00, Chimu_30)
Chi_01.subTimesI(Chimu_01, Chimu_31)
Chi_02.subTimesI(Chimu_02, Chimu_32)
Chi_10.subTimesI(Chimu_10, Chimu_20)
Chi_11.subTimesI(Chimu_11, Chimu_21)
Chi_12.subTimesI(Chimu_12, Chimu_22)
asmclose()
debugall('XM_PROJ sub', group='Chi')
curlyclose()
newline()

d['factor'] = 1
d['cycles_RECON'] += 15 * d['factor']
write('// XM_RECON')
definemultiline(F'XM_RECON_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()

# only necessary if not zeroed before
if MOVPRFX == False:
    result_20.zero()
    result_21.zero()
    result_22.zero()
    result_30.zero()
    result_31.zero()
    result_32.zero()

    result_20.addTimesI(UChi_10) # <--
    result_21.addTimesI(UChi_11)
    result_22.addTimesI(UChi_12)
    result_30.addTimesI(UChi_00)
    result_31.addTimesI(UChi_01)
    result_32.addTimesI(UChi_02)
else:
    result_20.addTimesI(zero0, UChi_10, constructive=True) # <--
    result_21.addTimesI(zero0, UChi_11, constructive=True)
    result_22.addTimesI(zero0, UChi_12, constructive=True)
    result_30.addTimesI(zero0, UChi_00, constructive=True)
    result_31.addTimesI(zero0, UChi_01, constructive=True)
    result_32.addTimesI(zero0, UChi_02, constructive=True)

result_00.move(UChi_00)
result_01.move(UChi_01)
result_02.move(UChi_02)
result_10.move(UChi_10)
result_11.move(UChi_11)
result_12.move(UChi_12)
asmclose()
debugall('XM_RECON result', group='result')
newline()

d['factor'] = 1
# add/sub issue latency = 1, latency is 9
d['cycles_PROJ'] += 9 * d['factor']
write('// YM_PROJ')
definemultiline(F'YM_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0312_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.add(Chimu_00, Chimu_30)
Chi_01.add(Chimu_01, Chimu_31)
Chi_02.add(Chimu_02, Chimu_32)
Chi_10.sub(Chimu_10, Chimu_20)
Chi_11.sub(Chimu_11, Chimu_21)
Chi_12.sub(Chimu_12, Chimu_22)
asmclose()
debugall('YM_PROJ', group='Chi')
curlyclose()
newline()

d['factor'] = 1
# FCADD issue latency = 1, latency is 6+9
d['cycles_PROJ'] += 15 * d['factor']
write('// ZM_PROJ')
definemultiline(F'ZM_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0213_PLUG \\')
curlyopen()
asmopen()
#pg1.loadpredication()
Chi_00.subTimesI(Chimu_00, Chimu_20)
Chi_01.subTimesI(Chimu_01, Chimu_21)
Chi_02.subTimesI(Chimu_02, Chimu_22)
Chi_10.addTimesI(Chimu_10, Chimu_30)
Chi_11.addTimesI(Chimu_11, Chimu_31)
Chi_12.addTimesI(Chimu_12, Chimu_32)
asmclose()
debugall('ZM_PROJ', group='Chi')
curlyclose()
newline()

d['factor'] = 1
# add/sub issue latency = 1, latency is 9
d['cycles_PROJ'] += 9 * d['factor']
write('// TM_PROJ')
definemultiline(F'TM_PROJ_{PRECSUFFIX}')
if ALTERNATIVE_LOADS == True:
    write('    LOAD_CHIMU_0213_PLUG \\')
curlyopen()
asmopen()
pg1.loadpredication()
Chi_00.sub(Chimu_00, Chimu_20)
Chi_01.sub(Chimu_01, Chimu_21)
Chi_02.sub(Chimu_02, Chimu_22)
Chi_10.sub(Chimu_10, Chimu_30)
Chi_11.sub(Chimu_11, Chimu_31)
Chi_12.sub(Chimu_12, Chimu_32)
asmclose()
debugall('TM_PROJ', group='Chi')
curlyclose()
newline()

# does not occur in GridBench
d['factor'] = 0
# add/sub issue latency = 1, latency is 9
d['cycles_RECON'] += 15 * d['factor']
write('// XM_RECON_ACCUM')
definemultiline(F'XM_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
# result_20.addTimesI(UChi_10)
# result_21.addTimesI(UChi_11)
# result_22.addTimesI(UChi_12)
# result_30.addTimesI(UChi_00)
# result_31.addTimesI(UChi_01)
# result_32.addTimesI(UChi_02)
#
# # result_00.move(UChi_00)
# # result_01.move(UChi_01)
# # result_02.move(UChi_02)
# # result_10.move(UChi_10)
# # result_11.move(UChi_11)
# # result_12.move(UChi_12)
#
# # faster than move ?
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)

result_30.addTimesI(UChi_00)    # reordered
result_31.addTimesI(UChi_01)
result_32.addTimesI(UChi_02)

result_20.addTimesI(UChi_10)
result_21.addTimesI(UChi_11)
result_22.addTimesI(UChi_12)

result_00.add(UChi_00)
result_01.add(UChi_01)
result_02.add(UChi_02)
result_10.add(UChi_10)
result_11.add(UChi_11)
result_12.add(UChi_12)
asmclose()
debugall('XM_RECON_ACCUM', group='result')
newline()



d['factor'] = 1
d['cycles_RECON'] += 9 * d['factor']
write('// YP_RECON_ACCUM')
definemultiline(F'YP_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
# result_20.add(UChi_10)
# result_21.add(UChi_11)
# result_22.add(UChi_12)
# result_30.sub(UChi_00)
# result_31.sub(UChi_01)
# result_32.sub(UChi_02)

result_00.add(UChi_00)  # reordered
result_30.sub(UChi_00)

result_01.add(UChi_01)
result_31.sub(UChi_01)

result_02.add(UChi_02)
result_32.sub(UChi_02)

result_10.add(UChi_10)
result_20.add(UChi_10)

result_11.add(UChi_11)
result_21.add(UChi_11)

result_12.add(UChi_12)
result_22.add(UChi_12)
asmclose()
debugall('YP_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
d['cycles_RECON'] += 9 * d['factor']
write('// YM_RECON_ACCUM')
definemultiline(F'YM_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
# result_20.sub(UChi_10)
# result_21.sub(UChi_11)
# result_22.sub(UChi_12)
# result_30.add(UChi_00)
# result_31.add(UChi_01)
# result_32.add(UChi_02)

result_00.add(UChi_00)  # reordered
result_30.add(UChi_00)

result_01.add(UChi_01)
result_31.add(UChi_01)

result_02.add(UChi_02)
result_32.add(UChi_02)

result_10.add(UChi_10)
result_20.sub(UChi_10)

result_11.add(UChi_11)
result_21.sub(UChi_11)

result_12.add(UChi_12)
result_22.sub(UChi_12)
asmclose()
debugall('YM_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
d['cycles_RECON'] += 15 * d['factor']
write('// ZP_RECON_ACCUM')
definemultiline(F'ZP_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_20.subTimesI(UChi_00)
# result_21.subTimesI(UChi_01)
# result_22.subTimesI(UChi_02)
# result_30.addTimesI(UChi_10)
# result_31.addTimesI(UChi_11)
# result_32.addTimesI(UChi_12)
#
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
result_20.subTimesI(UChi_00)    # reordered
result_00.add(UChi_00)

result_21.subTimesI(UChi_01)
result_01.add(UChi_01)

result_22.subTimesI(UChi_02)
result_02.add(UChi_02)

result_30.addTimesI(UChi_10)
result_10.add(UChi_10)

result_31.addTimesI(UChi_11)
result_11.add(UChi_11)

result_32.addTimesI(UChi_12)
result_12.add(UChi_12)
asmclose()
debugall('ZP_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
d['cycles_RECON'] += 15 * d['factor']
write('// ZM_RECON_ACCUM')
definemultiline(F'ZM_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_20.addTimesI(UChi_00)
# result_21.addTimesI(UChi_01)
# result_22.addTimesI(UChi_02)
# result_30.subTimesI(UChi_10)
# result_31.subTimesI(UChi_11)
# result_32.subTimesI(UChi_12)
#
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
result_20.addTimesI(UChi_00)    # reordered
result_00.add(UChi_00)

result_21.addTimesI(UChi_01)
result_01.add(UChi_01)

result_22.addTimesI(UChi_02)
result_02.add(UChi_02)

result_30.subTimesI(UChi_10)
result_10.add(UChi_10)

result_31.subTimesI(UChi_11)
result_11.add(UChi_11)

result_32.subTimesI(UChi_12)
result_12.add(UChi_12)
asmclose()
debugall('ZM_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
d['cycles_RECON'] += 9 * d['factor']
write('// TP_RECON_ACCUM')
definemultiline(F'TP_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
# result_20.add(UChi_00)
# result_21.add(UChi_01)
# result_22.add(UChi_02)
# result_30.add(UChi_10)
# result_31.add(UChi_11)
# result_32.add(UChi_12)

result_00.add(UChi_00)  # reordered
result_20.add(UChi_00)

result_01.add(UChi_01)
result_21.add(UChi_01)

result_02.add(UChi_02)
result_22.add(UChi_02)

result_10.add(UChi_10)
result_30.add(UChi_10)

result_11.add(UChi_11)
result_31.add(UChi_11)

result_12.add(UChi_12)
result_32.add(UChi_12)
asmclose()
debugall('TP_RECON_ACCUM', group='result')
newline()

d['factor'] = 1
d['cycles_RECON'] += 9 * d['factor']
write('// TM_RECON_ACCUM')
definemultiline(F'TM_RECON_ACCUM_{PRECSUFFIX}')
asmopen()
#pg1.loadpredication()
# result_00.add(UChi_00)
# result_01.add(UChi_01)
# result_02.add(UChi_02)
# result_10.add(UChi_10)
# result_11.add(UChi_11)
# result_12.add(UChi_12)
# result_20.sub(UChi_00)
# result_21.sub(UChi_01)
# result_22.sub(UChi_02)
# result_30.sub(UChi_10)
# result_31.sub(UChi_11)
# result_32.sub(UChi_12)

result_00.add(UChi_00)  # reordered
result_20.sub(UChi_00)

result_01.add(UChi_01)
result_21.sub(UChi_01)

result_02.add(UChi_02)
result_22.sub(UChi_02)

result_10.add(UChi_10)
result_30.sub(UChi_10)

result_11.add(UChi_11)
result_31.sub(UChi_11)

result_12.add(UChi_12)
result_32.sub(UChi_12)
asmclose()
debugall('TM_RECON_ACCUM', group='result')
newline()

d['factor'] = 0
# have 12 instructions
# picking dual issue versions
d['cycles_ZERO_PSI'] += 6 * d['factor']
write('// ZERO_PSI')
definemultiline(F'ZERO_PSI_{PRECSUFFIX}')
asmopen()
pg1.loadpredication()
result_00.zero()
result_01.zero()
result_02.zero()
result_10.zero()
result_11.zero()
result_12.zero()
result_20.zero()
result_21.zero()
result_22.zero()
result_30.zero()
result_31.zero()
result_32.zero()
asmclose()
#debugall('ZERO_PSI', group='result')
newline()

# prefetch store spinors to L2 cache
d['factor'] = 0
d['cycles_PREFETCH_L2'] += 0 * d['factor']
write('// PREFETCH_RESULT_L2_STORE (prefetch store to L2)')
definemultiline(F'PREFETCH_RESULT_L2_STORE_INTERNAL_{PRECSUFFIX}(base)')
curlyopen()
fetch_base_ptr(F"base")
asmopen()
fetch_base_ptr(F"base", target='A')
prefetch_L2_store(F"base", 0)
prefetch_L2_store(F"base", 1)
prefetch_L2_store(F"base", 2)
asmclose()
curlyclose()
newline()

# prefetch store spinors to L1 cache
d['factor'] = 0
d['cycles_PREFETCH_L1'] += 0 * d['factor']
write('// PREFETCH_RESULT_L1_STORE (prefetch store to L1)')
definemultiline(F'PREFETCH_RESULT_L1_STORE_INTERNAL_{PRECSUFFIX}(base)')
curlyopen()
fetch_base_ptr(F"base")
asmopen()
fetch_base_ptr(F"base", target='A')
prefetch_L1_store(F"base", 0)
prefetch_L1_store(F"base", 1)
prefetch_L1_store(F"base", 2)
asmclose()
curlyclose()
newline()


d['factor'] = 0
write('// ADD_RESULT_INTERNAL')
definemultiline(F'ADD_RESULT_INTERNAL_{PRECSUFFIX}')
asmopen()
result_00.add(Chimu_00)
result_01.add(Chimu_01)
result_02.add(Chimu_02)
result_10.add(Chimu_10)
result_11.add(Chimu_11)
result_12.add(Chimu_12)
result_20.add(Chimu_20)
result_21.add(Chimu_21)
result_22.add(Chimu_22)
result_30.add(Chimu_30)
result_31.add(Chimu_31)
result_32.add(Chimu_32)
asmclose()
#debugall('ZERO_PSI', group='result')
newline()

# --------------------------------------------------------------------------------

# C
f = open('w.h', 'w')
f.write(d['C'])
f.close()

# intrin
f = open('wi.h', 'w')
f.write(d['I'])
f.close()

filename = ''
if PRECISION == 'double':
    filename = "Fujitsu_A64FX_intrin_double.h"
else:
    filename = "Fujitsu_A64FX_intrin_single.h"
f = open(filename, 'w')
f.write(LEGAL.format(filename))
f.write(d['I'])
f.close()


# asm
f = open('wa.h', 'w')
f.write(d['A'])
f.close()

filename = ''
if PRECISION == 'double':
    filename = "Fujitsu_A64FX_asm_double.h"
else:
    filename = "Fujitsu_A64FX_asm_single.h"
f = open(filename, 'w')
f.write(LEGAL.format(filename))
f.write(d['A'])
f.close()


# arithmetics instruction count, mul/mac = 2 instructions each
d['acount'] = d['add'] + d['sub'] + \
    d['mul'] + d['mac'] + d['addTimesI'] + d['subTimesI']

# permutations
d['permutes'] += 2*d['timesI'] + 1*d['timesMinusI']
d['neg'] = 1*d['timesI'] + 1*d['timesMinusI']

# instruction count,  mul/mac = 2 instructions each, +/- *i = 3 instructions each
d['icount'] = d['load'] + d['store'] + d['move'] + d['add'] + d['sub'] + \
    d['mul'] + d['mac'] + d['permutes'] + d['neg'] + \
    d['addTimesI'] + d['subTimesI'] + d['zero'] + d['movprfx']

# flops
d['flops'] = 4*d['mac'] + 3*d['mul'] + d['add'] + d['sub'] + \
    d['addTimesI'] + d['subTimesI']





print('Statistics')
print('')
print('Type                     Occurences      Total / Arith instructions')
print('-------------------------------------------------------------------')
print('Variables                {:4d}'.format(d['registers']))
print('')
print('load                     {:4d}'.format(d['load']))
print('store                    {:4d}'.format(d['store']))
print('move                     {:4d}'.format(d['move']))
print('movprfx                  {:4d}'.format(d['movprfx']))
print('zero                     {:4d}'.format(d['zero']))
print('negate                   {:4d}'.format(d['neg']))


print('add                      {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['add'], d['add'] /   d['icount'], d['add'] /   d['acount']))
print('sub                      {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['sub'], d['sub'] /   d['icount'], d['sub'] /   d['acount']))
print('mul                      {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['mul'], 2*d['mul'] / d['icount'], 2*d['mul'] /   d['acount']))
print('mac                      {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['mac'], 2*d['mac'] / d['icount'], 2*d['mac'] /   d['acount']))
print('addTimesI                {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['addTimesI'], 2*d['addTimesI'] / d['icount'], 2*d['addTimesI'] / d['acount']))
print('subTimesI                {:4d}              {:0.2f} / {:0.2f}'.\
    format(d['subTimesI'], 2*d['subTimesI'] / d['icount'], 2*d['subTimesI'] / d['acount']))

print('timesI                   {:4d}'.format(d['timesI']))
print('timesMinusI              {:4d}'.format(d['timesMinusI']))
print('permutes                 {:4d}              {:0.2f}'.\
    format(d['permutes'], d['permutes'] / d['icount']))
print('')
print('flops                    {:4d}'.format(d['flops']))
print('instruction count        {:4d}'.format(d['icount']))
print('arith. instruction count {:4d}              {:0.2f}'.\
    format(d['acount'], d['acount'] / d['icount']))


# ---- static pipeline resources consumption ----
FLA = 0
FLA += 2 * d['mac'] + 2 * d['mul']
FLA += 1 * d['addTimesI'] + 1 * d['subTimesI']
FLA += 1 * d['move']
FLA += 1 * d['permutes']
FLA += 1 * d['store']
FLA += 1 * d['zero']

FLB = 0
FLB += 1 * d['addTimesI'] + 1 * d['subTimesI']

FLAB = 0
FLAB += 1 * d['mac'] + 1 * d['mul']
FLAB += 1 * d['add'] + 1 * d['sub']
FLAB += 1 * d['neg'] + 1 * d['movprfx']
#FLAB += 1 * d['zero']


FL_slots = 2 * d['icount']
FL_micro_ops = FLA + FLB + FLAB

print('')
print('------------------------------------------------------------------')
print('')
print('Static FL slot usage')
print('')
print('  FLA                      {:4d}'.format(FLA))
print('  FLB                      {:4d}'.format(FLB))
print('  FLA/B                    {:4d}'.format(FLAB))

print('')
print('Static FL slot efficiency')
print('')
print('  Total FL slots           {:4d}'.format(FL_slots))
print('  FL slots occupied        {:4d}'.format(FL_micro_ops))
print('  FL slot efficiency       {:0.2f}'.format(FL_micro_ops / FL_slots))

cycles_total = d['cycles_ZERO_PSI'] + d['cycles_LOAD_CHIMU'] + \
    d['cycles_PROJ'] + d['cycles_PERM'] + d['cycles_MULT_2SPIN'] + \
    d['cycles_RECON'] + d['cycles_RESULT']
cycles_total_hidden = d['cycles_ZERO_PSI'] + \
    d['cycles_PROJ'] + d['cycles_MULT_2SPIN'] + \
    d['cycles_RECON']

# ---- dynamic estimate ----

print('')
print('Dynamic cycles estimate (incl. latencies)')
print('')
print('  ZERO_PSI                 {:4d}'.format(d['cycles_ZERO_PSI']))
print('  LOAD_CHIMU               {:4d}'.format(d['cycles_LOAD_CHIMU']))
print('  PROJ                     {:4d}'.format(d['cycles_PROJ']))
print('  PERM                     {:4d}'.format(d['cycles_PERM']))
print('  MULT_2SPIN               {:4d}'.format(d['cycles_MULT_2SPIN']))
print('  RECON                    {:4d}'.format(d['cycles_RECON']))
print('  STORE                    {:4d}'.format(d['cycles_RESULT']))
print('')
print('  Sum                      {:4d}'.format(cycles_total))
print('')
print('  Sum*                     {:4d}'.format(cycles_total_hidden))
print('  Total FL slots*          {:4d}'.format(cycles_total_hidden * 2))
print('  FL slots occupied*       {:4d}'.format(FL_micro_ops))
print('  FL slot efficiency*      {:0.2f}'.format(FL_micro_ops / (2*cycles_total_hidden)))
print('')
print('  *load/store/PERM hidden')

estimated_cycles = cycles_total_hidden
# Estimate percent peak DP; dual issue, fma
pp = 100 * 4 * d['flops'] / (2*2*8*estimated_cycles)
print('')
print('Model prediction')
print('')
print('  Cycles*                  {:4d}'.format(estimated_cycles))
print('  Percent peak*            {:4.1f} %'.format(pp))

# estimated RF throughput in GB/s @ 2.2 GHz
tp10 = (d['load'] + d['store']) * 64 * 2.2 / estimated_cycles
tp2  = (d['load'] + d['store']) * 64 * 1000.**3 * 2.2 / 1024.**3 / estimated_cycles
print('')
print('  Estimated RF throughput* {:4.1f}      GB/s'.\
    format(tp10))
print('  Estimated RF throughput* {:4.1f}      GiB/s'.\
    format(tp2))

# ---- dynamic pipeline resources consumption ----

runtime = measured_cycles  # runtime in cycles
pp_runtime = 100 * 4 * d['flops'] / (2*2*8*runtime)
runtime_FL_slots = 2 * runtime
delta = runtime - estimated_cycles


print('')
print('------------------------------------------------------------------')
print('')
print('Dynamic runtime analysis (cycles from measurements)')
print('')
print('  Cycles                   {:4d}'.format(runtime))
print('  Percent peak             {:4.1f} %'.format(pp_runtime))
print('  Deviation from estimate  {:4d}               {:4.2f} %'.\
    format(delta, 100. * abs(delta/runtime)))
print('  Deviation per direction  {:4.1f}'.format(delta/8))

# estimated RF throughput in GB/s @ 2.2 GHz
tp10_rt = (d['load'] + d['store']) * 64 * 2.2 / runtime
tp2_rt  = (d['load'] + d['store']) * 64 * 1000.**3 * 2.2 / 1024.**3 / runtime
print('')
print('  RF throughput            {:4.1f}      GB/s'.\
    format(tp10_rt))
print('  RF throughput            {:4.1f}      GiB/s'.\
    format(tp2_rt))
print('')
print('  Total FL slots           {:4d}'.format(runtime_FL_slots))
print('  FL slots occupied        {:4d}'.format(FL_micro_ops))
print('  FL slot efficiency       {:0.2f}'.format(FL_micro_ops / runtime_FL_slots))
print('')
