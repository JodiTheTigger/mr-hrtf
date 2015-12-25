#!/usr/bin/python3
# mr-hrtf. A simple 3d sound head related transfer function (HRTF) filter
# Copyright (C) 2015 Richard Maxwell <jodi.the.tigger@gmail.com>
# This file is part of mr-hrtf
# mr-hrtf is free software: you can redistribute it and/or modify
# it under the terms of the GNU Affero General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Affero General Public License for more details.
# You should have received a copy of the GNU Affero General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>

import scipy.io as io
import os
import textwrap
import numpy

from struct import *

cipic_root = os.path.expanduser('.')

anthro = io.loadmat(cipic_root + '/anthropometry/anthro.mat')

subject_root = cipic_root + '/standard_hrir_database'
subjects = os.listdir(subject_root)

anthros = 17 + (2 * (8 + 2))

with open('cipic.hrim', 'bw+') as meta_data:

    size_in_bytes = 4 + (4 * len(anthro['id']) * (1 + anthros))

    meta_data.write(pack('bbbb', ord('H'),ord('R'),ord('I'),ord('M')))
    meta_data.write(pack('<L', size_in_bytes))
    meta_data.write(pack('<L', 0))
    # TIFF style header, hrim version 0.

    for i in range(0, len(anthro['id'])):
        meta_data.write(pack('<L', *anthro['id'][i]))
        meta_data.write(pack('<17f', *anthro['X'][i]))
        meta_data.write(pack('<8f', *anthro['D'][i][:8]))
        meta_data.write(pack('<2f', *anthro['theta'][i][:2]))
        meta_data.write(pack('<8f', *anthro['D'][i][8:]))
        meta_data.write(pack('<2f', *anthro['theta'][i][2:]))

for hrir_dir in subjects:

    if (hrir_dir.startswith('subject_')):

        print(hrir_dir)

        subject_id = int(hrir_dir[-3:])

        full_path = subject_root + '/' + hrir_dir + '/hrir_final.mat'

        hrir = io.loadmat(full_path)

        i = 0
        count = 0;
        for s in anthro['id']:
            if s == subject_id:
                i = count
            count = count + 1

        if i == len(anthro['id']):
            print("can't find subject " + subject_id + " in anthro.mat")
            sys.exit(1)

        with open(hrir_dir + '.hrir', 'bw+') as data:

            size_in_bytes = 4 + (200 * 25 * 50 * 2 * 4) + (anthros * 4);

            data.write(pack('bbbb', ord('H'),ord('R'),ord('I'),ord('R')))
            data.write(pack('<L', size_in_bytes))
            data.write(pack('<L', 0))
            # TIFF style header, hrir pack version 0.

            data.write(pack('<17f', *anthro['X'][i]))
            data.write(pack('<8f', *anthro['D'][i][:8]))
            data.write(pack('<2f', *anthro['theta'][i][:2]))
            data.write(pack('<8f', *anthro['D'][i][8:]))
            data.write(pack('<2f', *anthro['theta'][i][2:]))
            # Anthro data, x, D and theta

            for az in range(0, 25):
                for e in range(0, 50):
                    data.write(pack('<200f', *hrir['hrir_l'][az][e]))
                    data.write(pack('<200f', *hrir['hrir_r'][az][e]))
                    # hrir data, interleaved one left, one right.
